#!/usr/bin/env python3
import os, sys
import os.path as path
import pysam
import argparse
import collections
import multiprocessing

SELFPATH = path.dirname(__file__)

parser_ar = argparse.ArgumentParser(prog='get_TMB_length.py',
                                    usage='get_TMB_length [OPTION] <bam_file> <bed_file> ',
                                    description='Show the coverage statistics for a bam file',
                                    epilog='',
                                    formatter_class=argparse.RawTextHelpFormatter)

parser_ar.add_argument('bam_file', help='FILE. A sorted and indexed bam file or a bam file list (one file each line)', metavar='bam_file')
parser_ar.add_argument('bed_file', help='FILE. a BED file', metavar='bed_file')

parser_ar.add_argument('-o', metavar='"length_results.txt"', default='length_results.txt', help='FILE. Result file', dest='output_file')
parser_ar.add_argument('-c', metavar='[30]', default=30, type=int, help='INT. Coverage cutoff for the enriched regions.', dest='coverage_cutoff')
parser_ar.add_argument('-t', metavar='[2]', default=2, type=int, help='INT. The processes ', dest='thread')
parser_ar.add_argument('-v', action='store_true', default=False, help='Print out counting message ?', dest='VERBOSE')

paramters = parser_ar.parse_args()

BAM_FILE = paramters.bam_file
BED_FILE = paramters.bed_file

OUTPUT_FILE = paramters.output_file
COVERAGE_CUTOFF = paramters.coverage_cutoff
CHUNK_SIZE = paramters.thread
VERBOSE = paramters.VERBOSE

def get_length(bamfile: str, bedfile: str, coverage_cutoff: int):
    '''
    This function does the coverage statistics for the regions in a bed file. Return a tuple:
    return region_length_int, effective_bases_int, coverage_length_lst

    Parameters:
        **bamfile**: string
            The bam file (must be in absolute path) to be calculated. We assume it has been coordinate-sorted and indexed

        **bedfile**: string
            The regions in this bed file will be checked. Only first three columns will be used.

        **cutoff_tu...**: arbitrary argumants, must be integers
            Contains the coverage cutoff, such as 50, 100, 200.

    Returns:
        **tuple**:
            Return a tuple, (region_length_int, effective_bases_int, coverage_length_lst)
    '''
    bamfile_af = pysam.AlignmentFile(path.realpath(path.expanduser(bamfile)), 'rb')
    bam_references_tu = bamfile_af.references
    bamfile_af.close()

    region_length_int = 0
    loci = collections.namedtuple('loci', ['chrom', 'pos'])  # loci(str chrom, int pos)
    amplicon = collections.namedtuple('amplicon', ['start', 'end', 'info'])   # amplicon(loci start, loci end, info=all strings after third columns)
    amplicons_lst = []
    # read all BED lines into an amplicons list
    with open(path.realpath(path.expanduser(bedfile)), 'rt') as bedfile_f:
        i = 0
        for line_str in bedfile_f.readlines():
            i += 1
            try:
                if line_str.strip() != '':
                    line_lst = line_str.strip().split('\t')
                    if len(line_lst) <= 2:
                        raise ValueError
                else:
                    raise ValueError
            except:
                message = 'Line {i}, {line_str} has wrong format for BED files.'.format(i=i, line_str=line_str)
                print(message)
                continue

            try:
                # sometime the sequence names in banm are 'chr1', 'chr2'.... but in BED file are '1', '2'...
                if line_lst[0] in bam_references_tu:
                    chrom_str = line_lst[0]
                elif 'chr' + line_lst[0] in bam_references_tu:
                    chrom_str = 'chr' + line_lst[0]
                elif line_lst[0][3:] in bam_references_tu:
                    chrom_str = line_lst[0][3:]
                else:
                    raise ValueError

                if len(line_lst) == 3:
                    amplicons_lst.append(amplicon(start=loci(chrom=chrom_str, pos=int(line_lst[1])), end=loci(chrom=chrom_str, pos=int(line_lst[2])), info=''))
                elif len(line_lst) > 3:
                    amplicons_lst.append(amplicon(start=loci(chrom=chrom_str, pos=int(line_lst[1])), end=loci(chrom=chrom_str, pos=int(line_lst[2])), info='\t'.join(line_lst[3:])))
                else:
                    pass
                region_length_int = region_length_int + int(line_lst[2]) - int(line_lst[1]) + 1
            except IndexError:
                message = 'BED file format error. BED file must have at least three columns. Skip\n{}'.format(line_str)
                print(message)
                continue
            except ValueError:
                message = 'BED file format error. The second and third column must be integers. Skip\n{}'.format(line_str)
                print(message)
                continue
            except:
                message = 'An error occur when parse BED file. Skip\n{}'.format(line_str)
                print(message)
                continue

    if len(amplicons_lst) == 0:
        message = 'BED file contains no amplicon region. Check the BED file.'
        raise ValueError(message)

    # calculate length where coverage >= coverage_cutoff
    coverage_length_int = 0
    bamfile_af = pysam.AlignmentFile(bamfile, 'rb')
    j = 0
    for amp in amplicons_lst:
        for pileupcolumn in bamfile_af.pileup(contig=amp.start.chrom, start=amp.start.pos, stop=amp.end.pos, truncate=True, max_depth=9999999, stepper='nofilter', ignore_overlaps=False, min_base_quality=0, ignore_orphans=False):
            effective_coverage_int = 0
            for pileupread in pileupcolumn.pileups:
                segment = pileupread.alignment
                if segment.is_duplicate or segment.is_qcfail or segment.is_secondary or segment.is_supplementary or segment.is_unmapped:
                    pass
                else:
                    effective_coverage_int += 1

            # print(pileupcolumn.reference_name, pileupcolumn.reference_pos, effective_coverage_int)
            if effective_coverage_int >= coverage_cutoff:
                coverage_length_int += 1

            j += 1
            if j % 100 == 0 and VERBOSE:
                print(j, flush=True)

    bamfile_af.close()

    return coverage_length_int


def __wrapper(bamfile: str, bedfile: str, coverage_cutoff: int, output: str, lock):

    bamfile_str = path.realpath(path.expanduser(bamfile))
    bedfile_str = path.realpath(path.expanduser(bedfile))
    print('Counting length {0} ...'.format(path.split(bamfile)[1]), flush=True)
    length_int = 0
    length_int = get_length(bamfile_str, bedfile_str, COVERAGE_CUTOFF)

    message = '{0}\t{1}'.format(path.split(bamfile)[1], length_int)
    print(message, flush=True)
    lock.acquire()
    with open(path.realpath(path.expanduser(output)), 'at', buffering=1) as output_f:
        output_f.writelines(message + '\n')
    lock.release()
    return 0

def main(argvList = sys.argv, argv_int = len(sys.argv)):

    bamfile_lst = []
    try:
        bamfile_af = pysam.AlignmentFile(path.realpath(path.expanduser(BAM_FILE)), 'rb')
        bamfile_lst.append(path.realpath(path.expanduser(BAM_FILE)))
        bamfile_af.close()
    except ValueError:  # maybe this is a list
        with open(path.realpath(path.expanduser(BAM_FILE)), 'rt') as f:
            for line_str in f.readlines():
                if line_str.strip()[0] == '#':
                    continue
                bamfile_lst.append(path.realpath(path.expanduser(line_str.strip())))

        if len(bamfile_lst) == 0:
            message = ''
            sys.exit(message)

    #output_f = open(path.realpath(path.expanduser(OUTPUT_FILE)), 'wt', buffering=1)
    #header = 'File\tLength\n'
    #output_f.writelines(header)
    #output_f.close()

    lock = multiprocessing.Manager().Lock()
    arguments_lst = []
    for bamfile in bamfile_lst:
        # __wrapper(bamfile, BED_FILE, COVERAGE_CUTOFF, OUTPUT_FILE, lock)
        arguments_lst.append((bamfile, BED_FILE, COVERAGE_CUTOFF, OUTPUT_FILE, lock))

    p = multiprocessing.Pool(processes=CHUNK_SIZE)
    p.starmap_async(__wrapper, arguments_lst)
    p.close()
    p.join()

    print('Done')
    return

if __name__ == '__main__':
    main()

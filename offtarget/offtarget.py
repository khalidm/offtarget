#!/usr/bin/env python

'''
A program to measure on/off target behaviour of a given set of
Hi-Plex primers.

Requies two fastq files, a bam file and primer coordinates as input.

Primers are inserted into a Trie.

Each read in the fastq files is processed by trying to find a matching
primer in the Trie.

Each read in the BAM file is associated with its primer (if possible)
and then checked to see if it aligns in the correct location.

There are a few situations to consider:

- the read is unmapped (thus we do not know its coordinates, but
  we _might_ know its primers)
- the read is mapped:
     - we could not find any primers matching the read
     - we could only find one primer matching the read
     - we found two primers matching the read:
             - the two primers do not belong together
             - the two primers belong together:
                       - the read does not sufficiently overlap the
                         expected location of the primers
                       - (***) the read overlaps the expected location
                         of the primers (the read is on-target)

We allow the user to specify prefix length which says how long a prefix of
a read (from the fastq file) we should test against the primers. Making this
short gives more chance for a match. Making this long makes it more specific.

Example usage:

offtarget --prefix 15 --primers primer_coords.tsv \
          --fastq1 R1.fastq --fastq2 R2.fastq --bam example.bam > example.out

'''

import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv
from argparse import ArgumentParser
import datrie
from collections import namedtuple, defaultdict
from Bio import SeqIO
import logging
import sys
import re
import os

DESCRIPTION = 'Plot offtarget reads from Hi-Plex sequencing data'
DEFAULT_OUTPUT_FILE = 'out.png'
DEFAULT_LOG_FILE = 'offtarget.log'

def start_log(log):
    '''Initiate program logging. If no log file is specified then
    log output goes to DEFAULT_LOGFILE.'''
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    # Log the command line that was used to run the program
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

def parse_args():
    parser = ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--primers', type=str, required=True,
        help='Name of input primer file, in CSV format: CHROM POS_FWD POS_REV NAME_FWD NAME_REV SEQ_FWD PRIM_FWD_START PRIM_FWD_END SCORE SEQ_REV'),
    parser.add_argument('--fastq1', type=str,
        help='Fastq file representing all read 1s in all read pairs')
    parser.add_argument('--fastq2', type=str,
        help='Fastq file representing all read 2s in all read pairs')
    parser.add_argument('--bam', type=str,
        help='Name of input BAM file')
    parser.add_argument('--prefix', type=int, required=False,
        help='Length of primer prefix to test against reads, defaults to unbounded')
    parser.add_argument(
        '--log', metavar='FILE', type=str, default=DEFAULT_LOG_FILE,
        help='Log progress in FILENAME.')
    return parser.parse_args()

BASE_PAIRS = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

def reverse_complement(seq):
    return ''.join([BASE_PAIRS.get(base, base) for base in reversed(seq)])

def reverse_seq(seq):
    return ''.join([base for base in reversed(seq)])

def complement(seq):
    return ''.join([BASE_PAIRS.get(base, base) for base in seq])

PrimerInfo = namedtuple('PrimerInfo', ['chrom', 'dir', 'name', 'pos', 'seq'])

def read_primers(args):
    trie = datrie.Trie(u'ATGCN')
    primer_names = []
    with open(args.primers) as primer_file:
        reader = csv.reader(primer_file, delimiter='\t')
        for row in reader:
            try:
                chrom, pos_fwd, pos_rev, name_fwd, name_rev, seq_fwd, prim_fws_start, prim_fwd_end, prim_score, seq_rev = row[:10]
                pos_fwd = int(pos_fwd)
                pos_rev = int(pos_rev)
            except ValueError:
                print("Can't parse row: {}".format(row))
            else:
                primer_names.extend([name_fwd, name_rev])
                if args.prefix is not None:
                    # use only a prefix of each primer in the trie
                    seq_fwd = unicode(seq_fwd.upper()[:args.prefix])
                    seq_rev = unicode(seq_rev.upper()[:args.prefix])
                else:
                    # use the whole primer in the trie
                    seq_fwd = unicode(seq_fwd.upper())
                    seq_rev = unicode(seq_rev.upper())
                trie[unicode(seq_fwd)] = PrimerInfo(chrom, 'fwd', name_fwd, pos_fwd, seq_fwd)
                trie[unicode(seq_rev)] = PrimerInfo(chrom, 'rev', name_rev, pos_rev, seq_rev)
    return trie, primer_names

class ReadPrimerMapper(object):
    def __init__(self, primers):
        # (read_id, read number) -> [PrimerInfo]
        self.read_id_to_primers = {}
        self.primers = primers

    def read_fastq(self, fastq_filename, read_number):
        with open(fastq_filename) as handle:
            num_reads = 0
            logging.info("Processing fastq file: {}".format(fastq_filename))
            for record in SeqIO.parse(handle, "fastq") :
                num_reads += 1
                read_id = record.id
                seq = unicode(str(record.seq).upper())
                matches = [primer_info for (primer_seq, primer_info) in self.primers.prefix_items(seq)]
                self.read_id_to_primers[(read_id, read_number)] = matches
            logging.info("Number of reads = {}".format(num_reads))

    def read_primer_stats(self):
        stats = {}
        for read_id, primers in self.read_id_to_primers.items():
            num_primers = len(primers)
            if num_primers in stats:
                stats[num_primers] += 1
            else:
                stats[num_primers] = 1
        return stats

    def show_read_primer_map(self):
        for read_id, primers in self.read_id_to_primers.items():
            print("{}: {}".format(read_id, ','.join([p.dir + ' ' + p.name for (_seq, p) in primers])))


def matching_primer_names(primer1, primer2):
    '''Primers are matching if they are on the same chrom and
    have complementary names such as:
    CHEK2_46_F1 CHEK2_46_R1
    '''
    regex = re.compile(r'(?P<region>.+_.+)_(?P<direction>[F|R])(?P<number>.+)')
    match1 = regex.match(primer1)
    match2 = regex.match(primer2)
    return match1.group is not None and match2.group is not None and \
           match1.group('region') == match2.group('region') and \
           (match1.group('direction') == 'F' and match2.group('direction') == 'R' or \
           match1.group('direction') == 'R' and match2.group('direction') == 'F') and \
           match1.group('number') == match2.group('number')

class PrimerCount(object):
    def __init__(self):
        # unmapped + mapped_offtarget = off_target
        self.on_target = 0
        self.off_target = 0
        self.unmapped = 0
        self.mapped_offtarget = 0

def process_bam(args, read_id_to_primers):
    bam = pysam.Samfile(args.bam)
    read_info = {}
    primer_counts = defaultdict(PrimerCount)
    reads_seen_primers = defaultdict(set)
    num_reads = 0
    num_mapped_reads = 0
    num_unmapped_reads = 0
    for read in bam:
        num_reads += 1
        read_id = read.qname
        read_number = None
        if read.is_read1:
            read_number = 1
        if read.is_read2:
            read_number = 2
        # do we have primers for this read?
        if (read_id, read_number) in read_id_to_primers:
            primers = read_id_to_primers[(read_id, read_number)]
            # check that the read is aligned to the reference
            if not read.is_unmapped:
                num_mapped_reads += 1
                read_start_pos = read.pos
                read_end_pos = read.aend - 1
                try:
                    chrom = bam.getrname(read.rname)
                except:
                    chrom = "?"

                # there should really only be one matching primer
                # but we are a bit more generous just in case the
                # test for a primer found multiple possible solutions
                for primer in primers:
                    overlap = read_start_pos <= primer.pos <= read_end_pos
                    if chrom == primer.chrom and overlap:
                        primer_counts[primer.name].on_target += 1
                    else:
                        # read is mapped somewhere but in the wrong location for this primer
                        primer_counts[primer.name].off_target += 1
                        primer_counts[primer.name].mapped_offtarget += 1
            else:
                # An unmapped read
                num_unmapped_reads += 1
                for primer in primers:
                    # An umapped read is technically offtarget but to
                    # avoid double counting we will not increment the
                    # off_target count and only increment the unmapped count
                    #primer_counts[primer.name].off_target += 1
                    primer_counts[primer.name].unmapped += 1
        else:
            primer_counts[primer.name].unmapped += 1

    logging.info("Number reads {}".format(num_reads))
    logging.info("Number mapped reads {}".format(num_mapped_reads))
    logging.info("Number unmapped reads {}".format(num_unmapped_reads))
    return primer_counts

def print_primer_counts(args, primer_names, primer_counts):
    print("primer name, on target, off target, mapped off target, unmapped")
    on_target_counts = []
    off_target_counts = []
    off_target_mapped_counts = []
    unmapped_counts = []
    for primer in primer_names:
        if primer in primer_counts:
            info = primer_counts[primer]
            line = [primer, str(info.on_target), str(info.off_target),
                    str(info.mapped_offtarget), str(info.unmapped)]
            print(','.join(line))
            on_target_counts.append(info.on_target)
            off_target_counts.append(info.off_target)
            off_target_mapped_counts.append(info.mapped_offtarget)
            unmapped_counts.append(info.unmapped)
        else:
            print('{},,,,'.format(primer))
            on_target_counts.append(0)
            off_target_counts.append(0)
            off_target_mapped_counts.append(0)
            unmapped_counts.append(0)
    plot_counts(args, 'On target', 'on_target', on_target_counts)
    plot_counts(args, 'Off target', 'off_target', off_target_counts)
    plot_counts(args, 'Off target mapped', 'on_target_mapped', off_target_mapped_counts)
    plot_counts(args, 'Unmapped', 'unmapped', unmapped_counts)

def plot_counts(args, title, extension, counts):
    bamfile_basename = os.path.basename(args.bam)
    plt.ylabel("read count")
    plt.xlabel("primer")
    plt.title('{} count for {}'.format(title, bamfile_basename))
    x_axis = range(len(counts))
    plt.bar(x_axis, counts)
    output_filename = bamfile_basename + '.' + extension + '.png'
    plt.savefig(output_filename)
    plt.close()

def main():
    args = parse_args()
    start_log(args.log)
    primer_trie, primer_names = read_primers(args)
    read_primer_mapper = ReadPrimerMapper(primer_trie)
    read_primer_mapper.read_fastq(args.fastq1, 1)
    read_primer_mapper.read_fastq(args.fastq2, 2)
    for num_primers, count in read_primer_mapper.read_primer_stats().items():
        logging.info("Reads with {} primers: {}".format(num_primers, count))
    primer_counts = process_bam(args, read_primer_mapper.read_id_to_primers)
    print_primer_counts(args, primer_names, primer_counts)

if __name__ == '__main__':
    main()

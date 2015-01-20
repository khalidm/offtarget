#!/usr/bin/env python

import pysam
import matplotlib.pyplot as plt
import csv
from argparse import ArgumentParser
import datrie
#from collections import namedtuple
#from operator import attrgetter

DESCRIPTION = 'Plot offtarget reads from Hi-Plex sequencing data'
DEFAULT_OUTPUT_FILE = 'out.png'

def parse_args():
    parser = ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--output', type=str, default=DEFAULT_OUTPUT_FILE,
        help='Name of output file for graph in PNG format, must end in .png'),
    parser.add_argument('--coords', type=str, required=True,
        help='Name of input primer coords, in CSV format'),
    parser.add_argument('--primers', type=str, required=True,
        help='Name of input primer file, in CSV format'),
    parser.add_argument('bam', type=str,
        help='Name of input BAM file')
    return parser.parse_args()

def read_coords(args):
    coord_info = {}
    with open(args.coords) as coord_file:
        reader = csv.reader(coord_file, delimiter='\t')
        for row in reader:
            chrom, fwd_coord, rev_coord, fwd_name, rev_name = row
            coord_info[fwd_name] = (chrom, int(fwd_coord))
            coord_info[rev_name] = (chrom, int(rev_coord))
    return coord_info 

def read_primers(args, coord_info):
    trie = datrie.Trie(u'ATGCN')
    with open(args.primers) as primer_file:
        reader = csv.reader(primer_file, delimiter='\t')
        for row in reader:
            name, primer = row[0:2]
            if name in coord_info:
                (chrom, pos) = coord_info[name]
                trie[unicode(primer)] = (name, chrom, pos)
            else:
                print("No coords for {}".format(name))
    return trie

BASE_PAIRS = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

def reverse_complement(seq):
    return ''.join([BASE_PAIRS.get(base, base) for base in reversed(seq)])

def reverse_seq(seq):
    return ''.join([base for base in reversed(seq)])

def complement(seq):
    return ''.join([BASE_PAIRS.get(base, base) for base in seq])


def process_reads(args, primer_trie):
    bam = pysam.Samfile(args.bam)
    read_info = {}
    for read in bam:
        bases = read.seq.upper()
        rev_complement_bases = reverse_complement(bases) 
        matches_norm = primer_trie.prefixes(unicode(bases))
        matches_rev_complement = primer_trie.prefixes(unicode(rev_complement_bases))
        try:
            chrom = bam.getrname(read.rname)
        except:
            chrom = "?"
        pos = read.pos
        qname = read.qname
        if qname in read_info:
            read_info[qname].append(bool(matches_norm))
        else:
            read_info[qname] = [bool(matches_norm)]

    match_types = {}

    for (qname, matches) in read_info.items():
        tuple_match = tuple(matches)
        if tuple_match in match_types:
            match_types[tuple_match] += 1
        else:
            match_types[tuple_match] = 1

    for (match_type, count) in match_types.items():
        print("{}: {}".format(match_type, count))

def main():
    args = parse_args()
    coord_info = read_coords(args)
    primer_trie = read_primers(args, coord_info)
    process_reads(args, primer_trie)

if __name__ == '__main__':
    main()

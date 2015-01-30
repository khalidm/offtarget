#!/usr/bin/env python

import csv
from argparse import ArgumentParser

DESCRIPTION = 'Merge block coordinates with primer sequences'

def parse_args():
    parser = ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--blocks', type=str, required=True,
        help='Name of input block coords, in CSV format'),
    parser.add_argument('--primers', type=str, required=True,
        help='Name of input primer sequence file, in CSV format'),
    return parser.parse_args()

def main():
    args = parse_args()

    with open(args.blocks) as block_coord_file, \
         open(args.primers) as primer_seq_file:
    
        dna_reader = csv.reader(primer_seq_file, delimiter='\t')
        coord_reader = csv.reader(block_coord_file, delimiter='\t')
        
        dna_dict = {}
        for row in dna_reader:
            primer_name, primer_seq = row[:2]
            dna_dict[primer_name] = primer_seq
        
        for row in coord_reader:
            fwd_seq = None
            rev_seq = None
            chrom, start, end, fwd_name = row
            if fwd_name in dna_dict:
                fwd_seq = dna_dict[fwd_name]
            fwd_name_parts = fwd_name.split('_')
            fwd_name_last_part = fwd_name_parts[-1]
            rev_name_last_part = 'R' + fwd_name_last_part[1:]
            rev_name = '_'.join(fwd_name_parts[:-1] + [rev_name_last_part])
            if rev_name in dna_dict:
                rev_seq = dna_dict[rev_name]

            if fwd_seq and rev_seq:
                print('\t'.join([chrom, start, end, fwd_name, rev_name, fwd_seq, rev_seq]))
            elif not fwd_seq:
                print("Could not find sequence for {}".format(fwd_name))
            else:
                print("Could not find sequence for {}".format(rev_name))


if __name__ == '__main__':
    main()

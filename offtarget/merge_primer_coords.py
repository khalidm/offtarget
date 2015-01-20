#!/usr/bin/env python

import csv

dna_filename = "745plex_names+seq.txt"
coord_filename = "250plex_primers.txt"

dna_file = open(dna_filename)
coord_file = open(coord_filename)

dna_reader = csv.reader(dna_file, delimiter='\t')
coord_reader = csv.reader(coord_file, delimiter='\t')

dna_dict = {}
for row in dna_reader:
    primer_name, primer_seq = row
    dna_dict[primer_name] = primer_seq

outfilename = "primer_coords.csv"
outfile = open(outfilename, "w")
writer = csv.writer(outfile)

for row in coord_reader:
    chrom, start, end, num, name = row
    if name in dna_dict:
        seq = dna_dict[name]
        writer.writerow([chrom, start, end, name, seq])
        print(len(seq))
    else:
        print("Could not find sequence for {}".format(name))

dna_file.close()
coord_file.close()
outfile.close()

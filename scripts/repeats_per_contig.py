#!/usr/bin/env python

import sys
from Bio import SeqIO

def calculate_repeat_percentage(fasta_file):
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = record.seq
            total_length = len(sequence)
            lower_count = sum(1 for base in sequence if base.islower())
            lower_percentage = (lower_count / total_length) * 100
            print(f"{record.id}: {lower_percentage:.2f}%")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        calculate_repeat_percentage(sys.argv[1])
    else:
        print("Please provide a filename as a command-line argument.")
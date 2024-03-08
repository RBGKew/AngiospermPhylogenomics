#!/usr/bin/env python
"""
# check_direction.py is a Python script that compares a sequence set to an alignment in 
# order to identify sequences which direction was changed
# 
# #### Version
# 
# - version = 1.0
# - date: 2022.04.08
# 
# #### Author
# 
# - Alexandre R. Zuntini
"""

import argparse
import os
import sys
import time

try:
    from Bio import SeqIO
except:
    print("\033[31m ERROR: Import library failed. \33[m")
    print("BioPython library needs to be installed.")
    sys.exit()

###
# Classes
###

class Sequence(list):
    """ Store sequence organised by samples and genes """
    def __init__(self, sample, seq):
        self.sample = sample
        self.seq = seq

###
# Functions
###

def read_seqs(gene, seq_dir, seq_ext, verbose=False):
    all_sequences = []
    seqs_file = os.path.join(seq_dir, gene + seq_ext)
    for seq_record in SeqIO.parse(seqs_file, "fasta"):
        seq = str(seq_record.seq).upper()
        all_sequences.append(Sequence(seq_record.name, seq))
    return(all_sequences)

def read_aligns(gene, align_dir, align_ext, verbose=False):
    all_sequences = []
    align_file = os.path.join(align_dir, gene + align_ext)
    for seq_record in SeqIO.parse(align_file, "fasta"):
        seq = str(seq_record.seq.reverse_complement()).replace("-", "").upper()
        all_sequences.append(Sequence(seq_record.name, seq))
    return(all_sequences)

def main():
    start = time.process_time()
    parser = argparse.ArgumentParser(description="Specify parameters for handling sequences.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Run in verbose mode")
    parser.add_argument("--seq_dir", dest="seq_dir", help="Gene to identify", default="01_seqs")
    parser.add_argument("--seq_ext", dest="seq_ext", help="Gene to identify", default=".fasta")
    parser.add_argument("--align_dir", dest="align_dir", help="Gene to identify", default="02_aligned_i1")
    parser.add_argument("--align_ext", dest="align_ext", help="Gene to identify", default="_i1_aligned.fasta")
    args = parser.parse_args()

    genes = [fn.replace(args.seq_ext, "") for fn in os.listdir(args.seq_dir) if fn.endswith(args.seq_ext)]
    error_file_destination = os.path.join(args.seq_dir, "fix")
    if not os.path.exists(error_file_destination):
        os.mkdir(error_file_destination)
    for gene in genes:
        try:
            all_seqs = read_seqs(gene, args.seq_dir, args.seq_ext, verbose = args.verbose)
        except FileNotFoundError:
            print(f"File {args.seq_dir}/{gene}{args.seq_ext} not found")
            continue
        all_seqs_list = [x.seq for x in all_seqs]
        try:
            all_aligns = read_aligns(gene, args.align_dir, args.align_ext, verbose = args.verbose)
        except FileNotFoundError:
            print(f"File {args.align_dir}/{gene}{args.align_ext} not found")
            continue
        inverted_seqs = [x.sample for x in all_aligns if x.seq in all_seqs_list]
        if inverted_seqs:
            with open (os.path.join(error_file_destination, gene + "_direction.txt"), "w") as f:
                [f.write(">_R_" + x + "\n") for x in inverted_seqs]

    print(time.process_time() - start)

if __name__ == "__main__":
    main()

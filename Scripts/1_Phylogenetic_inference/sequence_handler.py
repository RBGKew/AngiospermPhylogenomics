#!/usr/bin/env python
"""
# sequence_handler.py is a Python script to manipulate sequences, allowing to rename
# sequences, read and save by genes or samples or combined as well as removing specific 
# sequences (e.g. from TreeShrink) or reverse complement (following mafft output)
# 
# #### Version
# 
# - version = 1.0.3
# - date: 2020.05.21
# 
# #### Author
# 
# - Alexandre R. Zuntini
# - a.zuntini@kew.org

"""

import argparse
import os
import sys
import re
import time
from itertools import groupby
from collections import Counter
from pathlib import Path

try:
    from Bio import SeqIO
except ModuleNotFoundError:
    print("\033[31m ERROR: Import library failed. \33[m")
    print("BioPython library must be installed.")
    sys.exit()
try:
    import pandas as pd
except ModuleNotFoundError:
    print("\033[31m ERROR: Import library failed. \33[m")
    print("Pandas library needs to be installed to output sequence sizes.")
    sys.exit()

start = time.perf_counter()

###
# Variables 
###

primary_data_repo_path = "."
secondary_data_repo_path = "~/users_area/data_repo/angio353_local/"
min_length_default = 10


###
# Classes
###

class Sequence(list):
    """ Store sequence organised by samples and genes """
    def __init__(self, sample, gene, seq):
        self.sample = sample
        self.gene = gene
        self.seq = seq

    @classmethod
    def from_concat(cls, seq_name, seq):
        """ Alternative constructor for concatenated names as sample-gene """
        sample, gene = seq_name.split("-")
        return cls(sample, gene, seq)

    def prepend_names(self):
        """ Prepend samples and gene names with 's' and 'g' and pad sample name to width 5 """
        self.sample = "s" + str(self.sample).replace("s", "", 1).zfill(5)
        self.gene = self.gene if self.gene.startswith("g") else "g" + self.gene

    def export_seq(self, naming_variable, f):
        """ Export sequence named with selected instance variable such as self.sample, self.gene or self.name """
        f.write(">" + naming_variable + "\n" + str(self.seq) + "\n")

    def create_concat_name(self):
        """ Concatenate sample and gene with dash """
        self.name = "-".join([self.sample, self.gene])

    def len(self):
        """ Return sequence length """
        return(len(self.seq))

    def ungapped_len(self):
        """ Return sequence length excluding gaps """
        return(len(self.seq.ungap("-")))

    def reverse_complementary(self):
        """ Return the reverse complementary of the sequence """
        self.seq = self.seq.reverse_complement()

    def remove_gaps(self):
        self.seq = self.seq.ungap("-")

###
# Functions
###

def load_sequences(files_set, in_type, data_type="exons", verbose=False):
    all_sequences = []
    for infile in files_set:
        #sample_filename = sample + "_" + data_type + ".fasta"
        infile_path = os.path.join(primary_data_repo_path, infile)
        #sample_file = os.path.join(primary_data_repo_path, data_type , sample_filename)
        if verbose:
            print("Reading " + infile_path)
        if not os.path.isfile(infile_path):
            infile_path = os.path.join(secondary_data_repo_path, data_type , infile)
            if not os.path.isfile(infile_path):
                with open("missing_files.csv", "a") as f:
                    f.write(infile + "\n")
                continue

        with open(infile_path, "r") as file:
            first_line = file.readline()
        elements_first_line = first_line.split(" ")[0].split("-") # Remove sequence description and try to split
        if len(elements_first_line) == 1:
            file_no_prefix, *_ = os.path.basename(infile).split(".")#[0]
            if in_type == "sample":
                #print("Reading as sample")
                for seq_record in SeqIO.parse(infile_path, "fasta"):
                    all_sequences.append(Sequence(file_no_prefix, seq_record.name, seq_record.seq))
            else:
                #print("Reading as gene")
                for seq_record in SeqIO.parse(infile_path, "fasta"):
                    all_sequences.append(Sequence(seq_record.name, file_no_prefix, seq_record.seq))
                    #print(seq_record)
                    #print(all_sequences[-1].name)
        else:
            for seq_record in SeqIO.parse(infile_path, "fasta"):
                all_sequences.append(Sequence.from_concat(seq_record.name, seq_record.seq))
    return(all_sequences)

def parse_drop_list(file, verbose):
    """ Parse drop list from TreeShrink
        Expected format: gene:sample1 sample2 ... """
    genes = []
    if verbose:
        print(f"Parsing {file}")
    with open(file, "r") as f:
        drop_list = list(f.read().splitlines())
        seqs_to_drop = []
    for line in drop_list:
        gene, samples = line.split(":")
        if samples:
            genes.append(gene)
            sample_list = samples.split("\t")
            [seqs_to_drop.append("-".join([s,gene])) for s in sample_list if s]
    gene_set = set(genes)
    if verbose:
        print(f"{len(seqs_to_drop)} sequences to discard in {len(gene_set)} genes")
    return([gene_set, seqs_to_drop])

def get_seq_length(seqs, verbose):
    """ Tabulate sequence length by sample and gene """
    gene_list = [x.gene for x in seqs]
    gene_set = set(gene_list)
    if verbose:
        print(f"{len(gene_set)} genes found")
    sample_list = [x.sample for x in seqs]
    sample_set = set(sample_list)
    if verbose:
        print(f"{len(sample_set)} samples found")
    df_melt = pd.DataFrame(
        data={"sample": sample_list, "gene": gene_list, "length": [x.len() for x in seqs]})
    return(df_melt)

def is_fasta(file):
    """" Test if a file is fasta by starting with '>' """
    try:
        with open(file, "r") as file:
            first_line = file.readline()
    except OSError:
        print(f"File {file} not found")
        return False
    if first_line.startswith(">"):
        return True

def main():
    parser = argparse.ArgumentParser(description="Specify parameters for handling sequences")
    parser.add_argument("-v", "--verbose", action="store_true", help="Run in verbose mode")
    parser.add_argument("-f", "--files", help="List of input files, a single file, or a directory", default=[])
    parser.add_argument("-e", "--extension", help="List of input files or a directory", default="fasta")
    parser.add_argument("-s", "--source", help="Directory where sequence files are located", default=".")
    parser.add_argument("-d", "--destination", help="Directory to save the files.", default=".")
    parser.add_argument("--in_type", dest="in_type", help="List of input files", default="sample")
    parser.add_argument("--prepend", dest="prepend", action="store_true",
    					help="Prepend samples' and genes' names with 's' and 's' if absent.")
    parser.add_argument("-r", "--reverse", help="List of sequences to be reverse", default=[])
    parser.add_argument("--no_gap", dest="no_gap", action="store_true", help="Remove gap from sequences.")
    parser.add_argument("--drop_list", dest="drop_list", help="List of sequences to drop from TreeShrink", default=[])
    parser.add_argument("--dropped_only", dest="dropped_only", action="store_true",
                        help="Export only genes with dropped sequences. Takes precedence over --genes")
    parser.add_argument("--samples", dest="samples", help="List of samples to include", default=[])
    parser.add_argument("--genes", dest="genes", help="List of genes to include", default=[])
    parser.add_argument("--min_length", dest="min_length", help="Mininum lenght of a sequence to be kept",
    					default=min_length_default)
    parser.add_argument("--by_sample", dest="by_sample", action="store_true", help="Save sequences grouped by sample.")
    parser.add_argument("--by_gene", dest="by_gene", action="store_true", help="Save sequences grouped by gene.")
    parser.add_argument("--single_file", dest="single_file", action="store_true",
    					help="Save all sequences in a single file, named according to --out argument.")
    parser.add_argument("--relative_stats", dest="relative_stats", action="store_true",
                        help="Save relative sequences' length by gene.")
    parser.add_argument("-l", "--length", action="store_true", help="Save sequences' length by gene.")
    parser.add_argument("--summary", dest="summary", action="store_true", help="Save sequences' length summary.")
    parser.add_argument("-t", "--target", help="Target file used to calculate average recovery per gene", default=[])
    parser.add_argument("--save_sub_dir", dest="save_sub_dir", action="store_true", help="Save by gene.")
    parser.add_argument("-o", "--out", help="Prefix used to save stat files.", default="")
    parser.add_argument("--full_name", dest="full_name", action="store_true",
    					help="Export sequence name as sample-gene")
    parser.add_argument("--append", dest="append", action="store_true",
    					help="Append stat files instead of overwriting them.")
    args = parser.parse_args()

    try:
        min_length = int(args.min_length)
    except ValueError:
        print("min length argument must be an integer\n")
        print("Value set to {}".format(str(min_length_default)))
        min_length = min_length_default

    if not args.by_sample and not args.by_gene and not args.single_file and not args.length and not args.summary:
        args.single_file = True
        print("Enforcing export as single file")

    export_mode = "a" if args.append else "w"
    data_type = "exons"
    in_type = args.in_type if args.in_type in ["sample", "gene"] else "sample"
    destination = args.destination if args.destination else "."
    out_prefix = args.out + "_" if args.out != "" else ""
    in_extension  = args.extension

    if args.source:
        source_dir = args.source
    else:
        source_dir = "."

    if args.files:
        path_to_files = os.path.abspath(args.files)
    else:
        path_to_files = source_dir
        print(f"Reading all files in the directory: {source_dir}")

    if os.path.isfile(path_to_files):
        if is_fasta(path_to_files):
            file_paths = [path_to_files]
        else:
            with open(path_to_files, "r") as f:
                files_list = list(f.read().splitlines())
                file_paths = [os.path.join(source_dir, x) for x in files_list]
    elif os.path.isdir(path_to_files):
        file_paths = []
        for root, folders, files in os.walk(path_to_files):
            for file in files:
                if in_extension in file and not os.path.basename(file).startswith("."):
                    file_paths.append(os.path.join(root, file))
        print("Directory search done")
    else:
        print("File or directory not found. \nExiting...")
        sys.exit()

    files_set = set(file_paths)

    if args.verbose:
         print(files_set)

    Path(destination).mkdir(parents=True, exist_ok=True)

    all_sequences = load_sequences(files_set, in_type, data_type=data_type, verbose=args.verbose)

    # Append samples and gene with "s" and "g" if not present
    if args.prepend:
        [x.prepend_names() for x in all_sequences]

    # Create full name as sample-gene if to be used in exports
    # Provisional do later clear the TreeShrink results
    if args.full_name or args.single_file or args.reverse or args.drop_list:
        [x.create_concat_name() for x in all_sequences]

    if args.reverse:
        path_to_reverse_file = os.path.abspath(args.reverse)
        if os.path.isfile(path_to_reverse_file):
            with open(path_to_reverse_file, "r") as f:
                seqs_to_reverse = list(f.read().splitlines())
            [x.reverse_complementary() for x in all_sequences if x.name in seqs_to_reverse]
        else:
            print(f"List of sequences to revert ({path_to_reverse_file}) not found. \nExiting...")
            sys.exit()

    if args.drop_list:
        if os.path.isfile(args.drop_list):
            genes_in_drop_list, seqs_to_drop = parse_drop_list(os.path.abspath(args.drop_list), args.verbose)
            all_sequences = [x for x in all_sequences if x.name not in seqs_to_drop]
        else:
            print(f"List of genes to drop ({args.drop_list}) not found\nExiting...")
            sys.exit()

    if args.samples:
        all_sequences_orig = all_sequences
        samples_fn = os.path.abspath(args.samples)
        samples = []
        with open(samples_fn, "r") as f:
            samples = list(f.read().splitlines())
        sample_set = set(samples)
        original_samples = [x.sample for x in all_sequences_orig]
        original_sample_set = set(original_samples)
        all_sequences = [x for x in all_sequences_orig if x.sample in sample_set]
        if len(all_sequences) == 0:
            print(f"All original {len(original_sample_set)} samples were dropped. Please revise your selection!")
            print("Exiting...")
            sys.exit()
        new_samples = [x.sample for x in all_sequences]
        new_sample_set = set(new_samples)
        difference_sets = len(original_sample_set) - len(new_sample_set)
        if difference_sets == 0:
            print("All samples were kept.")
        else:
            print(
             f"{difference_sets} sample(s) dropped. Original dataset had {len(original_sample_set)} sample(s) and filtered dataset has {len(new_sample_set)}.")

    if args.dropped_only or args.genes:
        if args.dropped_only:
            if args.drop_list and genes_in_drop_list:
                gene_set = genes_in_drop_list
            elif not args.drop_list:
                print("Drop list not available.\nExiting...")
                sys.exit()
            else:
                print("Drop list empty.\nExiting...")
                with open("seq_output.txt", "w") as f:
                    genes = f.write("Empty drop list")
                sys.exit()
        else:
            genes_fn = os.path.abspath(args.genes)
            genes = []
            with open(genes_fn, "r") as f:
                genes = list(f.read().splitlines())
            gene_set = set(genes)
        all_sequences_orig = all_sequences
        original_genes = [x.gene for x in all_sequences_orig]
        original_gene_set = set(original_genes)
        all_sequences = [x for x in all_sequences_orig if x.gene in gene_set]
        if len(all_sequences) == 0:
            print(f"All original {len(original_gene_set)} samples were dropped. Please revise your selection!")
            print("Exiting...")
            sys.exit()
        new_genes = [x.gene for x in all_sequences]
        new_gene_set = set(new_genes)
        difference_sets = len(original_gene_set) - len(new_gene_set)
        if difference_sets == 0:
            print("All genes were kept.")
        else:
            print(
                f"{difference_sets} genes(s) dropped. Original dataset had {len(original_gene_set)} gene(s) and filtered dataset has {len(new_gene_set)}.")

    if args.no_gap:
        [x.remove_gaps() for x in all_sequences]

    if args.length or args.summary:
        print("Exporting sequence length")
        df_melt = get_seq_length(seqs=all_sequences, verbose=args.verbose)
        df = df_melt.pivot(index="sample", columns="gene").fillna(0)
        if args.target:
            path_to_target = [os.path.abspath(args.target)]
            target_sequences = load_sequences(path_to_target, in_type, data_type=data_type, verbose=args.verbose)
            [x.prepend_names() for x in target_sequences]
            target_df_melt = get_seq_length(seqs=target_sequences, verbose=args.verbose)
            test_seq = str(target_sequences[0].seq).upper()
            target_df = target_df_melt.pivot_table(index=None, columns="gene")
            if sum(map(test_seq.count, ("A", "C", "G", "T")))/len(test_seq) < .9:
                print("Target is aminoacid")
                target_df = target_df.multiply(3)
            #print(target_df)
            #target_df.to_csv("target_4_stats.csv", mode=export_mode)
        #df_melt = pd.DataFrame(
        #    data={"sample": sample_list, "gene": gene_list, "length": [x.len() for x in all_sequences]})
        #df_melt.pivot(index="sample", columns="gene").fillna(0)


        # start = time.process_time()
        # df = df_melt.pivot_table(index="sample", columns="gene", aggfunc="max", fill_value=0)
        # print(time.process_time() - start)

        # print(df.shape)

        if args.length:
            df.to_csv(out_prefix + "stats.csv", mode=export_mode)

        # NOT WORKING!!!
        if args.relative_stats:
            pass
            #if args.target:
            #    ref_df = target_df.round()
            #else:
            #    ref_df = df_melt.pivot_table(index=None, columns="gene").round()
            #print(ref_df.length)
#
            ##ref_dict = ref_df.to_dict()
            ##print(ref_dict)
            #print(df)
            #df_rel = df.div(ref_df.length)
            #print(df_rel)
            #df_rel.to_csv(out_prefix + "rel_stats.csv", mode=export_mode)

        if args.summary:
            n_genes = df[df > 0].count(1)
            total_length = df[df > 0].sum(1)
            summary_df = pd.concat([n_genes, total_length], axis=1)
            summary_df.to_csv(out_prefix + "summary.csv", mode=export_mode)

    if args.by_sample:
    	sample_list = [x.sample for x in all_sequences]
    	sample_set = set(sample_list)

    	for sample in sample_set:
    	    filename = sample + ".fasta"
    	    if args.save_sub_dir:
    	        sample_prefix = sample[0:3]
    	        Path(destination + "/" + sample_prefix).mkdir(parents=True, exist_ok=True)
    	        file_path = os.path.join(destination, sample_prefix, filename)
    	    else:
    	        file_path = os.path.join(destination, filename)

    	    with open(file_path, "w") as f:
    	        if args.full_name:
    	            [x.export_seq(x.name, f) for x in all_sequences if x.sample == sample and x.ungapped_len() >= min_length]
    	        else:
    	            [x.export_seq(x.gene, f) for x in all_sequences if x.sample == sample and x.ungapped_len() >= min_length]
    	    if args.verbose:
    	        print(file_path + " saved")

    if args.by_gene:
    	gene_list = [x.gene for x in all_sequences]
    	gene_count = Counter(gene_list)
    	gene_set = set(gene_list)

    	for gene in gene_set:
    	    filename = gene + ".fasta"
    	    file_path = os.path.join(destination, filename)
    	    with open(file_path, "w") as f:
    	        if args.full_name:
    	            [x.export_seq(x.name, f) for x in all_sequences if x.gene == gene and x.ungapped_len() >= min_length]
    	        else:
    	            [x.export_seq(x.sample, f) for x in all_sequences if x.gene == gene and x.ungapped_len() >= min_length]
    	    if args.verbose:
    	        print(filename + " exported")

    if args.single_file:
        basename = args.out if args.out else "seqs"
        filename = basename + ".fasta"
        file_path = os.path.join(destination, filename)
        with open(file_path, "w") as f:
            [x.export_seq(x.name, f) for x in all_sequences if x.ungapped_len() >= min_length]

    end = time.perf_counter()
    print(f"Total runtime: {end - start:0.2f} seconds")
# End main

if __name__ == "__main__":
    main()
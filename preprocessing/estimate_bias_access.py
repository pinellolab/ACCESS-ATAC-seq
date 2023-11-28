import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pandas as pd
import pyranges as pr
import pysam
import logging
from itertools import product
from tqdm import tqdm

from utils import get_chrom_size_from_bam, revcomp


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script estimates bias for Ddd1 based on naked DNA library",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bam_file",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--ref_fasta",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--k_nb", type=int, default=5
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help=(
            "If specified all output files will be written to that directory. \n"
            "Default: the current working directory"
        ),
    )
    parser.add_argument(
        "--name",
        type=str,
        default="counts",
        help=("Names for output file. Default: counts"),
    )

    return parser.parse_args()


def get_access_pwm(bam, grs, fasta, window_size):
    pwm = dict(
        [
            ("A", [0.0] * window_size),
            ("C", [0.0] * window_size),
            ("G", [0.0] * window_size),
            ("T", [0.0] * window_size),
            ("N", [0.0] * window_size),
        ]
    )

    for chrom, start, end in tqdm(zip(grs.Chromosome, grs.Start, grs.End)):
        for read in bam.fetch(reference=chrom, start=start, end=end):
            
            # get reference and query sequence. 
            # for forward reads, it's from 5' to 3', for reverse reads, it's 3' to 5'
            refer_seq = read.get_reference_sequence().upper()
            query_seq = read.query_sequence.upper()
            
            # we only look at reads with substitution
            if len(refer_seq) != len(query_seq):
                continue
            
            for i in range(len(refer_seq)):
                edit_site = read.reference_start + i
                p1 = edit_site - window_size // 2
                p2 = p1 + window_size
            
                motif_seq = None
                if refer_seq[i] == 'C' and query_seq[i] == 'T':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue

                elif refer_seq[i] == 'G' and query_seq[i] == 'A':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue
                    motif_seq = revcomp(motif_seq)
                
                if motif_seq:
                    for j in range(len(motif_seq)):
                        if j == (window_size // 2):
                                continue
                        pwm[motif_seq[j]][j] += 1

    return pwm


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")
    fasta = pysam.FastaFile(args.ref_fasta)

    if args.outdir is None:
        args.outdir = os.getcwd()

    grs = get_chrom_size_from_bam(bam=bam)
    
    alphabet = ["A", "C", "G", "T"]
    kmer_comb = ["".join(e) for e in product(alphabet, repeat=args.k_nb)]
    bias_table = dict([(e, 0.0) for e in kmer_comb])
    
    for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
        for read in bam.fetch(reference=chrom, start=start, end=end):

            # get reference and query sequence. 
            # for forward reads, it's from 5' to 3', for reverse reads, it's 3' to 5'
            refer_seq = read.get_reference_sequence().upper()
            query_seq = read.query_sequence.upper()

            # we only look at reads with substitution
            if len(refer_seq) != len(query_seq):
                continue

            for i in range(len(refer_seq)):
                edit_site = read.reference_start + i
                p1 = edit_site - args.k_nb // 2
                p2 = p1 + args.k_nb
                
                motif_seq = None
                # C -> T at forward strand
                if refer_seq[i] == 'C' and query_seq[i] == 'T':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue

                # C -> T and reverse strand
                elif refer_seq[i] == 'G' and query_seq[i] == 'A':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                        motif_seq = revcomp(motif_seq)
                    except:
                        continue
                     
                if motif_seq and 'N' not in motif_seq:
                    bias_table[motif_seq] += 1

if __name__ == "__main__":
    main()

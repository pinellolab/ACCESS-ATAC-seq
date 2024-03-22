import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging
from itertools import product
import pyranges as pr

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
        "--bed_file",
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
    parser.add_argument("--k_mer", type=int, default=5)
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help=(
            "If specified all output files will be written to that directory. \n"
            "Default: the current working directory"
        ),
    )
    parser.add_argument(
        "--out_name",
        type=str,
        default="counts",
        help=("Names for output file. Default: counts"),
    )

    return parser.parse_args()


def revcomp(s):
    rev_dict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
    return "".join([rev_dict[e] for e in s[::-1]])


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")
    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.bed_file)

    logging.info(f"Total of {len(grs)} regions")

    alphabet = ["A", "C", "G", "T"]
    kmer_comb = ["".join(e) for e in product(alphabet, repeat=args.k_mer)]
    obs_kmer_dict = dict([(e, 0.0) for e in kmer_comb])
    exp_kmer_dict = dict([(e, 0.0) for e in kmer_comb])
    bias_dict = dict([(e, 0.0) for e in kmer_comb])

    for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
        # ger observed kmer
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
                p1 = edit_site - args.k_mer // 2
                p2 = p1 + args.k_mer

                obs_kmer = None
                # C -> T at forward strand
                if refer_seq[i] == "C" and query_seq[i] == "T":
                    try:
                        obs_kmer = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue

                # C -> T and reverse strand
                elif refer_seq[i] == "G" and query_seq[i] == "A":
                    try:
                        obs_kmer = str(fasta.fetch(chrom, p1, p2)).upper()
                        obs_kmer = revcomp(obs_kmer)
                    except:
                        continue
                
                if obs_kmer in obs_kmer_dict:
                    obs_kmer_dict[obs_kmer] += 1
                    
        # get expected kmer
        start = max(start - (args.k_mer // 2), 0)
        ref_seq = str(fasta.fetch(chrom, start, end + (args.k_mer // 2))).upper()
        
        for i in range(len(ref_seq)):
            exp_kmer_forward = ref_seq[i: i + args.k_mer]
            exp_kmer_reverse = revcomp(exp_kmer_forward)
            
            if exp_kmer_forward in exp_kmer_dict.keys():
                exp_kmer_dict[exp_kmer_forward] += 1
                exp_kmer_dict[exp_kmer_reverse] += 1

    # Normalized the observed and expected kmer table as a distribution
    total_count = sum(obs_kmer_dict.values())
    for kmer in obs_kmer_dict.keys():
        obs_kmer_dict[kmer] = obs_kmer_dict[kmer] / total_count
        
    total_count = sum(exp_kmer_dict.values())
    for kmer in exp_kmer_dict.keys():
        exp_kmer_dict[kmer] = exp_kmer_dict[kmer] / total_count

    # Normalize the kmer table as obs/exp
    for kmer in bias_dict.keys():
        if exp_kmer_dict[kmer] > 0:
            bias_dict[kmer] = round(obs_kmer_dict[kmer] / exp_kmer_dict[kmer], 6) 

    # Write the dictionary to a text file
    output_filename = os.path.join(args.out_dir, f"{args.out_name}.obs.txt")
    with open(output_filename, "w") as f:
        for key, value in obs_kmer_dict.items():
            f.write(f"{key}\t{value}\n")
        
    output_filename = os.path.join(args.out_dir, f"{args.out_name}.exp.txt")
    with open(output_filename, "w") as f:
        for key, value in exp_kmer_dict.items():
            f.write(f"{key}\t{value}\n")
            
    output_filename = os.path.join(args.out_dir, f"{args.out_name}.txt")
    with open(output_filename, "w") as f:
        for key, value in bias_dict.items():
            f.write(f"{key}\t{value}\n")
                
    logging.info(f"Done")


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging
import numpy as np
import pyranges as pr
import subprocess


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
        "--peak_file",
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
        "--bias_table_file",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
    )
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


def get_bias_signal_access(
    fasta: pysam.FastaFile, chrom: str, start: int, end: int, kmer_dict: dict
) -> np.array:
    """
    Generate bias signal using k-mer model.

    Parameters
    ----------
    fasta : pysam.FastaFile
        Reference genome
    chrom : str
        Chromosome name
    start : int
        Start site
    end : int
        End site
    kmer_dict : dict
        Dictionary object including bias for each k-mer

    Returns
    -------
    np.array
        Bias signal
    """
    k_nb = len(list(kmer_dict.keys())[0])
    signal_forward = np.zeros(shape=(end - start))
    signal_reverse = np.zeros(shape=(end - start))

    start = max(start - (k_nb // 2), 0)

    seq = str(fasta.fetch(chrom, start, end + (k_nb // 2))).upper()

    for i in range(len(signal_forward)):
        kmer_seq = seq[i : i + k_nb]

        if "N" in kmer_seq:
            continue

        signal_forward[i] = kmer_dict.get(kmer_seq, 0)
        signal_reverse[i] = kmer_dict.get(revcomp(kmer_seq), 0)

    signal = signal_forward + signal_reverse

    return signal


def main():
    args = parse_args()

    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.peak_file)

    logging.info(f"Total of {len(grs)} regions")

    # read bias table
    logging.info(f"Loading bias table")
    kmer_dict = {}
    with open(args.bias_table_file, "r") as file:
        for line in file:
            key, value = line.strip().split("\t")
            kmer_dict[key] = float(value)

    wig_filename = os.path.join(args.out_dir, f"{args.out_name}.wig")
    bw_filename = os.path.join(args.out_dir, f"{args.out_name}.bw")

    # Open a new bigwig file for writing
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            bias_signal = get_bias_signal_access(
                fasta=fasta, chrom=chrom, start=start, end=end, kmer_dict=kmer_dict
            )

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in bias_signal))
            f.write("\n")

    # convert to bigwig file
    subprocess.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

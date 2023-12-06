import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging
import numpy as np
import pyranges as pr
from utils import revcomp, wig_to_bw


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
    parser.add_argument("--ext", type=int, default=50)
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


def main():
    args = parse_args()

    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.peak_file)

    # grs = grs.extend(args.ext)
    # grs = grs.merge()
    logging.info(f"Total of {len(grs)} regions")

    # read bias table
    logging.info(f"Loading bias table")
    kmer_dict = {}
    with open(args.bias_table_file, "r") as file:
        for line in file:
            key, value = line.strip().split("\t")
            kmer_dict[key] = float(value)

    k_nb = len(list(kmer_dict.keys())[0])

    wig_forward_filename = os.path.join(
        args.out_dir, "{}.forward.wig".format(args.out_name)
    )
    wig_reverse_filename = os.path.join(
        args.out_dir, "{}.reverse.wig".format(args.out_name)
    )
    bw_forward_filename = os.path.join(
        args.out_dir, "{}.forward.bw".format(args.out_name)
    )
    bw_reverse_filename = os.path.join(
        args.out_dir, "{}.reverse.bw".format(args.out_name)
    )

    # Open a new bigwig file for writing
    with open(wig_forward_filename, "a") as forward_file, open(
        wig_reverse_filename, "a"
    ) as reverse_file:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            bias_signal_forward = np.zeros(shape=(end - start))
            bias_signal_reverse = np.zeros(shape=(end - start))
            seq = str(
                fasta.fetch(chrom, start - (k_nb // 2), end + (k_nb // 2))
            ).upper()

            for i in range(len(bias_signal_forward)):
                kmer_seq = seq[i : i + k_nb]

                if "N" in kmer_seq:
                    continue

                bias_signal_forward[i] = kmer_dict[kmer_seq]
                bias_signal_reverse[i] = kmer_dict[revcomp(kmer_seq)]

            forward_file.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            forward_file.write("\n".join(str(e) for e in bias_signal_forward))
            forward_file.write("\n")

            reverse_file.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            reverse_file.write("\n".join(str(e) for e in bias_signal_reverse))
            reverse_file.write("\n")

    # convert to bigwig file
    wig_to_bw(
        wig_filename=wig_forward_filename,
        bw_filename=bw_forward_filename,
        chrom_size_file=args.chrom_size_file,
    )

    wig_to_bw(
        wig_filename=wig_reverse_filename,
        bw_filename=bw_reverse_filename,
        chrom_size_file=args.chrom_size_file,
    )


if __name__ == "__main__":
    main()

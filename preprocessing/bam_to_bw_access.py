import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pyranges as pr
import pysam
import logging


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


from utils import get_chrom_size_from_bam, wig_to_bw
from get_signal import get_raw_signal_access


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates BigWig file from a BAM file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bam_file",
        type=str,
        default=None,
        help=("BAM file containing the aligned reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--ref_fasta",
        type=str,
        default=None,
        help=("FASTA file containing reference genome. \n" "Default: None"),
    )
    parser.add_argument(
        "--peak_file",
        type=str,
        default=None,
        help=(
            "BED file containing genomic regions for generating signal. \n"
            "Default: None"
        ),
    )
    parser.add_argument("--signal", type=str, choices=["raw", "exp"], default="raw")
    parser.add_argument(
        "--window_size",
        type=int,
        default=11,
        help=("Number of base pairs for smoothing signal. Default: 11"),
    )
    parser.add_argument(
        "--bias_table_file",
        type=str,
        default=None,
        help=(
            "TXT file containing enzyme bias for each k-mer. This file is required to for bias correction. \n"
            "Default: None"
        ),
    )
    parser.add_argument("--ext", type=int, default=50)
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
    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    if args.peak_file:
        logging.info(f"Loading genomic regions from {args.peak_file}")
        grs = pr.read_bed(args.peak_file)
        grs = grs.merge()
    else:
        logging.info(f"Using whole genome")
        grs = get_chrom_size_from_bam(bam=bam)

    logging.info(f"Total of {len(grs)} regions")

    if args.ref_fasta:
        fasta = pysam.FastaFile(args.ref_fasta)

    if args.bias_table_file:
        logging.info(f"Loading bias table")
        kmer_dict = {}
        with open(args.bias_table_file, "r") as file:
            for line in file:
                key, value = line.strip().split("\t")
                kmer_dict[key] = float(value)

    wig_forward = os.path.join(
        args.out_dir, "{}.wig".format(args.out_name)
    )
    wig_forward_filename = os.path.join(
        args.out_dir, "{}.forward.wig".format(args.out_name)
    )
    wig_reverse_filename = os.path.join(
        args.out_dir, "{}.reverse.wig".format(args.out_name)
    )
    bw_filename = os.path.join(
        args.out_dir, "{}.bw".format(args.out_name)
    )
    bw_forward_filename = os.path.join(
        args.out_dir, "{}.forward.bw".format(args.out_name)
    )
    bw_reverse_filename = os.path.join(
        args.out_dir, "{}.reverse.bw".format(args.out_name)
    )

    with open(wig_forward_filename, "a") as forward_file, open(
        wig_reverse_filename, "a"
    ) as reverse_file:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            if args.signal == "raw":
                signal_forward, signal_reverse = get_raw_signal_access(
                    chrom, start, end, bam
                )
                signal = signal_forward + signal_reverse

            forward_file.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            forward_file.write("\n".join(str(e) for e in signal_forward))
            forward_file.write("\n")

            reverse_file.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            reverse_file.write("\n".join(str(e) for e in signal_reverse))
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

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
        "--peak_file",
        type=str,
        default=None,
        help=(
            "BED file containing genomic regions for generating signal. \n"
            "If none, will use the whole genome as input regions. \n"
            "Default: None"
        ),
    )
    parser.add_argument(
        "--signal", type=str, choices=["raw", "smooth", "bias_correct"], default="raw"
    )

    parser.add_argument(
        "--smooth_window",
        type=int,
        default=11,
        help=("Number of base pairs for smoothing signal. Default: 11"),
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
    


def output_raw_smooth(bam, grs, wig_filename, smooth_window):
    half_smooth_window = smooth_window // 2

    with open(wig_filename, "a") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            signal = get_raw_edit_counts(
                chrom=chrom,
                start=start - half_smooth_window,
                end=end + half_smooth_window,
                bam=bam,
            )

            w = np.ones(half_smooth_window, "d")
            signal_ = np.convolve(w, signal, mode="valid")
            signal_ = signal_[half_smooth_window:-half_smooth_window]

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in signal_))
            f.write("\n")


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

    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))

    if args.signal == "raw":
        output_raw(bam=bam, grs=grs, wig_filename=wig_filename)
    elif args.signal == "smooth":
        output_raw_smooth(
            bam=bam,
            grs=grs,
            wig_filename=wig_filename,
            smooth_window=args.smooth_window,
        )

    # convert to bigwig file
    wig_to_bw(
        wig_filename=wig_filename,
        bw_filename=bw_filename,
        chrom_size_file=args.chrom_size_file,
    )


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pyranges as pr
import pyBigWig
import logging
import subprocess


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generated normalized signal",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bw_raw_file",
        type=str,
        default=None,
        help=("Bigwig file containing the raw signal. \n" "Default: None"),
    )
    parser.add_argument(
        "--bw_bias_file",
        type=str,
        default=None,
        help=("Bigwig file containing the bias signal. \n" "Default: None"),
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

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    bw_raw = pyBigWig.open(args.bw_raw_file)
    bw_bias = pyBigWig.open(args.bw_bias_file)

    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))

    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            raw_signal = np.array(bw_raw.values(chrom, start, end))
            raw_bias = np.array(bw_bias.values(chrom, start, end))

            raw_signal[np.isnan(raw_signal)] = 0
            raw_bias[np.isnan(raw_bias)] = 0

            raw_bias = raw_bias / np.sum(raw_bias)
            exp_signal = np.sum(raw_signal) * raw_bias
            norm_signal = np.divide(raw_signal, exp_signal)

            norm_signal = np.log2(norm_signal + 0.1)

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in norm_signal))
            f.write("\n")

    # convert to bigwig file
    subprocess.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    # os.remove(wig_filename)

if __name__ == "__main__":
    main()

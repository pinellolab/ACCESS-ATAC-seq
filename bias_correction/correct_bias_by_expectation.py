import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pyranges as pr
import pyBigWig
import logging
import subprocess as sp


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
    parser.add_argument("--bw_raw", type=str, default=None)
    parser.add_argument("--bw_bias", type=str, default=None)
    parser.add_argument(
        "--peak_file",
        type=str,
        default=None,
        help=(
            "BED file containing genomic regions for generating signal. \n"
            "Default: None"
        ),
    )
    parser.add_argument("--peak_extend", type=int, default=100)
    parser.add_argument("--window", type=int, default=100)
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
        default=None,
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
    grs = grs.extend(args.peak_extend)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    # load raw and bias signal
    bw_raw = pyBigWig.open(args.bw_raw)
    bw_bias = pyBigWig.open(args.bw_bias)

    # generate expected signal based on 100bp bins
    logging.info("Generating expected signal!")
    exp_wig_filename = f"{args.out_dir}/{args.out_name}.exp.wig"
    norm_wig_filename = f"{args.out_dir}/{args.out_name}.norm.wig"
    exp_bw_filename = f"{args.out_dir}/{args.out_name}.exp.bw"
    norm_bw_filename = f"{args.out_dir}/{args.out_name}.norm.bw"

    with open(exp_wig_filename, "w") as exp_f, open(norm_wig_filename, "w") as norm_f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            n = (end - start) // args.window

            for i in range(n):
                _start, _end = start + i * args.window, start + (i + 1) * args.window

                if _end > end:
                    _end = end

                # read raw and bias signal
                signal_raw = np.array(bw_raw.values(chrom, _start, _end))
                signal_bias = np.array(bw_bias.values(chrom, _start, _end))

                # NAN to zero if needed
                signal_raw[np.isnan(signal_raw)] = 0
                signal_bias[np.isnan(signal_bias)] = 0

                # normalize bias signal
                if np.sum(signal_bias) > 0:
                    signal_bias = signal_bias / np.sum(signal_bias)

                signal_exp = np.sum(signal_raw) * signal_bias

                signal_norm = np.divide(signal_raw + 1, signal_exp + 1)
                # signal_norm = np.log2(signal_norm)

                # save signal to wig files
                exp_f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                exp_f.write("\n".join(str(e) for e in signal_exp))
                exp_f.write("\n")

                norm_f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                norm_f.write("\n".join(str(e) for e in signal_norm))
                norm_f.write("\n")

    # convert to bigwig file
    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", exp_wig_filename, args.chrom_size_file, exp_bw_filename])
    sp.run(["wigToBigWig", norm_wig_filename, args.chrom_size_file, norm_bw_filename])
    os.remove(exp_wig_filename)
    os.remove(norm_wig_filename)

    logging.info("Done!")


if __name__ == "__main__":
    main()

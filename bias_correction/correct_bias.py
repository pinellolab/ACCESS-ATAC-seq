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
    parser.add_argument("--window", type=int, default=101)
    parser.add_argument("--pseudo_count", type=float, default=1)
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

    half_window = args.window // 2

    with open(exp_wig_filename, "w") as exp_f, open(norm_wig_filename, "w") as norm_f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            signal_raw = np.array(
                bw_raw.values(chrom, start - half_window, end + half_window)
            )
            signal_bias = np.array(
                bw_bias.values(chrom, start - half_window, end + half_window)
            )

            # NAN to zero if needed
            signal_raw[np.isnan(signal_raw)] = 0
            signal_bias[np.isnan(signal_bias)] = 0

            # get expected values
            signal_exp = np.zeros(shape=(end - start))
            for i in range(half_window, len(signal_raw) - half_window):
                total_signal = np.sum(signal_raw[i - half_window : i + half_window])
                total_bias = np.sum(signal_bias[i - half_window : i + half_window])

                if total_bias == 0:
                    signal_exp[i - half_window] = 0
                else:
                    signal_exp[i - half_window] = (
                        total_signal * signal_bias[i] / total_bias
                    )

            signal_norm = np.divide(
                signal_raw[half_window:-half_window] + args.pseudo_count,
                signal_exp + args.pseudo_count,
            )

            signal_exp = np.round(signal_exp, decimals=3)
            signal_norm = np.round(signal_norm, decimals=3)

            # save signal to wig files
            exp_f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            exp_f.write("\n".join(str(e) for e in signal_exp))
            exp_f.write("\n")

            norm_f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
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

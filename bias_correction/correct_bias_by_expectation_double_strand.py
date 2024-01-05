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
    parser.add_argument("--bw_raw_forward", type=str, default=None)
    parser.add_argument("--bw_raw_reverse", type=str, default=None)
    parser.add_argument("--bw_bias_forward", type=str, default=None)
    parser.add_argument("--bw_bias_reverse", type=str, default=None)
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
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    # load raw and bias signal for both strands
    bw_raw_forward = pyBigWig.open(args.bw_raw_forward)
    bw_raw_reverse = pyBigWig.open(args.bw_raw_reverse)

    bw_bias_forward = pyBigWig.open(args.bw_bias_forward)
    bw_bias_reverse = pyBigWig.open(args.bw_bias_reverse)

    # generate expected signal based on 100bp bins
    logging.info("Generating expected signal!")
    wig_filename = f"{args.out_dir}/{args.out_name}.exp.wig"
    wig_forward_filename = f"{args.out_dir}/{args.out_name}.exp.forward.wig"
    wig_reverse_filename = f"{args.out_dir}/{args.out_name}.exp.reverse.wig"

    bw_filename = f"{args.out_dir}/{args.out_name}.exp.bw"
    bw_forward_filename = f"{args.out_dir}/{args.out_name}.exp.forward.bw"
    bw_reverse_filename = f"{args.out_dir}/{args.out_name}.exp.reverse.bw"

    with open(wig_forward_filename, "w") as forward_file, open(
        wig_reverse_filename, "w"
    ) as reverse_file, open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            n = (end - start) // 100

            for i in range(n):
                _start, _end = start + i * 100, start + (i + 1) * 100

                if _end > end:
                    _end = end

                # read raw signal
                signal_raw_forward = np.array(
                    bw_raw_forward.values(chrom, _start, _end)
                )
                signal_raw_reverse = np.array(
                    bw_raw_reverse.values(chrom, _start, _end)
                )

                # read bias signal
                signal_bias_forward = np.array(
                    bw_bias_forward.values(chrom, _start, _end)
                )
                signal_bias_reverse = np.array(
                    bw_bias_reverse.values(chrom, _start, _end)
                )

                # NAN to zero if needed
                signal_raw_forward[np.isnan(signal_raw_forward)] = 0
                signal_raw_reverse[np.isnan(signal_raw_reverse)] = 0
                signal_bias_forward[np.isnan(signal_bias_forward)] = 0
                signal_bias_reverse[np.isnan(signal_bias_reverse)] = 0

                # normalize bias signal
                if np.sum(signal_bias_forward) > 0:
                    signal_bias_forward = signal_bias_forward / np.sum(
                        signal_bias_forward
                    )

                if np.sum(signal_bias_reverse) > 0:
                    signal_bias_reverse = signal_bias_reverse / np.sum(
                        signal_bias_reverse
                    )

                exp_signal_forward = np.sum(signal_raw_forward) * signal_bias_forward
                exp_signal_reverse = np.sum(signal_raw_reverse) * signal_bias_reverse
                exp_signal = exp_signal_forward + exp_signal_reverse

                # save signal to wig files
                forward_file.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                forward_file.write("\n".join(str(e) for e in exp_signal_forward))
                forward_file.write("\n")

                reverse_file.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                reverse_file.write("\n".join(str(e) for e in exp_signal_reverse))
                reverse_file.write("\n")

                f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                f.write("\n".join(str(e) for e in exp_signal))
                f.write("\n")

    # convert to bigwig file
    subprocess.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    subprocess.run(
        ["wigToBigWig", wig_forward_filename, args.chrom_size_file, bw_forward_filename]
    )
    subprocess.run(
        ["wigToBigWig", wig_reverse_filename, args.chrom_size_file, bw_reverse_filename]
    )

    os.remove(wig_filename)
    os.remove(wig_forward_filename)
    os.remove(wig_reverse_filename)

    logging.info("Generating normalized signal!")
    wig_filename = f"{args.out_dir}/{args.out_name}.norm.wig"
    wig_forward_filename = f"{args.out_dir}/{args.out_name}.norm.forward.wig"
    wig_reverse_filename = f"{args.out_dir}/{args.out_name}.norm.reverse.wig"

    bw_filename = f"{args.out_dir}/{args.out_name}.norm.bw"
    bw_forward_filename = f"{args.out_dir}/{args.out_name}.norm.forward.bw"
    bw_reverse_filename = f"{args.out_dir}/{args.out_name}.norm.reverse.bw"

    


if __name__ == "__main__":
    main()

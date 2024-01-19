import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import subprocess

import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates true labels for TF binding sites based on ChIP-seq peaks",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bw_obs_file", type=str, default=None)
    parser.add_argument("--bw_exp_file", type=str, default=None)
    parser.add_argument("--fp_window", type=int, default=None)
    parser.add_argument("--flank_window", type=int, default=None)
    parser.add_argument("--smooth_window", type=int, default=None)
    parser.add_argument(
        "--peak_file",
        type=str,
        default=100,
        help=("Flanking regions used to estimate footprint score. Default: 20"),
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

    bw_obs = pyBigWig.open(args.bw_obs_file)
    bw_exp = pyBigWig.open(args.bw_exp_file)

    extend = args.fp_window // 2 + args.flank_window

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    # compute footprint score per-nucleotide
    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))

    logging.info(f"Calculating footprint score per-nucleotide")
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            _start = start - extend - args.smooth_window // 2
            _end = end + extend + args.smooth_window // 2

            signal_obs = np.array(bw_obs.values(chrom, _start, _end))
            signal_exp = np.array(bw_exp.values(chrom, _start, _end))

            # smooth the observed and expected signal by moving average
            w = np.ones(args.smooth_window)
            signal_obs_smooth = np.convolve(w / w.sum(), signal_obs, mode="valid")
            signal_exp_smooth = np.convolve(w / w.sum(), signal_exp, mode="valid")

            signal_norm = signal_obs_smooth - signal_exp_smooth

            # compute footprint score per-nucleotide
            # fp_socre = mean(flank signal) - mean(footprint signal)
            fp_scores = np.zeros(end - start)
            for i in range(extend, len(signal_norm) - extend):
                fp_start, fp_end = i - args.fp_window // 2, i + args.fp_window // 2
                left_flank_start = i - args.fp_window // 2 - args.flank_window
                left_flank_end = i - args.fp_window // 2
                right_blank_start = i + args.fp_window // 2
                right_blank_end = i + args.fp_window // 2 + args.flank_window

                fp_signal = signal_norm[fp_start:fp_end]
                left_flank_signal = signal_norm[left_flank_start:left_flank_end]
                right_flank_signal = signal_norm[right_blank_start:right_blank_end]

                flank_signal = np.concatenate([left_flank_signal, right_flank_signal])

                fp_scores[i - extend] = np.mean(flank_signal) - np.mean(fp_signal)

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in fp_scores))
            f.write("\n")

    # convert to bigwig file
    subprocess.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    
    # os.remove(wig_filename)
    
    logging.info(f"Done!")

if __name__ == "__main__":
    main()

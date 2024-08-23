import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import subprocess as sp
from scipy.stats import norm, false_discovery_control
import warnings

warnings.filterwarnings("ignore")

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
    parser.add_argument("--peak_file", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    parser.add_argument("--chrom_size_file", type=str, default=None)
    return parser.parse_args()

def main():
    args = parse_args()

    bw_obs = pyBigWig.open(args.bw_obs_file)
    bw_exp = pyBigWig.open(args.bw_exp_file)

    extend = args.fp_window // 2 + args.flank_window
    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    # for test, we only run the process for chromosome 1
    # grs = grs[grs.Chromosome.isin(["chr1"])]

    logging.info(f"Total of {len(grs)} regions")

    # compute footprint score per-nucleotide
    obs_wig_filename = os.path.join(
        args.out_dir, "{}.obs.wig".format(args.out_name))
    obs_bw_filename = os.path.join(
        args.out_dir, "{}.obs.bw".format(args.out_name))
    exp_wig_filename = os.path.join(
        args.out_dir, "{}.exp.wig".format(args.out_name))
    exp_bw_filename = os.path.join(
        args.out_dir, "{}.exp.bw".format(args.out_name))

    logging.info(f"Calculating footprint score per-nucleotide")

    with open(obs_wig_filename, "w") as obs_f, open(
        exp_wig_filename, "w"
    ) as exp_f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            _start = start - extend - args.smooth_window // 2
            _end = end + extend + args.smooth_window // 2

            signal_obs = np.array(bw_obs.values(chrom, _start, _end))
            signal_exp = np.array(bw_exp.values(chrom, _start, _end))

            # smooth the observed and expected signal by moving average
            w = np.ones(args.smooth_window)
            signal_obs_smooth = np.convolve(
                w / w.sum(), signal_obs, mode="valid")
            signal_exp_smooth = np.convolve(
                w / w.sum(), signal_exp, mode="valid")

            # normalize the signal by average
            signal_obs_norm = signal_obs_smooth / np.mean(signal_obs_smooth)
            signal_exp_norm = signal_exp_smooth / np.mean(signal_exp_smooth)

            obs_f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
            obs_f.write("\n".join(str(e) for e in signal_obs_norm))
            obs_f.write("\n")

            exp_f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
            exp_f.write("\n".join(str(e) for e in signal_exp_norm))
            exp_f.write("\n")

    # convert to bigwig file
    sp.run(["wigToBigWig", obs_wig_filename,
           args.chrom_size_file, obs_bw_filename])
    sp.run(["wigToBigWig", exp_wig_filename,
           args.chrom_size_file, exp_bw_filename])

    os.remove(obs_wig_filename)
    os.remove(exp_wig_filename)

    logging.info(f"Done!")


if __name__ == "__main__":
    main()

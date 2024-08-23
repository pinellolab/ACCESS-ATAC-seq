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


def compute_fp_score(
    signal: np.array = None,
    peak_len: int = None,
    extend: int = None,
    fp_window: int = None,
    flank_window: int = None,
) -> np.array:
    """
    Compute footprint score for each position based on a local window.
    Specifically, fp_score = mean(flank signal) - mean(footprint signal)

    Parameters
    ----------
    signal : np.array, optional
        Input signal, by default None
    peak_len : int, optional
        _description_, by default None
    extend : int, optional
        _description_, by default None
    fp_window : int, optional
        _description_, by default None
    flank_window : int, optional
        _description_, by default None

    Returns
    -------
    np.array
        _description_
    """

    fp_scores = np.zeros(peak_len)

    for i in range(extend, len(signal) - extend):
        fp_signal = signal[i - fp_window // 2: i + fp_window // 2]
        left_flank_signal = signal[
            i - fp_window // 2 - flank_window: i - fp_window // 2
        ]
        right_flank_signal = signal[
            i + fp_window // 2: i + fp_window // 2 + flank_window
        ]

        flank_signal = np.concatenate([left_flank_signal, right_flank_signal])

        fp_scores[i - extend] = np.mean(flank_signal) - np.mean(fp_signal)

    return fp_scores


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
        args.out_dir, "{}.obs.fs.wig".format(args.out_name))
    obs_bw_filename = os.path.join(
        args.out_dir, "{}.obs.fs.bw".format(args.out_name))
    exp_wig_filename = os.path.join(
        args.out_dir, "{}.exp.fs.wig".format(args.out_name))
    exp_bw_filename = os.path.join(
        args.out_dir, "{}.exp.fs.bw".format(args.out_name))
    bc_wig_filename = os.path.join(
        args.out_dir, "{}.bc.fs.wig".format(args.out_name))
    bc_bw_filename = os.path.join(
        args.out_dir, "{}.bc.fs.bw".format(args.out_name))
    bc_wig_filename = os.path.join(
        args.out_dir, "{}.bc.fs.wig".format(args.out_name))
    bc_bw_filename = os.path.join(
        args.out_dir, "{}.bc.fs.bw".format(args.out_name))

    logging.info(f"Calculating footprint score per-nucleotide")

    with open(obs_wig_filename, "w") as obs_f, open(
        exp_wig_filename, "w"
    ) as exp_f, open(bc_wig_filename, "w") as bc_f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            _start = start - extend 
            _end = end + extend

            signal_obs = np.array(bw_obs.values(chrom, _start, _end))
            signal_exp = np.array(bw_exp.values(chrom, _start, _end))

            # compute footprint score per-nucleotide
            obs_fp_scores = compute_fp_score(
                signal=signal_obs,
                peak_len=end - start,
                extend=extend,
                fp_window=args.fp_window,
                flank_window=args.flank_window,
            )

            # We also need to estimate background footprint score per nucleotide,
            # so we can perform statistical test to get p-value.
            # To do this, we randomly shuffle the input signal, and then compute the
            # footprint score again.
            exp_fp_scores = compute_fp_score(
                signal=signal_exp,
                peak_len=end - start,
                extend=extend,
                fp_window=args.fp_window,
                flank_window=args.flank_window,
            )

            bc_fp_score = obs_fp_scores - exp_fp_scores

            obs_f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            obs_f.write("\n".join(str(e) for e in obs_fp_scores))
            obs_f.write("\n")

            exp_f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            exp_f.write("\n".join(str(e) for e in exp_fp_scores))
            exp_f.write("\n")

            bc_f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            bc_f.write("\n".join(str(e) for e in bc_fp_score))
            bc_f.write("\n")

    # convert to bigwig file
    sp.run(["wigToBigWig", obs_wig_filename,
           args.chrom_size_file, obs_bw_filename])
    sp.run(["wigToBigWig", exp_wig_filename,
           args.chrom_size_file, exp_bw_filename])
    sp.run(["wigToBigWig", bc_wig_filename, args.chrom_size_file, bc_bw_filename])

    os.remove(obs_wig_filename)
    os.remove(exp_wig_filename)
    os.remove(bc_wig_filename)

    logging.info(f"Done!")


if __name__ == "__main__":
    main()

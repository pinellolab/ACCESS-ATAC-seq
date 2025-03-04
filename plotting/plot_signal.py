import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import argparse
import pyranges as pr
import pysam
import pyBigWig
import logging
import matplotlib.pyplot as plt
import logomaker

from utils import revcomp, get_motif_df

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates footprint plot using BW file and BED file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bw_file",
        type=str,
        default=None,
        help=("BIGWIG file containing the signal. \n" "Default: None"),
    )

    parser.add_argument(
        "--ref_fasta",
        type=str,
        default=None,
        help=("BIGWIG file containing the signal. \n" "Default: None"),
    )

    parser.add_argument(
        "--motif_file",
        type=str,
        default=None,
        help=("BED file containing genomic regions for plotting. \n" "Default: None"),
    )

    parser.add_argument(
        "--bed_file",
        type=str,
        default=None,
        help=("BED file containing genomic regions for plotting. \n" "Default: None"),
    )

    parser.add_argument(
        "--extend", type=int, default=500, help=("Extend the regions. Default: 100")
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


def get_signal(bw, grs):
    window_size = grs.End.values[0] - grs.Start.values[0]
    signal = np.zeros(shape=(len(grs), window_size))

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        signal[i] = bw.values(chrom, start, end)

    signal[np.isnan(signal)] = 0
    signal = np.mean(signal, axis=0)

    return signal


def get_pwm(fasta, grs):
    window_size = grs.End.values[0] - grs.Start.values[0]

    pwm = dict(
        [
            ("A", [0.0] * window_size),
            ("C", [0.0] * window_size),
            ("G", [0.0] * window_size),
            ("T", [0.0] * window_size),
            ("N", [0.0] * window_size),
        ]
    )

    for i, (chrom, start, end, strand) in enumerate(
        zip(grs.Chromosome, grs.Start, grs.End, grs.Strand)
    ):
        seq = str(fasta.fetch(chrom, start, end)).upper()

        if strand == "-":
            seq = revcomp(seq)

            for i in range(len(seq)):
                pwm[seq[i]][i] += 1
                
    return pwm


def main():
    args = parse_args()

    logging.info(f"Loading genomic regions from {args.bed_file}")
    grs = pr.read_bed(args.bed_file)

    # if args.peak_file:
    #     grs_peaks = pr.read_bed(args.peak_file)
    #     grs = grs.overlap(grs_peaks, strandedness=False, invert=False)

    logging.info(f"Total of {len(grs)} regions")

    # extend regions
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - args.extend
    grs.End = mid + args.extend

    window_size = grs.End.values[0] - grs.Start.values[0]

    bw = pyBigWig.open(args.bw_file)
    fasta = pysam.FastaFile(args.ref_fasta)

    logging.info("Generating signal and pwm")
    signal = get_signal(bw=bw, grs=grs)
    pwm = get_pwm(fasta=fasta, grs=grs)

    df = get_motif_df(pwm)

    logging.info("Plotting")
    start = -(window_size // 2)
    end = (window_size // 2) - 1
    x = np.linspace(start, end, num=window_size)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, signal, color="#4daf4a")

    ax.text(
        0.15,
        0.9,
        "n = {}".format(len(grs)),
        verticalalignment="bottom",
        horizontalalignment="right",
        transform=ax.transAxes,
        fontweight="bold",
    )

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 15))

    ax.tick_params(direction="out")
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    min_signal = np.min(signal)
    max_signal = np.max(signal)
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels(
        [str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90
    )

    # ax.set_title(mpbs_name, fontweight='bold')
    ax.set_xlim(start, end)
    ax.set_ylim([min_signal, max_signal])
    ax.legend(loc="upper right", frameon=False)
    ax.spines["bottom"].set_position(("outward", 70))

    ax = plt.axes([0.105, 0.085, 0.85, 0.2])
    logo = logomaker.Logo(df, ax=ax, show_spines=False, baseline_width=0)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()

    output_filename = os.path.join(args.out_dir, "{}.png".format(args.out_name))
    plt.savefig(output_filename)
    
    df = pd.DataFrame(data={"signal": signal, "position":range(-args.extend, args.extend)})
    df.to_csv(f"{args.out_dir}/{args.out_name}.csv", index=False)


if __name__ == "__main__":
    main()

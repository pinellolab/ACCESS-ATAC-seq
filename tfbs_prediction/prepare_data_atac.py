import argparse
import os
import sys
import numpy as np
import warnings
import pysam
import logging
import pyBigWig
import pyranges as pr

from utils import one_hot_encode

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

if not sys.warnoptions:
    warnings.simplefilter("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generated input data for CNN model",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bw_raw",
        type=str,
        default=None,
        help=(
            "Bigwig file containing the raw signal. \n"
            "For ATAC, this should be counts of Tn5 cutting. \n"
            "Default: None"
        ),
    )
    parser.add_argument(
        "--bw_bias",
        type=str,
        default=None,
        help=("Bigwig file containing the bias signal. \n" "Default: None"),
    )
    parser.add_argument(
        "--regions",
        type=str,
        default=None,
        help="BED file containing genomic regions",
    )
    parser.add_argument("--window", type=int, default=100)
    parser.add_argument("--pseudo_count", type=float, default=1)
    parser.add_argument(
        "--ref_fasta", type=str, default=None, help="FASTQ file for reference genome"
    )
    parser.add_argument("--out_dir", type=str,
                        default=None, help="Output directory")
    parser.add_argument("--out_name", type=str,
                        default=None, help="Output name")
    return parser.parse_args()


def get_norm_signal(bw_raw, bw_bias, chrom, start, end, half_window, pseudo_count):
    signal_raw = bw_raw.values(chrom, start - half_window, end + half_window)
    signal_bias = bw_bias.values(chrom, start - half_window, end + half_window)

    signal_raw = np.array(signal_raw)
    signal_bias = np.array(signal_bias)

    signal_raw[np.isnan(signal_raw)] = 0
    signal_bias[np.isnan(signal_bias)] = 0

    # get expected signal
    signal_exp = np.zeros(shape=(end - start))
    for i in range(half_window, len(signal_raw) - half_window):
        total_signal = np.sum(signal_raw[i - half_window: i + half_window])
        total_bias = np.sum(signal_bias[i - half_window: i + half_window])

        if total_bias == 0:
            signal_exp[i - half_window] = 0
        else:
            signal_exp[i - half_window] = total_signal * \
                signal_bias[i] / total_bias

    # normalized signal = obs/exp
    signal_norm = np.divide(
        signal_raw[half_window:-half_window] +
        pseudo_count, signal_exp + pseudo_count
    )
    signal_norm = np.log2(signal_norm + 0.1)

    return signal_norm


def main():
    args = parse_args()

    # logging.info("Loading input files")
    bw_raw = pyBigWig.open(args.bw_raw)
    bw_bias = pyBigWig.open(args.bw_bias)
    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.regions)

    # extend 100 bps to both sides
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - 100
    grs.End = mid + 100

    # logging.info("Generating data")
    dat = np.empty(shape=(len(grs), 200, 6), dtype=np.float32)

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        signal_raw = np.array(bw_raw.values(chrom, start, end))
        signal_bias = np.array(bw_bias.values(chrom, start, end))

        # make no NAN
        signal_raw[np.isnan(signal_raw)] = 0
        signal_bias[np.isnan(signal_bias)] = 0

        signal_raw = np.expand_dims(signal_raw, axis=1)
        signal_bias = np.expand_dims(signal_bias, axis=1)

        # signal_normal = get_norm_signal(
        #     bw_raw=bw_raw,
        #     bw_bias=bw_bias,
        #     chrom=chrom,
        #     start=start,
        #     end=end,
        #     half_window=half_window,
        #     pseudo_count=args.pseudo_count,
        # )
        # signal = np.concatenate([signal_normal, signal_total])
        # signal_bias = np.expand_dims(signal, axis=1)

        seq = str(fasta.fetch(chrom, start, end)).upper()
        x = one_hot_encode(seq=seq)

        # assemble data
        dat[i] = np.concatenate([x, signal_bias, signal_raw], axis=1)

    # save data
    np.savez_compressed(
        f"{args.out_dir}/{args.out_name}.npz", x=dat, y=np.array(grs.ThickStart)
    )

    # logging.info("Finished!")


if __name__ == "__main__":
    main()
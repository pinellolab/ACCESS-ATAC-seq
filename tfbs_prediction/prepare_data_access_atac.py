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
        "--atac_bw_raw",
        type=str,
        default=None,
        help=(
            "Bigwig file containing the raw signal. \n"
            "Default: None"
        ),
    )
    parser.add_argument(
        "--access_bw_raw",
        type=str,
        default=None,
        help=(
            "Bigwig file containing the raw edit fraction signal. \n"
            "Default: None"
        ),
    )
    parser.add_argument(
        "--atac_bw_bias",
        type=str,
        default=None,
        help=("Bigwig file containing the bias signal. \n" "Default: None"),
    )
    parser.add_argument(
        "--access_bw_bias",
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
    parser.add_argument(
        "--ref_fasta", type=str, default=None, help="FASTQ file for reference genome"
    )
    parser.add_argument("--out_dir", type=str,
                        default=None, help="Output directory")
    parser.add_argument("--out_name", type=str,
                        default=None, help="Output name")
    return parser.parse_args()


def main():
    args = parse_args()

    # logging.info("Loading input files")
    atac_bw_raw = pyBigWig.open(args.atac_bw_raw)
    access_bw_raw = pyBigWig.open(args.access_bw_raw)
    atac_bw_bias = pyBigWig.open(args.atac_bw_bias)
    access_bw_bias = pyBigWig.open(args.access_bw_bias)

    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.regions)

    # extend 100 bps to both sides
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - 100
    grs.End = mid + 100

    # logging.info("Generating data")
    dat = np.empty(shape=(len(grs), 200, 8), dtype=np.float32)

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        atac_signal_raw = np.array(atac_bw_raw.values(chrom, start, end))
        access_signal_raw = np.array(access_bw_raw.values(chrom, start, end))

        atac_bias = np.array(atac_bw_bias.values(chrom, start, end))
        access_bias = np.array(access_bw_bias.values(chrom, start, end))

        # remove NAN
        atac_signal_raw[np.isnan(atac_signal_raw)] = 0
        access_signal_raw[np.isnan(access_signal_raw)] = 0
        atac_bias[np.isnan(atac_bias)] = 0
        access_bias[np.isnan(access_bias)] = 0

        atac_signal_raw = np.expand_dims(atac_signal_raw, axis=1)
        atac_bias = np.expand_dims(atac_bias, axis=1)
        access_signal_raw = np.expand_dims(access_signal_raw, axis=1)
        access_bias = np.expand_dims(access_bias, axis=1)

        seq = str(fasta.fetch(chrom, start, end)).upper()
        x = one_hot_encode(seq=seq)

        # assemble data
        dat[i] = np.concatenate(
            [x, atac_bias, atac_signal_raw, access_bias, access_signal_raw], axis=1)

    # save data
    np.savez_compressed(
        f"{args.out_dir}/{args.out_name}.npz", x=dat, y=np.array(grs.ThickStart)
    )

    logging.info("Finished!")


if __name__ == "__main__":
    main()

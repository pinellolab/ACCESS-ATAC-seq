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
    parser.add_argument("--regions", type=str, default=None)
    parser.add_argument("--ref_fasta", type=str, default=None)
    parser.add_argument("--out_dir", type=str,
                        default=None, help="Output directory")
    parser.add_argument("--out_name", type=str,
                        default=None, help="Output name")
    return parser.parse_args()


def main():
    args = parse_args()

    fasta = pysam.FastaFile(args.ref_fasta)
    grs = pr.read_bed(args.regions)

    # extend 64 bps to both sides
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - 64
    grs.End = mid + 64

    logging.info("Reading cutting counts")

    # logging.info("Generating data")
    dat = np.empty(shape=(len(grs), 128, 4), dtype=np.float32)
    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        seq = str(fasta.fetch(chrom, start, end)).upper()
        dat[i] = one_hot_encode(seq=seq)

    # save data
    np.savez_compressed(
        f"{args.out_dir}/{args.out_name}.npz",
        x=dat,
        y=np.array(grs.ThickStart)
    )

    logging.info("Finished!")


if __name__ == "__main__":
    main()

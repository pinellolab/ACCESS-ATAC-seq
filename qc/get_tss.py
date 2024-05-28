import os
import pyranges as pr
import numpy as np
import argparse
import logging
import warnings
import pysam
import pyBigWig

warnings.filterwarnings("ignore")

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates TSS enrichment plot",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bw_file", type=str, default=None)
    parser.add_argument("--tss", type=str, default=None)
    parser.add_argument("-flank", type=int, default=2000)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    grs = pr.read_bed(args.tss)
    # extend 100 bps to both sides
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - args.flank
    grs.End = mid + args.flank

    logging.info(f"Total of {len(grs)} TSS")
    bw = pyBigWig.open(args.bw_file)

    signal = np.empty(shape=(len(grs), args.flank * 2))

    logging.info(f"Reading data")
    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        try:
          _signal = np.array(bw.values(chrom, start, end))
          _signal[np.isnan(_signal)] = 0
        except:
          _signal = np.zeros(args.flank * 2)
          
        signal[i] = _signal

    np.save(f'{args.out_dir}/{args.out_name}.npy', signal)
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

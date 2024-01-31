import os
import pyranges as pr
import numpy as np
import argparse
import logging
import warnings
import pysam

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
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--tss", type=float, default=None)
    parser.add_argument("-flank", type= int, default=1000)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()

def main():
    args = parse_args()

    grs = pr.read_bigwig(args.tss)
    grs.extend(args.flank)
    
    logging.info(f"Total of {len(grs)} TSS")
    bam = pysam.Samfile(args.bam_file, "rb")
    
    signal = np.zeros()
    
    for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
        for read in bam.fetch(reference=chrom, start=start, end=end):
        # cut counts
        if read.is_reverse:
            cut_site = read.reference_end - 5
        else:
            cut_site = read.reference_start + 4

        if start <= cut_site < end:
            signal[cut_site - start] += 1


    logging.info(f"Done!")


if __name__ == "__main__":
    main()

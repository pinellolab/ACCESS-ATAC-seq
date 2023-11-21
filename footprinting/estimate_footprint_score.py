import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging


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
    parser.add_argument(
        "--bw_file",
        type=str,
        default=None,
        help=(
            "BED file containing ChIP-seq peaks for a specific TF. \n" "Default: None"
        ),
    )

    parser.add_argument(
        "--bed_file",
        type=str,
        default=None,
        help=("BED file containing predicted TF binding sites. \n" "Default: None"),
    )
    
    parser.add_argument(
        "--flanking_region", type=int, default=20, 
        help=("Flanking regions used to estimate footprint score. Default: 20")
    )
        
    parser.add_argument(
        "--outdir",
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


def main():
    args = parse_args()
    
    grs = pr.read_bed(args.bed_file)
    grs = grs.sort()
    
    # get left and right flanking regions
    grs_left = grs.copy()
    grs_left.Start = grs_left.Start - args.flanking_region
    grs_left.End = grs_left.End - args.flanking_region

    grs_right = grs.copy()
    grs_right.Start = grs_right.Start + args.flanking_region
    grs_right.End = grs_right.End + args.flanking_region
    
    bw = pyBigWig.open(args.bw_file)
    
    # compute signal in TF binding site center
    signal_center = np.zeros(len(grs))
    signal_left = np.zeros(len(grs))
    signal_right = np.zeros(len(grs))

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        signal_center[i] = np.sum(bw.values(chrom, start, end))
        
    for i, (chrom, start, end) in enumerate(zip(grs_left.Chromosome, grs_left.Start, grs_left.End)):
        signal_left[i] = np.sum(bw.values(chrom, start, end))

    for i, (chrom, start, end) in enumerate(zip(grs_right.Chromosome, grs_right.Start, grs_right.End)):
        signal_right[i] = np.sum(bw.values(chrom, start, end))
        
    # compute footprint score
    grs.score = np.divide((signal_left + signal_right), signal_center + 1)
    out_filename = os.path.join(args.outdir, "{}.bed".format(args.out_name))
    grs.to_bed(out_filename)
    
    
if __name__ == "__main__":
    main()
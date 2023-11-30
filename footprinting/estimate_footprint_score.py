import os
import sys
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import pyfaidx

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
        "--pos_bed_file",
        type=str,
        default=None,
        help=("BED file containing motif predicted binding sites. \n" "Default: None"),
    )
    parser.add_argument(
        "--neg_bed_file",
        type=str,
        default=None,
        help=("BED file containing motif predicted binding sites. \n" "Default: None"),
    )

    parser.add_argument(
        "--extend",
        type=int,
        default=100,
        help=("Flanking regions used to estimate footprint score. Default: 20"),
    )

    parser.add_argument(
        "--peak_file",
        type=str,
        default=100,
        help=("Flanking regions used to estimate footprint score. Default: 20"),
    )
    parser.add_argument(
        "--ref_fasta",
        type=str,
        default=None,
        help=("FASTA file containing reference genome. \n" "Default: None"),
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


def main():
    args = parse_args()

    grs_pos = pr.read_bed(args.pos_bed_file)
    grs_neg = pr.read_bed(args.neg_bed_file)

    grs_pos.Name = grs_pos.Name + '.Pos'
    grs_neg.Name = grs_neg.Name + '.Neg'

    grs = pr.concat([grs_pos, grs_neg])
    grs = grs.sort()

    grs_peak = pr.read_bed(args.peak_file)
    grs = grs.overlap(grs_peak, strandedness=False, invert=False)

    grs = grs.extend(args.extend)
    bw = pyBigWig.open(args.bw_file)

    pyf = pyfaidx.Fasta(args.ref_fasta)
    grs = pr.genomicfeatures.genome_bounds(grs, chromsizes=pyf, clip=True)

    # compute signal in TF binding site center
    score = np.zeros(len(grs))

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        score[i] = np.sum(bw.values(chrom, start, end))

    # compute footprint score
    grs.Score = score
    out_filename = os.path.join(args.out_dir, "{}.bed".format(args.out_name))
    grs.to_bed(out_filename)


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pyranges as pr
import pyBigWig
import logging
import subprocess as sp


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script estimates bias for Ddd1 based on naked DNA library",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bw_coverage",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--bw_edit",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--peak_file",
        type=str,
        default=None,
        help=(
            "BED file containing genomic regions for generating signal. \n"
            "If none, will use the whole genome as input regions. \n"
            "Default: None"
        ),
    )
    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
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

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    bw_coverage = pyBigWig.open(args.bw_coverage)
    bw_edit_counts = pyBigWig.open(args.bw_edit)

    wig_filename = f"{args.out_dir}/{args.out_name}.wig"
    bw_filename = f"{args.out_dir}/{args.out_name}.bw"

    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            coverage = np.array(bw_coverage.values(chrom, start, end))
            edit_counts = np.array(bw_edit_counts.values(chrom, start, end))

            coverage[np.isnan(coverage)] = 0
            edit_counts[np.isnan(edit_counts)] = 0

            # compute edit fraction
            edit_fraction = np.divide(
                edit_counts,
                coverage,
                out=np.zeros_like(edit_counts),
                where=coverage != 0,
            )
            
            # make no NAN and INF values in the results
            assert not np.isnan(edit_fraction).any(), "Find NAN values"
            assert not np.isinf(edit_fraction).any(), "Find INF values"
                      
 #           edit_fraction[np.isnan(edit_fraction)] = 0

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in edit_fraction))
            f.write("\n")

    # convert to bigwig file
    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()

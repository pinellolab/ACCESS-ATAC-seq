import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pyranges as pr
import pysam
import logging


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

from utils import get_chrom_size_from_bam


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates BigWig file from a BAM file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bam_file",
        type=str,
        default=None,
        help=("BAM file containing the aligned reads. \n" "Default: None"),
    )

    parser.add_argument(
        "--bed_file",
        type=str,
        default=None,
        help=(
            "BED file containing genomic regions for generating signal. \n"
            "If none, will use the whole genome as input regions. \n"
            "Default: None"
        ),
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

    parser.add_argument("--forward_shift", type=int, default=4)

    parser.add_argument("--reverse_shift", type=int, default=4)

    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
    )

    return parser.parse_args()


def get_atac_counts(
    chrom: str = None,
    start: int = None,
    end: int = None,
    forward_shift: int = None,
    reverse_shift: int = None,
    bam: pysam.Samfile = None,
) -> np.array:
    """
    Get Tn5 cutting sites from specific genomic region

    Parameters
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    bam : pysam.Samfile
        BAM file
    """

    signal = np.zeros(shape=(end - start))

    for read in bam.fetch(reference=chrom, start=start, end=end):
        # cut counts
        if read.is_reverse:
            cut_site = read.reference_end - reverse_shift
        else:
            cut_site = read.reference_start + forward_shift

        if start <= cut_site < end:
            signal[cut_site - start] += 1

    return signal


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    if args.bed_file:
        logging.info(f"Loading genomic regions from {args.bed_file}")
        grs = pr.read_bed(args.bed_file)
    else:
        logging.info(f"Using whole genome")
        grs = get_chrom_size_from_bam(bam=bam)

    logging.info(f"Total of {len(grs)} regions")

    output_fname = os.path.join(args.out_dir, "{}.wig".format(args.out_name))

    # Open a new bigwig file for writing
    with open(output_fname, "a") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            signal = get_atac_counts(
                chrom=chrom,
                start=start,
                end=end,
                bam=bam,
                forward_shift=args.forward_shift,
                reverse_shift=args.reverse_shift,
            )

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in signal))
            f.write("\n")

    # convert to bigwig file
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))
    os.system(
        " ".join(
            [
                "wigToBigWig",
                output_fname,
                args.chrom_size_file,
                bw_filename,
                "-verbose=0",
            ]
        )
    )
    os.remove(output_fname)


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pysam
import pyranges as pr
import logging
import subprocess as sp

from utils import get_chrom_size_from_bam

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
    parser.add_argument("--bam_file", type=str, default=None, help="input BAM file")
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
    parser.add_argument("--out_type", type=str, default="fraction")
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
        default=None,
        help=("Names for output file. Default: output"),
    )

    return parser.parse_args()


def get_edit_count(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
) -> np.array:
    """
    Get Ddda editing counts for each position in a genomic region

    Parameters
    ----------
    chrom : str
        Chromosome anme
    start : int
        Start position
    end : int
        End position
    bam : pysam.Samfile
        BAM file
    """

    signal_forward = np.zeros(shape=(end - start))
    signal_reverse = np.zeros(shape=(end - start))

    for read in bam.fetch(reference=chrom, start=start, end=end):
        refer_seq = read.get_reference_sequence().upper()
        query_seq = read.query_sequence.upper()

        # we only look at reads with substitution
        if len(refer_seq) != len(query_seq):
            continue

        for i in range(len(refer_seq)):
            edit_site = read.reference_start + i

            # C -> T at forward strand
            if refer_seq[i] == "C" and query_seq[i] == "T":
                if start <= edit_site < end:
                    signal_forward[edit_site - start] += 1

            # C -> T at reverse strand
            elif refer_seq[i] == "G" and query_seq[i] == "A":
                if start <= edit_site < end:
                    signal_reverse[edit_site - start] += 1

    return signal_forward, signal_reverse


def output_count(bam, grs, chrom_size_file, out_dir, out_name):
    wig_filename = os.path.join(out_dir, "{}.wig".format(out_name))
    bw_filename = os.path.join(out_dir, "{}.wig".format(out_name))

    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            signal_forward, signal_reverse = get_edit_count(chrom, start, end, bam)
            signal = signal_forward + signal_reverse

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in signal))
            f.write("\n")

    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, chrom_size_file, bw_filename])
    os.remove(wig_filename)


def output_fraction(bam, grs, chrom_size_file, out_dir, out_name):
    wig_filename = os.path.join(out_dir, "{}.wig".format(out_name))
    bw_filename = os.path.join(out_dir, "{}.wig".format(out_name))

    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            # get coverage
            coverage = bam.count_coverage(chrom, start, end, quality_threshold=0)
            coverage[np.isnan(coverage)] = 0

            # get edit counts
            signal_forward, signal_reverse = get_edit_count(chrom, start, end, bam)
            edit_count = signal_forward + signal_reverse
            edit_count[np.isnan(edit_count)] = 0

            # compute edit fraction
            edit_fraction = np.divide(
                edit_count,
                coverage,
                out=np.zeros_like(edit_count),
                where=coverage != 0,
            )

            # make sure there are no NAN and INF values in the results
            assert not np.isnan(edit_fraction).any(), "Find NAN values"
            assert not np.isinf(edit_fraction).any(), "Find INF values"

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in edit_fraction))
            f.write("\n")

    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, chrom_size_file, bw_filename])
    os.remove(wig_filename)


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")
    if args.peak_file:
        logging.info(f"Loading genomic regions from {args.peak_file}")
        grs = pr.read_bed(args.peak_file)
        grs = grs.merge()
    else:
        logging.info(f"Using whole genome")
        grs = get_chrom_size_from_bam(bam=bam)

    logging.info(f"Total of {len(grs)} regions")

    # if output type is counts, then compute the raw edit count for each position
    if args.out_type == "counts":
        output_count(
            bam=bam,
            grs=grs,
            chrom_size_file=args.chrom_size_file,
            out_dir=args.out_dir,
            out_name=args.out_name,
        )
    elif args.out_type == "fraction":
        output_fraction(
            bam=bam,
            grs=grs,
            chrom_size_file=args.chrom_size_file,
            out_dir=args.out_dir,
            out_name=args.out_name,
        )

    logging.info("Done!")


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging

from encode_lib_common import run_shell_cmd
from encode_lib_genomic import subsample_ta_pe

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script converts BAM to TAGALIGN file",
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
        "--type",
        type=str,
        choices=["atac", "access", "both"],
        default="both",
        help=(
            "How to quantify chromatin accessibility.\n"
            "atac: only use Tn5 cutting sites\n"
            "access: only use Ddda editting sites\n"
            "both: use both Tn5 cutting and Ddda editing sites"
        ),
    )
    parser.add_argument("--mito-chr-name", default="chrM", help="Mito chromosome name.")
    parser.add_argument(
        "--subsample",
        type=int,
        default=0,
        help="Subsample TAGALIGN. This affects all downstream analysis.",
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


def get_atac(read: pysam.AlignedSegment = None) -> str:
    strand = "-" if read.is_reverse else "+"
    res = f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\tN\t1000\t{strand}\n"

    return res


def get_access(read: pysam.AlignedSegment = None) -> str:
    refer_seq = read.get_reference_sequence().upper()
    query_seq = read.query_sequence.upper()

    res = ""
    # we only look at reads with substitution
    if len(refer_seq) != len(query_seq):
        return res
    else:
        strand = "-" if read.is_reverse else "+"
        for i in range(len(refer_seq)):
            start_site = read.reference_start + i
            end_site = start_site + read.query_length

            # C -> T or G -> A
            if (refer_seq[i] == "C" and query_seq[i] == "T") or (
                refer_seq[i] == "G" and query_seq[i] == "A"
            ):
                res += f"{read.reference_name}\t{start_site}\t{end_site}\tN\t1000\t{strand}\n"

        return res


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    logging.info("Converting BAM to TAGALIGN...")
    out_filename = os.path.join(args.out_dir, f"{args.out_name}.unsort.tagAlign")

    if args.type == "atac":
        logging.info("Using ATAC")
        with open(out_filename, "w") as f:
            for read in bam.fetch():
                res = get_atac(read)
                f.write(res)

    elif args.type == "access":
        logging.info("Using ACCESS")
        with open(out_filename, "w") as f:
            for read in bam.fetch():
                res = get_access(read=read)
                f.write(res)

    elif args.type == "both":
        logging.info("Using both ACCESS and ATAC")
        for read in bam.fetch():
            for read in bam.fetch():
                res_atac = get_atac(read)
                res_access = get_access(read=read)
                f.write(res_atac)
                f.write(res_access)

    else:
        print(f"Unkown data type {args.type}")

    # sort the output by position
    ta = os.path.join(args.out_dir, f"{args.out_name}.tagAlign.gz")
    cmd = f"sort -k1,1 -k2,2n {out_filename} | gzip -nc > {ta}"
    run_shell_cmd(cmd)

    if args.subsample:
        logging.info("Subsampling TAGALIGN...")
        subsample_ta_pe(
            ta, args.subsample, False, args.mito_chr_name, False, args.out_dir
        )

    os.remove(out_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()

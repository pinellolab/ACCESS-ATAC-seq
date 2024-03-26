import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging
import subprocess as sp


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

    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    logging.info("Converting BAM to TAGALIGN...")
    out_filename = os.path.join(args.out_dir, 
                                f"{args.out_name}.unsort.tagAlign")
    
    if args.type == "atac":
        logging.info("Using ATAC")
        with open(out_filename, 'w') as f:
            for read in bam.fetch():
                if read.is_reverse:
                    strand = "-"
                else:
                    strand = "+"

                f.write(
                    f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\tN\t1000\t{strand}\n"
                )


    # sort the output by position
    outfile = os.path.join(args.out_dir, 
                           f"{args.out_name}.tagAlign")
    sp.run(['sort', '-k1,1', '-k2,2n', out_filename, '-o', outfile])
    sp.run(['gzip', '-f', outfile])
    os.remove(out_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()

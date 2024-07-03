import pysam
import argparse
import logging
import warnings

warnings.filterwarnings("ignore")

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script filters reads by alignment score",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    infile = pysam.AlignmentFile(args.bam_file, "rb")
    outfile = pysam.AlignmentFile(
        f"{args.out_dir}/{args.out_name}.bam", "wb", template=infile
    )

    logging.info("Filtering reads by alignment score")
    iter = infile.fetch(until_eof=True)
    for read in iter:

        try:
            align_score = read.get_tag('AS')
        except KeyError:
            align_score = -10000

        # only output best alignment score
        if align_score == 0:
            outfile.write(read)

    infile.close()
    outfile.close()
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

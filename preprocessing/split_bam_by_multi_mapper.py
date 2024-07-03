from collections import defaultdict
import pysam
import argparse
import logging
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script split the input bam file as single and multi mapper",
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
    unique_map_outfile = pysam.AlignmentFile(
        f"{args.out_dir}/{args.out_name}.unique.bam", "wb", template=infile
    )
    multi_map_outfile = pysam.AlignmentFile(
        f"{args.out_dir}/{args.out_name}.multi.bam", "wb", template=infile
    )

    logging.info("Splitting bam file")

    # Dictionary to store read names and their counts
    read_name_counts = defaultdict(int)

    # First pass: Count occurrences of each read name
    for read in infile:
        read_name_counts[read.query_name] += 1

    # Reset file pointer to the beginning
    infile.reset()
    for read in infile:
        if read_name_counts[read.query_name] > 2:
            multi_map_outfile.write(read)
        elif read_name_counts[read.query_name] == 2:
            unique_map_outfile.write(read)

    infile.close()
    unique_map_outfile.close()
    multi_map_outfile.close()
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

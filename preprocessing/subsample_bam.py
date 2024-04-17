import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import logging
import random

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

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
    )

    parser.add_argument("--n_paired_reads", type=int, default=50)
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--out_name",
        type=str,
        default=None,
    )

    return parser.parse_args()


def main():
    args = parse_args()

    logging.info(f"Subsampling")

    bamfile = pysam.AlignmentFile(args.bam_file, "rb")
    total_pairs = sum(1 for read in bamfile if read.is_proper_pair and read.is_read1)
    bamfile.close()
    
    if total_pairs < args.n_paired_reads:
        raise ValueError(
            f"Requested number of paired reads ({args.n_paired_reads}) is more than available paired reads."
        )

    # Generate a random set of read indices to keep
    random.seed(42)
    selected_indices = set(random.sample(range(total_pairs), args.n_paired_reads))

    bamfile = pysam.AlignmentFile(args.bam_file, "rb")
    output_bam = os.path.join(args.out_dir, f"{args.out_name}.bam")
    with pysam.AlignmentFile(output_bam, "wb", template=bamfile) as outbam:
        pair_index = 0
        for read in bamfile:
            if read.is_proper_pair and read.is_read1:
                if pair_index in selected_indices:
                    # Fetch the mate of the current read
                    mate = bamfile.mate(read)
                    outbam.write(read)
                    outbam.write(mate)
                
                pair_index += 1
    
    bamfile.close()
    pysam.index(output_bam)
    logging.info("Done!")


if __name__ == "__main__":
    main()

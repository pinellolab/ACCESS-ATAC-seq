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
        description="This script filters multi-mapped reads",
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

    # Dictionary to store the best alignment for each read pair
    best_alignments = defaultdict(
        lambda: {'read': None, 'score': -float('inf')})

    for read in tqdm(infile):
        qname = read.query_name

        try:
            align_score = read.get_tag('AS')
        except KeyError:
            align_score = -10000

        if best_alignments[qname]['read'] is None:
            best_alignments[qname]['read'] = read
            best_alignments[qname]['score'] = align_score

        else:
            # update the best alignment
            if align_score > best_alignments[qname]['score']:
                best_alignments[qname]['read'] = read
                best_alignments[qname]['score'] = align_score

    # Write the best alignments to the output BAM file
    for alignment in best_alignments.values():
        outfile.write(alignment['read'])

    infile.close()
    outfile.close()
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

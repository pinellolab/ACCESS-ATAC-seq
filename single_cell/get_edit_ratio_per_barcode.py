import pysam
import pandas as pd
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
        description="This script adds barcode to bam file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--bc_tag", type=str, default="CB")
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam_file, "rb")

    barcode_dict = dict()
    logging.info("Adding barcode to bam file")
    iter = bam.fetch(until_eof=True)
    for read in iter:
        # check if there are any edit events
        refer_seq = read.get_reference_sequence().upper()
        query_seq = read.query_sequence.upper()

        # we only look at reads with substitution
        if len(refer_seq) != len(query_seq):
            continue

        barcode = read.get_tag(args.bc_tag)

        if barcode not in barcode_dict:
            barcode_dict[barcode] = {"n_reads": 0, "n_edited_reads": 0}

        # no editting
        if refer_seq == query_seq:
            barcode_dict[barcode]["n_reads"] += 1
        else:
            for i in range(len(refer_seq)):
                if refer_seq[i] == "C" and query_seq[i] == "T":
                    barcode_dict[barcode]["n_edited_reads"] += 1
                elif refer_seq[i] == "G" and query_seq[i] == "A":
                    barcode_dict[barcode]["n_edited_reads"] += 1

    bam.close()

    # output
    df = pd.DataFrame(data=barcode_dict).transpose()
    df.to_csv(f'{args.out_dir}/{args.out_name}.csv')
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

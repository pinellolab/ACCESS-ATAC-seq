import logging
import pyranges as pr
import pysam
import argparse
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")


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
    parser.add_argument("--bam_file", type=str, default=None,)
    parser.add_argument("--bed_file", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)

    return parser.parse_args()


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    logging.info(f"Loading genomic regions from {args.bed_file}")
    df = pd.read_csv(args.bed_file)

    # only output top 10
    df = df.head(10)
    
    print(df)
    
    grs = pr.PyRanges(df)
    print(grs)
    
    logging.info(f"Total of {len(grs)} regions")

    data = {}
    for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
        signal = list()
        
        for read in bam.fetch(reference=chrom, start=start, end=end):
            _signal = np.zeros(shape=(end - start))

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
                        _signal[edit_site - start] = 1

                # C -> T at reverse strand
                elif refer_seq[i] == "G" and query_seq[i] == "A":
                    if start <= edit_site < end:
                        _signal[edit_site - start] += 1
                        
            if _signal.sum() > 0:
                signal.append(_signal)
        
        if len(signal) > 0:
            signal = np.vstack(signal)
            np.savez(f'{args.out_dir}/{chrom}_{start}_{end}.npz', data=signal)
    
    logging.info("Done!")


if __name__ == "__main__":
    main()

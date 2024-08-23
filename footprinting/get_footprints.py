import os
import pyranges as pr
import numpy as np
import argparse
import logging
from scipy.stats import false_discovery_control
import warnings

warnings.filterwarnings("ignore")

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates true labels for TF binding sites based on ChIP-seq peaks",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bw_file", type=str, default=None)
    parser.add_argument("--fdr", type=float, default=0.1)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()

def main():
    args = parse_args()

    grs = pr.read_bigwig(args.bw_file)

    logging.info(f"Filtering footprints by FDR: {args.fdr}")

    # grs.pvalue = 10 ** (-grs.Value)
    
    # grs.fdr = false_discovery_control(grs.pvalue, method="bh")
    grs = grs[grs.Value > -np.log10(args.fdr)]
    # grs = grs.extend(5)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} footprints")
    grs.to_bed(os.path.join(args.out_dir, f"{args.out_name}.bed"))

    logging.info(f"Done!")


if __name__ == "__main__":
    main()

import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pyranges as pr
import pysam
import logging
import pandas as pd
import numpy as np
import pyBigWig
import logomaker
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns

# plt.rcParams['pdf.fonttype'] = 42

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

def get_signal(bw, grs):
    window_size = grs.End.values[0] - grs.Start.values[0]
    signal = np.zeros(shape=(len(grs), window_size))

    for i, (chrom, start, end) in enumerate(zip(grs.Chromosome, grs.Start, grs.End)):
        signal[i] = bw.values(chrom, start, end)

    signal[np.isnan(signal)] = 0
    signal = np.mean(signal, axis=0)

    return signal

def get_all_signal(grs):
    # extend regions
    mid = (grs.End + grs.Start) // 2
    grs.Start = mid - 100
    grs.End = mid + 100

    bw = pyBigWig.open('/data/pinello/PROJECTS/2023_10_ACCESS/results/44_tfbs_pred/03_bam_to_bw/K562_ATAC.bw')
    signal_atac = get_signal(bw=bw, grs=grs)
    
    bw = pyBigWig.open('/data/pinello/PROJECTS/2023_10_ACCESS/results/44_tfbs_pred/03_bam_to_bw/K562_ACCESS_ATAC_10M_cutting_counts.bw')
    signal_access_atac_tn5 = get_signal(bw=bw, grs=grs)
    
    bw = pyBigWig.open('/data/pinello/PROJECTS/2023_10_ACCESS/results/44_tfbs_pred/03_bam_to_bw/K562_ACCESS_ATAC_10M_editing_counts.bw')
    signal_access_atac_ddd1 = get_signal(bw=bw, grs=grs)
    
    df1 = pd.DataFrame(data={"position":range(-100, 100), "signal": signal_atac, "data": "atac"})
    df2 = pd.DataFrame(data={"position":range(-100, 100), "signal": signal_access_atac_tn5, "data": "access_atac_tn5"})
    df3 = pd.DataFrame(data={"position":range(-100, 100), "signal": signal_access_atac_ddd1, "data": "access_atac_ddd1"})
    
    df = pd.concat([df1, df2, df3])
    
    return df

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates BigWig file from a BAM file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "--bed_file",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--name",
        type=str,
        default="counts",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    logging.info(f"Reading regions from {args.bed_file}")
    grs = pr.read_bed(args.bed_file)

    df = get_all_signal(grs)
    
    fig, ax = plt.subplots(1, 1, figsize=[6, 4])
    sns.lineplot(data=df, x="position", y="signal", hue="data")
    plt.legend(loc='upper right')
    plt.title(args.name)

    fig.tight_layout()
    plt.savefig(f'{args.outdir}/{args.name}.png')
    df.to_csv(f"{args.outdir}/{args.name}.csv", index=False)

    logging.info("Done")


if __name__ == "__main__":
    main()

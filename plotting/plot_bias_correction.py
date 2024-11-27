import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import pandas as pd
import argparse
import pyranges as pr
import pysam
import pyBigWig
import seaborn as sns
import logging
import matplotlib.pyplot as plt
import logomaker

from utils import revcomp, get_motif_df

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates footprint plot using BW file and BED file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--obs", type=str, default=None)
    parser.add_argument("--bias", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default="counts")

    return parser.parse_args()


def main():
    args = parse_args()

    df_bias = pd.read_csv(args.bias)
    df_obs = pd.read_csv(args.obs)

    # normalize bias signal
    signal_raw = df_obs['bias'].values
    signal_bias = df_bias['bias'].values

    if np.sum(signal_bias) > 0:
        signal_bias = signal_bias / np.sum(signal_bias)
        
    signal_exp = np.sum(signal_raw) * signal_bias
    signal_norm = np.divide(signal_raw + 0.01, signal_exp + 0.01)
    
    plt.close("all")
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    
    sns.lineplot(signal_bias, ax=axes[0, 0])
    sns.lineplot(signal_raw, ax=axes[0, 1])
    sns.lineplot(signal_exp, ax=axes[1, 0])
    sns.lineplot(signal_norm, ax=axes[1, 1])

    axes[0, 0].set_title("Bias signal")
    axes[0, 1].set_title("Raw signal")
    axes[1, 0].set_title("Exp signal")
    axes[1, 1].set_title("Corrected signal")
    
    fig.tight_layout()

    output_filename = os.path.join(args.out_dir, "{}.png".format(args.out_name))
    plt.savefig(output_filename)
    
    df1 = pd.DataFrame(data={"signal": signal_bias, "data": "signal_bias", "position":range(-50, 50)})
    df2 = pd.DataFrame(data={"signal": signal_raw, "data": "signal_raw", "position":range(-50, 50)})
    df3 = pd.DataFrame(data={"signal": signal_exp, "data": "signal_exp", "position":range(-50, 50)})
    df4 = pd.DataFrame(data={"signal": signal_norm, "data": "signal_norm", "position":range(-50, 50)})
    df = pd.concat([df1, df2, df3, df4])
    df.to_csv(f"{args.out_dir}/{args.out_name}.csv", index=False)


if __name__ == "__main__":
    main()

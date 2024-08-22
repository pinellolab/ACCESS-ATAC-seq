import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import subprocess as sp
from scipy.stats import norm, false_discovery_control
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
    parser.add_argument("--bw_obs_file", type=str, default=None)
    parser.add_argument("--bw_exp_file", type=str, default=None)
    parser.add_argument("--peak_file", type=str, default=None)
    parser.add_argument("--fdr", type=float, default=0.01)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    parser.add_argument("--chrom_size_file", type=str, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    bw_obs = pyBigWig.open(args.bw_obs_file)
    bw_exp = pyBigWig.open(args.bw_exp_file)

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    # logging.info(f"Building background model for footprint score")
    # fp_exp_list = []
    # for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
    #     fp_exp = np.array(bw_exp.values(chrom, start, end))
    #     fp_exp_list.append(fp_exp)

    # fp_exp = np.concatenate(fp_exp_list)
    # fp_exp[np.isnan(fp_exp)] = 0
    # mu, std = norm.fit(fp_exp)

    wig_filename = os.path.join(args.out_dir, "{}.pvalue.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.pvalue.bw".format(args.out_name))

    logging.info(f"Estimating p-values per-nucleotide")
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            fp_obs = np.array(bw_obs.values(chrom, start, end))
            fp_obs[np.isnan(fp_obs)] = 0
            
            fp_exp = np.array(bw_exp.values(chrom, start, end))
            fp_exp[np.isnan(fp_exp)] = 0
            mu, std = norm.fit(fp_exp)
            
            p_values = norm.sf(fp_obs, loc=mu, scale=std)

            masked_array = np.ma.masked_equal(p_values, 0)
            min_pvalue = masked_array.min()
            p_values[p_values == 0] = min_pvalue
            p_values = -np.log10(p_values)

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in p_values))
            f.write("\n")

    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

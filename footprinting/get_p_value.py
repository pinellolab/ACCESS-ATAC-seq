import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import subprocess as sp
from scipy.stats import norm

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
    parser.add_argument("-q", "--qvalue", type=float, default=0.05)
    parser.add_argument("-p", "--pvalue", type=float, default=0.05)
    parser.add_argument("--peak_file", type=str, default=None)
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help=(
            "If specified all output files will be written to that directory. \n"
            "Default: the current working directory"
        ),
    )
    parser.add_argument(
        "--out_name",
        type=str,
        default="counts",
        help=("Names for output file. Default: counts"),
    )
    parser.add_argument(
        "--chrom_size_file",
        type=str,
        default=None,
        help="File including chromosome size. Default: None",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    bw_obs = pyBigWig.open(args.bw_obs_file)
    bw_exp = pyBigWig.open(args.bw_exp_file)

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()
    
    

    chrom_list = chrom_list = grs.Chromosome.unique()

    for chrom in chrom_list:
        grs = grs[grs.Chromosome.isin(["chr1"])]

    logging.info("Building background model!")
    fp_exp_list = []
    for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
        fp_exp = np.array(bw_exp.values(chrom, start, end))
        fp_exp_list.append(fp_exp)

    fp_exp = np.concatenate(fp_exp_list)
    
    mu, std = norm.fit(fp_exp)
    
    logging.info(f"Null distribution: mean: {mu}, std: {std}")

    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))
    
    logging.info("Performing test")
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            fp_obs = np.array(bw_obs.values(chrom, start, end))        
            p_values = norm.sf(fp_obs, loc=mu, scale=std)
        
            p_values = -np.log10(p_values)
            p_values[np.isnan(p_values)] = 100
            
            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in p_values))
            f.write("\n")


    # convert to bigwig file
    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])

    # os.remove(wig_filename)
    logging.info(f"Done!")


if __name__ == "__main__":
    main()

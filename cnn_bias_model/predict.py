import argparse
import os
import sys
import numpy as np
import subprocess
import warnings
import torch
import logging
import pyranges as pr
import pysam

from model import BiasNet
from utils import pad_and_split, one_hot_encode

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

if not sys.warnoptions:
    warnings.simplefilter("ignore")


def parse_args():
    parser = argparse.ArgumentParser()

    # Required parameters
    parser.add_argument(
        "--peak_file", type=str, default=None, help="BED file of input genomic regions"
    )
    parser.add_argument(
        "--ref_fasta", type=str, default=None, help="FASTQ file for reference genome"
    )
    parser.add_argument("--model_path", type=str, default=None, help="Model path")
    parser.add_argument(
        "--chrom_size_file", type=str, default=None, help="chrom_size_file"
    )
    parser.add_argument("--out_dir", type=str, default=None, help="Output directory")
    parser.add_argument("--out_name", type=str, default=None, help="Output name")
    return parser.parse_args()


def main():
    args = parse_args()

    logging.info("Loading input files")
    grs = pr.read_bed(args.peak_file)
    fasta = pysam.FastaFile(args.ref_fasta)

    logging.info("Preprocessing input regions")
    # grs = pad_and_split(grs, k=128)

    # Setup model
    logging.info("Creating model")
    model = BiasNet(seq_len=128)
    state_dict = torch.load(args.model_path)
    model.load_state_dict(state_dict["state_dict"])
    device = torch.device("cuda")
    model.to(device)
    model.eval()

    logging.info("Predicting values for input regions")
    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))

    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            # split the regions
            n = (end - start) // 128
            mid = (start + end) // 2
            start_new = int(mid - (n + 1) / 2 * 128)
            end_new = int(mid + (n + 1) / 2 * 128)

            preds = []
            for i in range(n + 1):
                seq = str(
                    fasta.fetch(chrom, start_new + i * 128, start_new + (i + 1) * 128)
                ).upper()
                x = torch.tensor(one_hot_encode(seq=seq))
                x = x.unsqueeze(dim=0)
                pred = model(x.to(device)).detach().cpu().view(-1).numpy()
                preds.append(pred)

            pred = np.concatenate(preds).tolist()
            pred = pred[(start - start_new) : -(end_new - end)]

            # convert from log-scale to original scale
            pred = np.exp(pred) - 0.01

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in pred))
            f.write("\n")

    logging.info(f"Predicting finished")

    subprocess.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)

if __name__ == "__main__":
    main()

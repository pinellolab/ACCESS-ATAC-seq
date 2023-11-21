import warnings

warnings.filterwarnings("ignore")

import os
import numpy as np
import argparse
import pandas as pd
import pyranges as pr
import pysam
import logging
import logomaker
import matplotlib.pyplot as plt
from tqdm import tqdm


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

# from . import utils
from utils import get_chrom_size_from_bam, revcomp


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
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--ref_fasta",
        type=str,
        default=None,
        help=("BAM file containing reads. \n" "Default: None"),
    )
    parser.add_argument(
        "--bed_file",
        type=str,
        default=None,
        help=("BED file containing open regions. \n" "Default: None"),
    )
    parser.add_argument(
        "--window_size", type=int, default=6, help=("Extend the regions. Default: 100")
    )
    parser.add_argument(
        "--acc",
        type=str,
        choices=["atac", "access"],
        default="atac",
        help=(
            "How to quantify chromatin accessibility.\n"
            "atac: only use Tn5 cutting sites\n"
            "access: only use Ddda editting sites\n"
            "both: use both Tn5 cutting and Ddda editing sites"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help=(
            "If specified all output files will be written to that directory. \n"
            "Default: the current working directory"
        ),
    )
    parser.add_argument(
        "--name",
        type=str,
        default="counts",
        help=("Names for output file. Default: counts"),
    )
    parser.add_argument(
        "--out_format",
        type=str,
        default="png",
        help=("Names for output file. Default: png"),
    )


    return parser.parse_args()


def get_atac_pwm(bam, grs, fasta, window_size):
    pwm = dict(
        [
            ("A", [0.0] * window_size),
            ("C", [0.0] * window_size),
            ("G", [0.0] * window_size),
            ("T", [0.0] * window_size),
            ("N", [0.0] * window_size),
        ]
    )

    for chrom, start, end in tqdm(zip(grs.Chromosome, grs.Start, grs.End)):
        for read in bam.fetch(reference=chrom, start=start, end=end):
            if read.is_reverse:
                cut_size = read.reference_end - 4
            else:
                cut_size = read.reference_start + 4

            p1 = cut_size - window_size // 2
            p2 = p1 + window_size

            try:
                seq = str(fasta.fetch(chrom, p1, p2)).upper()
            except:
                continue

            if read.is_reverse:
                seq = revcomp(seq)

            for i in range(len(seq)):
                pwm[seq[i]][i] += 1

    return pwm


def get_access_pwm(bam, grs, fasta, window_size):
    pwm = dict(
        [
            ("A", [0.0] * window_size),
            ("C", [0.0] * window_size),
            ("G", [0.0] * window_size),
            ("T", [0.0] * window_size),
            ("N", [0.0] * window_size),
        ]
    )

    for chrom, start, end in tqdm(zip(grs.Chromosome, grs.Start, grs.End)):
        for read in bam.fetch(reference=chrom, start=start, end=end):
            
            # get reference and query sequence. 
            # for forward reads, it's from 5' to 3', for reverse reads, it's 3' to 5'
            refer_seq = read.get_reference_sequence().upper()
            query_seq = read.query_sequence.upper()
            
            # we only look at reads with substitution
            if len(refer_seq) != len(query_seq):
                continue
            
            for i in range(len(refer_seq)):
                edit_site = read.reference_start + i
                p1 = edit_site - window_size // 2
                p2 = p1 + window_size
            
                motif_seq = None
                if refer_seq[i] == 'C' and query_seq[i] == 'T':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue

                elif refer_seq[i] == 'G' and query_seq[i] == 'A':
                    try:
                        motif_seq = str(fasta.fetch(chrom, p1, p2)).upper()
                    except:
                        continue
                    motif_seq = revcomp(motif_seq)
                
                if motif_seq:
                    for j in range(len(motif_seq)):
                        if j == (window_size // 2):
                                continue
                        pwm[motif_seq[j]][j] += 1

    return pwm


def get_motif_df(pwm, pseudo_count=1):
    pwm = {k: pwm[k] for k in ("A", "C", "G", "T")}
    pwm = pd.DataFrame(data=pwm)
    pwm = pwm.add(pseudo_count)
    pwm_prob = (pwm.T / pwm.T.sum()).T
    pwm_prob_log = np.log2(pwm_prob)
    pwm_prob_log = pwm_prob_log.mul(pwm_prob)
    info_content = pwm_prob_log.T.sum() + 2
    df = pwm_prob.mul(info_content, axis=0)

    window_size = df.shape[0]

    start = -(window_size // 2)
    end = start + window_size - 1

    df.index = np.linspace(start, end, num=window_size)
    
    return df

def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")
    fasta = pysam.FastaFile(args.ref_fasta)

    if args.bed_file:
        logging.info(f"Reading regions from {args.bed_file}")
        grs = pr.read_bed(args.bed_file)
    else:
        logging.info(f"Using whole genome")
        grs = get_chrom_size_from_bam(bam=bam)

    logging.info(f"Total of {len(grs)} regions")

    if args.outdir is None:
        args.outdir = os.getcwd()

    # get pwm for Tn5
    if args.acc == "atac":
        pwm = get_atac_pwm(bam=bam, grs=grs, fasta=fasta, window_size=args.window_size)
    elif args.acc == 'access':
        pwm = get_access_pwm(bam=bam, grs=grs, fasta=fasta, window_size=args.window_size)

    df = get_motif_df(pwm)

    # plot motif logo
    fig, ax = plt.subplots(1, 1, figsize=[4, 2])
    ax.set_title(args.name)
    ax.set_xlabel('Distance from motif center')
    ax.set_ylabel('Bit score')
    
    logo = logomaker.Logo(df, ax=ax, baseline_width=0)
    
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.xaxis.set_ticks_position('none')
    logo.ax.xaxis.set_tick_params(pad=-1)
    logo.ax.xaxis.set_tick_params(pad=-1)
    
    fig.tight_layout()
    output_fname = os.path.join(args.outdir, "{}.{}".format(args.name, args.out_format))
    plt.savefig(output_fname)


if __name__ == "__main__":
    main()

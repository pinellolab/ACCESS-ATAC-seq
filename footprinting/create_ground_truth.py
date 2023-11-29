import os
import pyranges as pr

import argparse
import logging


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
    parser.add_argument(
        "--chip_seq_peaks",
        type=str,
        default=None,
        help=(
            "BED file containing ChIP-seq peaks for a specific TF. \n" "Default: None"
        ),
    )
    
    parser.add_argument(
        "--motif_match",
        type=str,
        default=None,
        help=("BED file containing predicted TF binding sites. \n" "Default: None"),
    )

    parser.add_argument("--motif_name", type=str, default=None, help=("Motif name"))

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
        "--out_name",
        type=str,
        default="counts",
        help=("Names for output file. Default: counts"),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    chip_peak_grs = pr.read_bed(args.chip_seq_peaks)
    motif_match_grs = pr.read_bed(args.motif_match)

    # remove chrY
    motif_match_grs = motif_match_grs[motif_match_grs.Chromosome != 'chrY']

    motif_match_grs.Name = [name.split('.')[2] for name in motif_match_grs.Name]

    # subset the moti matching results to select candidate TF binding sites
    motif_match_grs = motif_match_grs[motif_match_grs.Name == args.motif_name]

    # overlap with ChIP-seq peaks to get labels
    motif_match_tp_grs = motif_match_grs.overlap(
        chip_peak_grs, strandedness=False, invert=False
    )
    motif_match_tn_grs = motif_match_grs.overlap(
        chip_peak_grs, strandedness=False, invert=True
    )

    # motif_match_tp_grs.Name = [
    #     args.motif_name + ".Pos." + str(i) for i in range(len(motif_match_tp_grs))
    # ]
    # motif_match_tn_grs.Name = [
    #     args.motif_name + ".Neg." + str(i) for i in range(len(motif_match_tn_grs))
    # ]
    
    # motif_match_grs = pr.concat([motif_match_tp_grs, motif_match_tn_grs])
    
    motif_match_tp_grs.to_bed(os.path.join(args.outdir, "{}.Pos.bed".format(args.out_name)))
    motif_match_tn_grs.to_bed(os.path.join(args.outdir, "{}.Neg.bed".format(args.out_name)))


if __name__ == "__main__":
    main()

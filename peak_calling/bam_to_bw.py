import numpy as np
import pyranges as pr
import logging
import pysam
import argparse
import os
import warnings
import subprocess as sp

warnings.filterwarnings("ignore")


logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script converts BAM to TAGALIGN file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--bed_file", type=str, default=None)
    parser.add_argument("--type", type=str, default='atac')
    parser.add_argument("--bin_size", type=int, default=100)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    parser.add_argument("--chrom_size_file", type=str, default=None)

    return parser.parse_args()


def get_chrom_size_from_bam(bam: pysam.Samfile) -> pr.PyRanges:
    """
    Extract chromsome size from the input bam file

    Parameters
    ----------
    bam : pysam.Samfile
        Input bam file

    Returns
    -------
    pr.PyRanges
        A PyRanges object containing chromosome size. Note this is 0-based
    """

    chromosome = list(bam.references)
    start = [0] * len(chromosome)
    end = list(bam.lengths)
    end = [x - 1 for x in end]

    grs = pr.from_dict({"Chromosome": chromosome, "Start": start, "End": end})

    return grs


def get_access(
    bam: pysam.Samfile = None,
    chrom: str = None,
    start: int = None,
    end: int = None,
) -> np.array:
    """
    Get Ddda editing counts for each position in a genomic region

    Parameters
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    bam : pysam.Samfile
        BAM file
    """

    signal = np.zeros(shape=(end - start))

    for read in bam.fetch(reference=chrom, start=start, end=end):
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
                    signal[edit_site - start] += 1

            # C -> T at reverse strand
            elif refer_seq[i] == "G" and query_seq[i] == "A":
                if start <= edit_site < end:
                    signal[edit_site - start] += 1

    return signal


def get_atac(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
) -> np.array:
    """
    Get Tn5 cutting count from specific genomic region

    Parameters
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    bam : pysam.Samfile
        BAM file
    """

    signal = np.zeros(shape=(end - start))

    for read in bam.fetch(reference=chrom, start=start, end=end):
        # cut counts
        if read.is_reverse:
            cut_site = read.reference_end - 5
        else:
            cut_site = read.reference_start + 4

        if start <= cut_site < end:
            signal[cut_site - start] += 1

    return signal


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")
    if args.bed_file:
        logging.info(f"Loading genomic regions from {args.bed_file}")
        grs = pr.read_bed(args.bed_file)
        grs = grs.merge()
    else:
        logging.info(f"Using whole genome")
        grs = get_chrom_size_from_bam(bam=bam)

    logging.info("Converting BAM to WIG...")
    wig_filename = os.path.join(args.out_dir, f"{args.out_name}.wig")
    bw_filename = os.path.join(args.out_dir, f"{args.out_name}.bw")

    bs = args.bin_size

    # Open a new bigwig file for writing
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            if (end - start) % bs == 0:
                n_batchs = (end - start) // bs
            else:
                n_batchs = (end - start) // bs + 1

            for i in range(n_batchs):
                _start = start + i * bs
                _end = min(start + (i + 1) * bs, end)

                if args.type == 'atac':
                    signal = get_atac(
                        chrom=chrom, start=_start, end=_end, bam=bam)

                elif args.type == 'access':
                    signal = get_access(
                        chrom=chrom, start=_start, end=_end, bam=bam)

                elif args.type == 'both':
                    atac_signal = get_atac(
                        chrom=chrom, start=_start, end=_end, bam=bam)
                    access_signal = get_access(
                        chrom=chrom, start=_start, end=_end, bam=bam)

                    signal = atac_signal + access_signal

                _signal = np.empty(_end - _start)
                _signal.fill(np.sum(signal))

                f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                f.write("\n".join(str(e) for e in _signal))
                f.write("\n")

    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()

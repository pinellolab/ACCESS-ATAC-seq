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
        description="This script converts BAM to BEDGRAPH file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)

    return parser.parse_args()


def get_signal(
    chrom: str = None,
    start: int = None,
    end: int = None,
    windon_size: int = None,
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
            half_window_size = int(windon_size / 2)

            _start = max(0, cut_site - start - half_window_size)
            _end = min(cut_site - start + half_window_size, end)

            signal[_start:_end] += 1

    return signal


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    # cut counts
    if read.is_reverse:
        cut_site = read.reference_end - 5
    else:
        cut_site = read.reference_start + 4

    logging.info("Converting BAM to WIG...")
    out_filname = os.path.join(args.out_dir, f"{args.out_name}.bg")

    with open(out_filname, 'w') as f:
        for read in bam:
            if read.is_reverse:
                start_site = read.reference_end
            else:
                start_site = read.reference_start
                
        f.write(f'{read.reference_name}\t{start_site}')
            
            
        

    with open(f'args')

    # Open a new bigwig file for writing
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            # get signal using a batch size of 10kb
            if (end - start) % 10000 == 0:
                n_batchs = (end - start) // 10000
            else:
                n_batchs = (end - start) // 10000 + 1
                
            for i in range(n_batchs):
                _start = start + i * 10000
                _end = min(start + (i + 1) * 10000, end)

                signal = get_signal(chrom=chrom,
                                    start=_start,
                                    end=_end,
                                    windon_size=args.windon_size,
                                    bam=bam)

                f.write(f"fixedStep chrom={chrom} start={_start+1} step=1\n")
                f.write("\n".join(str(e) for e in signal))
                f.write("\n")

    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()

import warnings
warnings.filterwarnings('ignore')

import os
import numpy as np
import argparse
import pyranges as pr
import pysam
import logging


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script generates BigWig file from a BAM file',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('--bw_file', 
                        type=str, 
                        default=None,
                        help=('BIGWIG file containing the signal. \n'
                              'Default: None'))
    
    parser.add_argument('--bed_file', 
                        type=str, 
                        default=None,
                        help=('BED file containing genomic regions for plotting. \n'
                             'Default: None'))
    
    parser.add_argument('--extend',
                        type=int,
                        default=100,
                        help=('Extend the regions. Default: 100'))
    
    
    parser.add_argument('--outdir', 
                        type=str, 
                        default=None,
                        help=('If specified all output files will be written to that directory. \n'
                              'Default: the current working directory'))
        
    parser.add_argument('--name',
                        type=str,
                        default='counts',
                        help=('Names for output file. Default: counts'))
    
    return parser.parse_args()


def get_chrom_size(bam: pysam.Samfile) -> pr.PyRanges:
    """
    Extract chromsome size from the input bam file

    Parameters
    ----------
    bam : pysam.Samfile
        _description_

    Returns
    -------
    pr.PyRanges
        _description_
    """

    chromosome = list(bam.references)
    start = [1] * len(chromosome)
    end = list(bam.lengths)

    grs = pr.from_dict({"Chromosome": chromosome, "Start": start, "End": end})
    
    return grs


def main():
    args = parse_args()
    
    grs = pr.read_bed(args.bed_file)
    
    bam = pysam.Samfile(args.bam_file, "rb")
    
    if args.bed_file:
        logging.info(f'Loading genomic regions from {args.bed_file}')
        grs = pr.read_bed(args.bed_file)
    else:
        logging.info(f'Using whole genome')
        grs = get_chrom_size(bam=bam)
    
    logging.info(f'Total of {len(grs)} regions')
    
    if args.outdir is None:
        args.outdir = os.getcwd()
    
    logging.info(f'Counting {args.acc} sites')
    output_fname = os.path.join(args.outdir, "{}.wig".format(args.name))
    
    # Open a new bigwig file for writing
    with open(output_fname, "a") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            if args.acc == 'cut':
                signal = get_cut_sites(chrom=chrom, start=start, end=end, bam=bam)
                
            elif args.acc == 'edit':
                signal = get_edit_sites(chrom=chrom, start=start, end=end, bam=bam)
                
            elif args.acc == 'both':
                signal_cut = get_cut_sites(chrom=chrom, start=start, end=end, bam=bam)
                signal_edit = get_edit_sites(chrom=chrom, start=start, end=end, bam=bam)
                signal = signal_cut + signal_edit
                
            f.write(f'fixedStep chrom={chrom} start={start+1} step=1\n')
            f.write('\n'.join(str(e) for e in signal))
            f.write('\n')

    # convert to bigwig file
    bw_filename = os.path.join(args.outdir, "{}.bw".format(args.name))
    os.system(" ".join(["wigToBigWig", 
                        output_fname, 
                        args.chrom_size_file, 
                        bw_filename, 
                        "-verbose=0"]))
    os.remove(output_fname)

if __name__ == '__main__':
    main()

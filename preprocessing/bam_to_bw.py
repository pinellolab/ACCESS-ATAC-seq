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
    parser.add_argument('--bam_file', 
                        type=str, 
                        default=None,
                        help=('BAM file containing the aligned reads. \n'
                              'Default: None'))
    
    parser.add_argument('--bed_file', 
                        type=str, 
                        default=None,
                        help=('BED file containing genomic regions for generating signal. \n'
                             'If none, will use the whole genome as input regions. \n'
                             'Default: None'))
    
    parser.add_argument('--outdir', 
                        type=str, 
                        default=None,
                        help=('If specified all output files will be written to that directory. \n'
                              'Default: the current working directory'))
        
    parser.add_argument('--name',
                        type=str,
                        default='counts',
                        help=('Names for output file. Default: counts'))
    
    parser.add_argument('--acc',
                        type=str,
                        choices=['cut', 'edit', 'both'],
                        default='both',
                        help=('How to quantify chromatin accessibility.\n'
                              'cut: only use Tn5 cutting sites\n'
                              'edit: only use Ddda editting sites\n'
                              'both: use both Tn5 cutting and Ddda editing sites'))
    
    parser.add_argument('--chrom_size_file',
                        type=str,
                        default=None,
                        help="File including chromosome size. Default: None")     
        
    return parser.parse_args()


def get_cut_sites(chrom: str = None, 
                  start: int = None, 
                  end: int = None, 
                  bam: pysam.Samfile = None) -> np.array:
    """
    Get Tn5 cutting sites from specific genomic region

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
        if read.is_reverse:
            cut_site = read.reference_end - 5
        else:
            cut_site = read.reference_start + 4
                    
        if start <= cut_site < end:
            signal[cut_site - start] += 1
        
    return signal


def get_edit_sites(chrom: str = None, 
                  start: int = None, 
                  end: int = None, 
                  bam: pysam.Samfile = None) -> np.array:
    """
    Get Ddda editing sites from specific genomic region

    Parameters
    ----------
    chrom : str
        Chromosome anme
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
        
        if len(refer_seq) == len(query_seq):
            for i in range(len(query_seq)):
                if (refer_seq[i] == 'C' and query_seq[i] == 'T') or (refer_seq[i] == 'G' and query_seq[i] == 'A'):
                    if start <= read.reference_start + i < end:
                        signal[read.reference_start + i - start] += 1  
    
    return signal
    

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

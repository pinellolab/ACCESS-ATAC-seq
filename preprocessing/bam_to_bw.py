import warnings
warnings.filterwarnings('ignore')

import os
import numpy as np
import argparse
import pyranges as pr
import pysam
import pyBigWig
from tqdm import tqdm
import logging


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def parse_args():
    parser = argparse.ArgumentParser(
        description='help',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('--bed_file', 
                        type=str, 
                        default=None,
                        help='Input bed file containing genomic regions. Default: None')
    
    parser.add_argument('--bam_file', 
                        type=str, 
                        default=None,
                        help='Input bam file containing the aligned reads. Default: None')

    parser.add_argument('--outdir', 
                        type=str, 
                        default=None,
                        help='If specified all output files will be written to that directory. Default: the current working directory')
        
    parser.add_argument('-n', '--name',
                        type=str,
                        default="counts",
                        help="Names for output file. Default: counts")
    
    parser.add_argument('-a', '--acc',
                        type=str,
                        choices=['cut', 'edit', 'both'],
                        default="both",
                        help="How to quantify chromatin accessibility. Default: both")    
    
    
    parser.add_argument('-ref', '--ref_genome', 
                        type=str, 
                        default=None,
                        help='Reference genome. Default: None')
    
    parser.add_argument('-tb', '--to_bigwig', 
                        action='store_true',
                        help=('Convert wig to bigwig file\n'
                              'To do this, wigtobigwig needs to be installed\n'
                              'Default: False'))
    
    parser.add_argument('--chrom_size_file',
                        type=str,
                        default=None,
                        help=('A file containing chromosome size. Used to convert wig to bigwig'))
        
    return parser.parse_args()


def main():
    args = parse_args()
    
    bam = pysam.Samfile(args.bam_file, "rb")
    
    bam_chrom_info = dict(zip(bam.references, bam.lengths))

    logging.info('Reading bed file')    
    grs = pr.read_bed(args.bed_file)
    logging.info(f'Total of {len(grs)} regions')
    
    if args.outdir is None:
        args.outdir = os.getcwd()
    
    output_fname = os.path.join(args.outdir, "{}.wig".format(args.name))
    
    logging.info(f'Counting {args.acc} sites')
    with open(output_fname, "a") as f:
        for chrom, start, end in tqdm(zip(grs.Chromosome, grs.Start, grs.End)):            
            signal = np.zeros(shape=(end - start), dtype=np.int16)
            # get cut site
            if args.acc == 'cut':
                # for each regions, fetch the reads and get Tn5 cut site
                # BAM file uses 1-based coordinate system
                for read in bam.fetch(reference=chrom, start=start, end=end):
                    if read.is_reverse:
                        cut_site = read.reference_end - 5
                    else:
                        cut_site = read.reference_start + 4
                    
                    if start <= cut_site < end:
                        signal[cut_site - start] += 1
                        
            elif args.acc == 'edit':       
                # get edit site  
                for read in bam.fetch(reference=chrom, start=start, end=end):
                    refer_seq = read.get_reference_sequence().upper()
                    query_seq = read.query_sequence.upper()
                    
                    if len(refer_seq) == len(query_seq):
                        for i in range(len(query_seq)):
                            if (refer_seq[i] == 'C' and query_seq[i] == 'T') or (refer_seq[i] == 'G' and query_seq[i] == 'A'):
                                if start <= read.reference_start + i < end:
                                    signal[read.reference_start + i - start] += 1     
                
                # for pileupcolumn in bam.pileup(reference=chrom, start=start, end=end): 
                #     for read in pileupcolumn.pileups:
                #         if not read.is_del and not read.is_refskip:
                #             refer_seq = read.alignment.get_reference_sequence().upper()
                #             query_seq = read.alignment.query_sequence.upper()
            
                #             if len(refer_seq) == len(query_seq):
                #                 refer_base = refer_seq[read.query_position]
                #                 query_base = query_seq[read.query_position]
                                
                #                 if (refer_base == 'C' and query_base == 'T') or (refer_base == 'G' and query_base == 'A'):
                #                     if start <= pileupcolumn.pos < end:
                #                         signal[pileupcolumn.pos - start] += 1
                             
            
            #signal = cut_sites + edit_sites
            f.write(f'fixedStep chrom={chrom} start={start+1} step=1\n')
            f.write('\n'.join(str(e) for e in signal))
            f.write('\n')

    if args.to_bigwig:
        bw_filename = os.path.join(args.outdir, "{}.bw".format(args.name))
        os.system(" ".join(["wigToBigWig", 
                            output_fname, 
                            args.chrom_size_file, 
                            bw_filename, "-verbose=0"]))
        os.remove(output_fname)

if __name__ == '__main__':
    main()

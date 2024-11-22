import os
import pyBigWig
import pyranges as pr
import numpy as np

import argparse
import logging
import subprocess as sp
import warnings
from multiprocessing import Pool, cpu_count
from rgt.GenomicRegionSet import GenomicRegionSet

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
    parser.add_argument("--bw_files", type=str, default=None)
    
    parser.add_argument("--mpbs-files", metavar='FILE1,FILE2...', type=str,
                        help='Predicted motif binding sites for each condition.'
                             'Files should be separated with comma.')
    parser.add_argument("--conditions", metavar='STRING', type=str,
                        help='Name for each condition. DEFAULT: condition1,condition2, ...')
    parser.add_argument("--window_size", type=int, default=100)
    parser.add_argument("--nc", type=int, default=10)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    parser.add_argument("--chrom_size_file", type=str, default=None)
    return parser.parse_args()


def get_signal(arguments):
    (bw_file, regions, window_size) = arguments
    
    bw = pyBigWig.open(bw_file)
    signal = np.zeros(shape=(len(regions), window_size))

    for i, region in enumerate(regions):
        mid = (region.final + region.initial) // 2
        p1 = mid - window_size // 2
        p2 = mid + window_size // 2
        
        try:
          _signal = bw.values(region.chrom, p1, p2)
        except:  # noqa: E722
          _signal = np.zeros(window_size)
        
        signal[i] = _signal

    signal[np.isnan(signal)] = 0
    signal = np.mean(signal, axis=0)

    return signal
    
def compute_factors(signals):
    signals = np.sum(signals, axis=2)
    # avoid to get a NA value
    signals_log = np.log(signals + 0.01)
    signals_log = signals_log[:, ~np.isnan(signals_log).any(axis=0)]
    signals_log = signals_log - np.mean(signals_log, axis=0, keepdims=True)
    factors = np.exp(np.median(signals_log, axis=1))

    return factors

def main():
    args = parse_args()

    bw_files = args.bw_files.strip().split(",")
    mpbs_files = args.mpbs_files.strip().split(",")
    conditions = args.conditions.strip().split(",")
    
    for bw_file in bw_files:
        assert os.path.exists(bw_file), f"Cannot find file {bw_file}"
    
    for mpbs_file in mpbs_files:
        assert os.path.exists(mpbs_file), f"Cannot find file {mpbs_file}"
    
    logging.info("Read MPBS files")
    mpbs = GenomicRegionSet("Motif Predicted Binding Sites of All Conditions")
    for i, mpbs_file in enumerate(mpbs_files):
        mpbs.read(mpbs_file)
    
    mpbs.sort()
    mpbs_name_list = list(set(mpbs.get_names()))
    signals = np.zeros(shape=(len(conditions), len(mpbs_name_list), args.window_size), dtype=np.float32)
    
    for i, condition in enumerate(conditions):
        logging.info(f"Generate signal for {condition}")
        with Pool(processes=args.nc) as pool:
            arguments_list = list()
            
            for mpbs_name in mpbs_name_list:
                mpbs_regions = mpbs.by_names([mpbs_name])
                mpbs_regions.remove_duplicates()
                 
                arguments = (bw_files[i], mpbs_regions, args.window_size)
                arguments_list.append(arguments)
                
            res = pool.map(get_signal, arguments_list)
            signals[i] = np.array(res)
            
    logging.info("Compute normalization factor")
    factors = compute_factors(signals)
    
    print(factors)
    
    logging.info("Done!")


if __name__ == "__main__":
    main()

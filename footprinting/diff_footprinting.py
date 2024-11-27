import os
import pyBigWig
import pyranges as pr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
    signal = np.sum(signal, axis=0)

    return signal
    
def compute_factors(signals):
    signals = np.sum(signals, axis=2)
    # avoid to get a NA value
    signals_log = np.log(signals + 0.01)
    signals_log = signals_log[:, ~np.isnan(signals_log).any(axis=0)]
    signals_log = signals_log - np.mean(signals_log, axis=0, keepdims=True)
    factors = np.exp(np.median(signals_log, axis=1))

    return factors

def output_profiles(mpbs_name_list, signals, conditions, out_dir):
    for i, condition in enumerate(conditions):
        for j, mpbs_name in enumerate(mpbs_name_list):
            df = pd.DataFrame(data={'signal': signals[i][j]})
            df.to_csv(f"{out_dir}/{condition}_{mpbs_name}.csv")
            # with open(output_filename, "w") as f:
            #     f.write("\t".join(map(str, signals[i][j])) + "\n")

def output_lineplot(arguments):
    (mpbs_name, mpbs_num, signals, conditions, out_dir, window_size, colors) = arguments
    mpbs_name = mpbs_name.replace("(", "_").replace(")", "")

    # output signal
    # output_filename = os.path.join(out_dir, "{}.txt".format(mpbs_name))
    # with open(output_filename, "w") as f:
    #     f.write("\t".join(conditions) + "\n")
    #     for i in range(window_size):
    #         res = []
    #         for j, condition in enumerate(conditions):
    #             res.append(signals[j][i])

    #         f.write("\t".join(map(str, res)) + "\n")

    # # to create a motif loge, we only use A, C, G, T
    # pwm = {k: pwm[k] for k in ('A', 'C', 'G', 'T')}
    # pwm = pd.DataFrame(data=pwm)
    # pwm = pwm.add(1)
    # pwm_prob = (pwm.T / pwm.T.sum()).T
    # pwm_prob_log = np.log2(pwm_prob)
    # pwm_prob_log = pwm_prob_log.mul(pwm_prob)
    # info_content = pwm_prob_log.T.sum() + 2
    # icm = pwm_prob.mul(info_content, axis=0)

    start = -(window_size // 2)
    end = (window_size // 2) - 1
    x = np.linspace(start, end, num=window_size)

    plt.close('all')
    fig, ax = plt.subplots()
    for i, condition in enumerate(conditions):
        ax.plot(x, signals[i], color=colors[i], label=condition)

    ax.text(0.15, 0.9, 'n = {}'.format(mpbs_num), 
            verticalalignment='bottom', 
            horizontalalignment='right',
            transform=ax.transAxes, 
            fontweight='bold')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    min_signal = np.min(signals)
    max_signal = np.max(signals)
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels([str(round(min_signal, 3)), 
                        str(round(max_signal, 3))], rotation=90)

    ax.set_title(mpbs_name, fontweight='bold')
    ax.set_xlim(start, end)
    ax.set_ylim([min_signal, max_signal])
    ax.legend(loc="upper right", frameon=False)
    ax.spines['bottom'].set_position(('outward', 70))

    # ax = plt.axes([0.105, 0.085, 0.85, .2])
    # logo = logomaker.Logo(icm, ax=ax, show_spines=False, baseline_width=0)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()

    output_filename = os.path.join(out_dir, "{}.png".format(mpbs_name))
    plt.savefig(output_filename)


def main():
    args = parse_args()

    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf",
              "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
              "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
              "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
              "#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f", "#bf5b17", "#666666"]


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
    # for test
    # mpbs_name_list = mpbs_name_list[:50]
    logging.info(f"Number of TFs: {len(mpbs_name_list)}")
    
    signals = np.zeros(shape=(len(conditions), len(mpbs_name_list), args.window_size), 
                       dtype=np.float32)
    
    motif_len = list()
    motif_num = list()
    motif_pwm = list()
    
    for i, condition in enumerate(conditions):
        logging.info(f"Generate signal for {condition}")
        with Pool(processes=args.nc) as pool:
            arguments_list = list()
            
            for mpbs_name in mpbs_name_list:
                mpbs_regions = mpbs.by_names([mpbs_name])
                mpbs_regions.remove_duplicates()
                 
                arguments = (bw_files[i], mpbs_regions, args.window_size)
                arguments_list.append(arguments)
                
                motif_num.append(len(mpbs_regions))
                
            res = pool.map(get_signal, arguments_list)
            signals[i] = np.array(res)
            
    logging.info("Compute normalization factor")
    factors = compute_factors(signals)
    
    print(factors)
    
    # normalize signals by factor and number of motifs
    for i in range(len(conditions)):
        for j in range(len(mpbs_name_list)):
            signals[i, j, :] = signals[i, j, :] / (factors[i] * motif_num[j])
    
    logging.info("Save signal")
    output_profiles(mpbs_name_list, signals, conditions, args.out_dir)
    
    logging.info("Plot figure")
    for i, mpbs_name in enumerate(mpbs_name_list):
        output_lineplot((mpbs_name, 
                        motif_num[i], 
                        signals[:, i, :], 
                        conditions, 
                        args.out_dir,
                        args.window_size, colors))
    
    logging.info("Done!")


if __name__ == "__main__":
    main()

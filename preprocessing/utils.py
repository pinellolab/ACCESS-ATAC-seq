import os
import numpy as np
import pandas as pd
import pysam
import pyranges as pr


def wig_to_bw(wig_filename, bw_filename, chrom_size_file):
    os.system(
        " ".join(
            [
                "wigToBigWig",
                wig_filename,
                chrom_size_file,
                bw_filename,
                "-verbose=0",
            ]
        )
    )
    os.remove(wig_filename) 


def revcomp(s):
    rev_dict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
    return "".join([rev_dict[e] for e in s[::-1]])


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
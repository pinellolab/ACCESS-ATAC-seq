import warnings

warnings.filterwarnings("ignore")

import pysam
import numpy as np
from utils import revcomp


def get_bias_signal_access(
    fasta: pysam.FastaFile, chrom: str, start: int, end: int, kmer_dict: dict
) -> np.array:
    """
    Generate bias signal using k-mer model.

    Parameters
    ----------
    fasta : pysam.FastaFile
        Reference genome
    chrom : str
        Chromosome name
    start : int
        Start site
    end : int
        End site
    kmer_dict : dict
        Dictionary object including bias for each k-mer

    Returns
    -------
    np.array
        Bias signal
    """
    k_nb = len(list(kmer_dict.keys())[0])
    bias_signal_forward = np.zeros(shape=(end - start))
    bias_signal_reverse = np.zeros(shape=(end - start))
    seq = str(fasta.fetch(chrom, start - (k_nb // 2), end + (k_nb // 2))).upper()

    for i in range(len(bias_signal_forward)):
        kmer_seq = seq[i : i + k_nb]

        if "N" in kmer_seq:
            continue

        bias_signal_forward[i] = kmer_dict[kmer_seq]
        bias_signal_reverse[i] = kmer_dict[revcomp(kmer_seq)]

    signal = bias_signal_forward + bias_signal_forward

    return signal


def get_raw_signal_access(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
) -> np.array:
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

    signal_forward = np.zeros(shape=(end - start))
    signal_reverse = np.zeros(shape=(end - start))

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
                    signal_forward[edit_site - start] += 1

            # C -> T at reverse strand
            elif refer_seq[i] == "G" and query_seq[i] == "A":
                if start <= edit_site < end:
                    signal_reverse[edit_site - start] += 1

    # signal = signal_forward + signal_reverse

    return signal_forward, signal_reverse


def get_raw_signal_atac(
    chrom: str = None,
    start: int = None,
    end: int = None,
    forward_shift: int = None,
    reverse_shift: int = None,
    bam: pysam.Samfile = None,
) -> np.array:
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
        # cut counts
        if read.is_reverse:
            cut_site = read.reference_end - reverse_shift
        else:
            cut_site = read.reference_start + forward_shift

        if start <= cut_site < end:
            signal[cut_site - start] += 1

    return signal


# def get_exp_signal_access(
#     chrom: str = None,
#     start: int = None,
#     end: int = None,
#     bam: pysam.Samfile = None,
#     window_size: int = 11,
#     fasta: pysam.FastaFile = None,
#     kmer_dict: dict = None,
# ):
#     half_window_size = window_size // 2

#     # get raw signal extended by window_size // 2 to both sides
#     raw_access_signal = get_raw_signal_access(
#         chrom=chrom, start=start - half_window_size, end=end + half_window_size, bam=bam
#     )

#     # smooth signal within a window
#     smooth_access_signal = np.zeros(end - start)
#     for i in range(len(smooth_access_signal)):
#         smooth_access_signal[i] = np.sum(raw_access_signal[i:i+window_size])
    
#     # kernel = np.ones(window_size, "d")
#     # smooth_access_signal = np.convolve(raw_access_signal, kernel, mode="valid")

#     # get bias signal
#     bias_signal = get_bias_signal_access(
#         fasta=fasta,
#         chrom=chrom,
#         start=start - half_window_size,
#         end=end + half_window_size,
#         kmer_dict=kmer_dict,
#     )

#     # smooth bias signal
#     smooth_bias_signal = np.zeros(end - start)
#     for i in range(len(smooth_access_signal)):
#         smooth_bias_signal[i] = np.sum(bias_signal[i:i+window_size])
        
#         # normalize bias signal to [0, 1]
#         if smooth_bias_signal[i] != 0:
#             smooth_bias_signal[i] = bias_signal[i+half_window_size] / smooth_bias_signal[i]

#     exp_signal = np.multiply(smooth_bias_signal, smooth_access_signal)

#     return smooth_access_signal


def get_exp_signal_access(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
    fasta: pysam.FastaFile = None,
    window_size: int = 11,
    kmer_dict: dict = None,
):
    # get raw signal extended by window_size // 2 to both sides
    raw_access_signal = get_raw_signal_access(
        chrom=chrom, start=start, end=end, bam=bam
    )

    # get bias signal
    bias_signal = get_bias_signal_access(
        fasta=fasta,
        chrom=chrom,
        start=start,
        end=end,
        kmer_dict=kmer_dict,
    )

    bias_signal = bias_signal / np.sum(bias_signal)

    exp_signal = np.sum(raw_access_signal) * bias_signal

    return exp_signal

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
    signal_forward = np.zeros(shape=(end - start))
    signal_reverse = np.zeros(shape=(end - start))
    seq = str(fasta.fetch(chrom, start - (k_nb // 2), end + (k_nb // 2))).upper()

    for i in range(len(signal_forward)):
        kmer_seq = seq[i : i + k_nb]

        if "N" in kmer_seq:
            continue

        signal_forward[i] = kmer_dict[kmer_seq]
        signal_reverse[i] = kmer_dict[revcomp(kmer_seq)]

    return signal_forward, signal_reverse


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

    return signal_forward, signal_reverse


def get_bc_signal_access(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
    fasta: pysam.FastaFile = None,
    kmer_dict: dict = None,
):
    # get raw signal
    raw_signal_forward, raw_signal_reverse = get_raw_signal_access(
        chrom=chrom, start=start, end=end, bam=bam
    )

    # get bias signal
    bias_signal_forward, bias_signal_reverse = get_bias_signal_access(
        fasta=fasta,
        chrom=chrom,
        start=start,
        end=end,
        kmer_dict=kmer_dict,
    )

    bias_signal_forward = bias_signal_forward / np.sum(bias_signal_forward)
    bias_signal_reverse = bias_signal_reverse / np.sum(bias_signal_reverse)

    exp_signal_forward = np.sum(raw_signal_forward) * bias_signal_forward
    exp_signal_reverse = np.sum(raw_signal_reverse) * bias_signal_reverse

    bc_signal_forward = np.divide(raw_signal_forward + 1, exp_signal_forward + 1)
    bc_signal_reverse = np.divide(raw_signal_reverse + 1, exp_signal_reverse + 1)

    # for A and T, both of the obseved and expected counts are zeors
    # so they will be nan after dividing, here we replace them with ones
    bc_signal_forward[np.isnan(bc_signal_forward)] = 1
    bc_signal_reverse[np.isnan(bc_signal_reverse)] = 1

    bc_signal_forward = np.log2(bc_signal_forward)
    bc_signal_reverse = np.log2(bc_signal_reverse)

    return bc_signal_forward, bc_signal_reverse


def get_exp_signal_access(
    chrom: str = None,
    start: int = None,
    end: int = None,
    bam: pysam.Samfile = None,
    fasta: pysam.FastaFile = None,
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


def get_raw_signal_atac(
    chrom: str = None,
    start: int = None,
    end: int = None,
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
            cut_site = read.reference_end - 5
        else:
            cut_site = read.reference_start + 4

        if start <= cut_site < end:
            signal[cut_site - start] += 1

    return signal

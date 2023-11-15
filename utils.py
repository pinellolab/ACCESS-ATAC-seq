import pysam
import pyranges as pr

def revcomp(s):
    rev_dict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
    return "".join([rev_dict[e] for e in s[::-1]])


def get_chrom_size_from_bam(bam: pysam.Samfile) -> pr.PyRanges:
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
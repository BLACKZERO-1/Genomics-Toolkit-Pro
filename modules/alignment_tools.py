from Bio import pairwise2
from Bio.Seq import Seq

def pairwise_align(seq1, seq2, match_score, mismatch_penalty, gap_penalty, align_type):
    """
    Performs pairwise alignment using BioPython's pairwise2 module.
    Returns the formatted alignment string.
    """
    if align_type == "Global":
        # Global alignment finds the best alignment over the entire length of the sequences.
        alignments = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)
    else: # Local
        # Local alignment finds the best matching subsequence.
        alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)

    if not alignments:
        return "No alignment found."

    # Format the first and best alignment for display
    formatted_alignment = pairwise2.format_alignment(*alignments[0])
    return formatted_alignment
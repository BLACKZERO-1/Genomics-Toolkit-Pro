import streamlit as st
from Bio import pairwise2, SeqIO, AlignIO
from Bio.Seq import Seq
import subprocess
from Bio.Align import AlignInfo
import math

def pairwise_align(seq1, seq2, match_score, mismatch_penalty, gap_penalty, align_type):
    if align_type == "Global":
        alignments = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)
    else:
        alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)
    if not alignments:
        return "No alignment found."
    return pairwise2.format_alignment(*alignments[0])

def multiple_align(sequences, method):
    input_file, output_file = "temp_input.fasta", "temp_output.aln"
    with open(input_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    if method == "ClustalW":
        command = ["clustalw", f"-INFILE={input_file}", f"-OUTFILE={output_file}", "-OUTPUT=CLUSTAL"]
    elif method == "MUSCLE":
        command = ["muscle", "-in", input_file, "-out", output_file]
    else:
        return None
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return AlignIO.read(output_file, "clustal")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        error_message = e.stderr.decode() if hasattr(e, 'stderr') else f"{method} command not found."
        st.error(f"Error running {method}: {error_message}")
        return None

def generate_consensus(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    return summary_align.dumb_consensus()

def calculate_conservation(alignment):
    scores = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i].upper()
        counts, valid_chars = {}, 0
        for char in column:
            if char in "ACGTU":
                counts[char] = counts.get(char, 0) + 1
                valid_chars += 1
        if valid_chars == 0:
            scores.append(0)
            continue
        entropy = -sum((c / valid_chars) * math.log2(c / valid_chars) for c in counts.values())
        scores.append(2.0 - entropy)
    return scores

def identity_matrix(alignment):
    num_sequences = len(alignment)
    matrix = [[0.0] * num_sequences for _ in range(num_sequences)]
    for i in range(num_sequences):
        for j in range(i, num_sequences):
            seq1, seq2 = alignment[i], alignment[j]
            matches, total_positions = 0, 0
            for k in range(len(seq1)):
                if seq1[k] != '-' and seq2[k] != '-':
                    total_positions += 1
                    if seq1[k] == seq2[k]:
                        matches += 1
            identity = (matches / total_positions) * 100 if total_positions > 0 else 0.0
            matrix[i][j] = matrix[j][i] = identity
    return matrix
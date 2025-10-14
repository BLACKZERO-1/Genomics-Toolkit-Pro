import streamlit as st
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
import subprocess
from Bio.Align import AlignInfo
from Bio import AlignIO
import math

def pairwise_align(seq1, seq2, match_score, mismatch_penalty, gap_penalty, align_type):
    """
    Performs pairwise alignment using BioPython's pairwise2 module.
    """
    if align_type == "Global":
        alignments = pairwise2.align.globalms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)
    else: # Local
        alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_penalty)
    
    if not alignments:
        return "No alignment found."
        
    return pairwise2.format_alignment(*alignments[0])

def multiple_align(sequences, method):
    """
    Performs multiple sequence alignment using an external tool.
    """
    input_file = "temp_input.fasta"
    with open(input_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")

    output_file = "temp_output.aln"
    
    if method == "ClustalW":
        command = ["clustalw", "-INFILE=" + input_file, "-OUTFILE=" + output_file, "-OUTPUT=CLUSTAL"]
    elif method == "MUSCLE":
        command = ["muscle", "-in", input_file, "-out", output_file]
    else:
        return None

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return AlignIO.read(output_file, "clustal")
    except subprocess.CalledProcessError as e:
        st.error(f"Error running {method}: {e.stderr.decode()}")
        return None
    except FileNotFoundError:
        st.error(f"Error: {method} command not found. Ensure it's installed and in your system's PATH.")
        return None

def generate_consensus(alignment):
    """Generates a simple consensus sequence from a multiple alignment."""
    summary_align = AlignInfo.SummaryInfo(alignment)
    return summary_align.dumb_consensus()

def calculate_conservation(alignment):
    """
    Manually calculates a conservation score for each column (Shannon Entropy),
    ignoring gaps and ambiguous characters. This is compatible with all Biopython versions.
    """
    scores = []
    num_sequences = len(alignment)
    
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        counts = {}
        
        # Count only valid characters (A, C, G, T, U)
        valid_chars = 0
        for char in column.upper():
            if char in "ACGTU":
                counts[char] = counts.get(char, 0) + 1
                valid_chars += 1
        
        if valid_chars == 0:
            scores.append(0)
            continue
            
        # Calculate Shannon Entropy
        entropy = 0.0
        for char in counts:
            p = counts[char] / valid_chars
            entropy -= p * math.log2(p)
            
        # Information content = Max Entropy - Observed Entropy
        # Max entropy for DNA/RNA is log2(4) = 2
        info_content = 2.0 - entropy
        scores.append(info_content)
        
    return scores

def identity_matrix(alignment):
    """Calculates a pairwise identity matrix from a multiple alignment."""
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
            
            if total_positions > 0:
                identity = (matches / total_positions) * 100
                matrix[i][j] = matrix[j][i] = identity
            else:
                matrix[i][j] = matrix[j][i] = 0.0
                
    return matrix
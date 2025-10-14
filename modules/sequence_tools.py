from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import Restriction
import numpy as np
import re

def calculate_gc_content(seq):
    """Calculates the GC content of a sequence."""
    return GC(seq)

def transcribe(dna_seq):
    """Transcribes a DNA sequence into RNA."""
    dna = Seq(dna_seq)
    return dna.transcribe()

def translate_six_frames(seq):
    """Translates a DNA sequence in all six reading frames."""
    dna = Seq(seq)
    translations = {}
    for i in range(3):
        translations[f"Frame {i+1}"] = dna[i:].translate(to_stop=True)
        rev_comp = dna.reverse_complement()
        translations[f"Frame -{i+1}"] = rev_comp[i:].translate(to_stop=True)
    return translations

def reverse_complement(seq):
    """Calculates the reverse complement of a sequence."""
    dna = Seq(seq)
    return dna.reverse_complement()

def find_orfs(seq, min_prot_len):
    """Finds all Open Reading Frames (ORFs) in a sequence."""
    seq_obj = Seq(seq)
    orfs = []
    for strand, nuc in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(to_stop=False)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_prot_len:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = frame + aa_end * 3 + 3
                    else:
                        start = len(seq_obj) - (frame + aa_end * 3 + 3)
                        end = len(seq_obj) - (frame + aa_start * 3)
                    orfs.append({
                        "Frame": f"{strand * (frame + 1)}",
                        "Start": start,
                        "End": end,
                        "Length (AA)": aa_end - aa_start,
                        "Protein": str(trans[aa_start:aa_end])
                    })
                aa_start = aa_end + 1
    return orfs

def sliding_window_gc(seq, window_size, step_size):
    """Calculates GC content over a sliding window."""
    positions = []
    gc_values = []
    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i+window_size]
        positions.append(i + window_size // 2)
        gc_values.append(GC(window))
    return positions, gc_values

def codon_usage(seq):
    """Calculates the codon usage frequency for a sequence."""
    seq_obj = Seq(seq)
    codon_counts = {}
    for i in range(0, len(seq_obj) - 2, 3):
        codon = str(seq_obj[i:i+3])
        if codon in codon_counts:
            codon_counts[codon] += 1
        else:
            codon_counts[codon] = 1
    return codon_counts

def find_restriction_sites(seq, enzyme):
    """Finds all cut sites for a given restriction enzyme."""
    batch = Restriction.RestrictionBatch([enzyme])
    search_results = batch.search(Seq(seq))
    cut_sites = search_results.get(enzyme, [])
    return cut_sites
def design_primers(seq, target_start, target_end, primer_len=20):
    """Designs basic forward and reverse primers for a target region."""
    seq_obj = Seq(seq)

    # Design forward primer
    forward_primer_seq = seq_obj[target_start : target_start + primer_len]

    # Design reverse primer
    reverse_primer_template = seq_obj[target_end - primer_len : target_end]
    reverse_primer_seq = reverse_primer_template.reverse_complement()

    # Calculate Melting Temperature (Tm) using the basic formula: 4(G+C) + 2(A+T)
    tm_forward = 4 * (forward_primer_seq.count('G') + forward_primer_seq.count('C')) + 2 * (forward_primer_seq.count('A') + forward_primer_seq.count('T'))
    tm_reverse = 4 * (reverse_primer_seq.count('G') + reverse_primer_seq.count('C')) + 2 * (reverse_primer_seq.count('A') + reverse_primer_seq.count('T'))

    primers = {
        "Forward Primer": str(forward_primer_seq),
        "Forward Tm (°C)": tm_forward,
        "Reverse Primer": str(reverse_primer_seq),
        "Reverse Tm (°C)": tm_reverse
    }
    return primers

def find_motif(seq, pattern):
    """Finds all occurrences of a motif/pattern in a sequence using regex."""
    # re.IGNORECASE makes the search case-insensitive
    matches = re.finditer(pattern, seq, re.IGNORECASE)
    # Store the start and end position of each match
    positions = [(match.start(), match.end()) for match in matches]
    return positions
def kmer_analysis(seq, k):
    """Counts the occurrences of each k-mer in a sequence."""
    kmer_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1
    return kmer_counts
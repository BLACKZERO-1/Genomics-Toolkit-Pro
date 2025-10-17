import re
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import Restriction
import numpy as np

def calculate_gc_content(seq):
    return GC(seq)

def transcribe(dna_seq):
    dna = Seq(dna_seq)
    return dna.transcribe()

def translate_six_frames(seq):
    dna = Seq(seq)
    translations = {}
    for i in range(3):
        translations[f"Frame {i+1}"] = dna[i:].translate(to_stop=True)
        rev_comp = dna.reverse_complement()
        translations[f"Frame -{i+1}"] = rev_comp[i:].translate(to_stop=True)
    return translations

def reverse_complement(seq):
    dna = Seq(seq)
    return dna.reverse_complement()

def find_orfs(seq, min_prot_len):
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
                        "Frame": f"{strand * (frame + 1)}", "Start": start, "End": end,
                        "Length (AA)": aa_end - aa_start, "Protein": str(trans[aa_start:aa_end])
                    })
                aa_start = aa_end + 1
    return orfs

def sliding_window_gc(seq, window_size, step_size):
    positions, gc_values = [], []
    for i in range(0, len(seq) - window_size + 1, step_size):
        window = seq[i:i+window_size]
        positions.append(i + window_size // 2)
        gc_values.append(GC(window))
    return positions, gc_values

def codon_usage(seq):
    seq_obj = Seq(seq)
    codon_counts = {}
    for i in range(0, len(seq_obj) - 2, 3):
        codon = str(seq_obj[i:i+3])
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
    return codon_counts

def find_restriction_sites(seq, enzyme):
    batch = Restriction.RestrictionBatch([enzyme])
    search_results = batch.search(Seq(seq))
    return search_results.get(enzyme, [])

def design_primers(seq, target_start, target_end, primer_len=20):
    seq_obj = Seq(seq)
    forward_primer_seq = seq_obj[target_start : target_start + primer_len]
    reverse_primer_template = seq_obj[target_end - primer_len : target_end]
    reverse_primer_seq = reverse_primer_template.reverse_complement()
    tm_forward = 4 * (forward_primer_seq.count('G') + forward_primer_seq.count('C')) + 2 * (forward_primer_seq.count('A') + forward_primer_seq.count('T'))
    tm_reverse = 4 * (reverse_primer_seq.count('G') + reverse_primer_seq.count('C')) + 2 * (reverse_primer_seq.count('A') + reverse_primer_seq.count('T'))
    return {
        "Forward Primer": str(forward_primer_seq), "Forward Tm (°C)": tm_forward,
        "Reverse Primer": str(reverse_primer_seq), "Reverse Tm (°C)": tm_reverse
    }

def find_motif(seq, pattern):
    matches = re.finditer(pattern, seq, re.IGNORECASE)
    return [(match.start(), match.end()) for match in matches]

def kmer_analysis(seq, k):
    kmer_counts = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts
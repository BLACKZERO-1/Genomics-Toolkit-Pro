import streamlit as st
from modules import sequence_tools as stoo
from Bio.Seq import Seq
import pandas as pd
from modules import visualization as viz

st.set_page_config(page_title="Sequence Analysis", layout="wide")

st.title("🔬 Sequence Analysis Tools")

# --- Input Section ---
st.header("1. Input Sequence")

input_sequence = st.text_area("Paste your sequence here (DNA, RNA, or Protein)", height=200, placeholder=">Example Sequence\nACGTACGTACGT...")
uploaded_file = st.file_uploader("Or upload a FASTA file", type=["fasta", "fa"])

if uploaded_file:
    # Read sequence from uploaded file
    try:
        fasta_string = uploaded_file.getvalue().decode("utf-8")
        # Simple parser for the first sequence in a FASTA file
        sequence_lines = fasta_string.splitlines()[1:]
        input_sequence = "".join(sequence_lines)
        st.success("FASTA file uploaded and parsed successfully!")
    except Exception as e:
        st.error(f"Error parsing FASTA file: {e}")

st.divider()

# --- Analysis Section ---
if input_sequence:
    st.header("2. Analysis Results")

    # --- Basic Analysis ---
    st.subheader("Basic Analysis")
    seq_len = len(input_sequence)
    gc_content = stoo.calculate_gc_content(input_sequence)

    col1, col2 = st.columns(2)
    col1.metric("Sequence Length", f"{seq_len} bp")
    col2.metric("GC Content", f"{gc_content:.2f} %")

    st.divider()

    # --- Transcription, Translation, and Complements ---
    st.subheader("Transcription, Translation & Complements")
    tab1, tab2, tab3 = st.tabs(["Transcription (DNA -> RNA)", "6-Frame Translation", "Complements"])

    with tab1:
        st.write("Converts a DNA sequence into its corresponding RNA sequence.")
        if st.button("Transcribe DNA to RNA"):
            rna_seq = stoo.transcribe(input_sequence)
            st.code(rna_seq, language="text")

    with tab2:
        st.write("Translates a DNA sequence into amino acid sequences for all six possible reading frames.")
        if st.button("Translate in 6 Frames"):
            translations = stoo.translate_six_frames(input_sequence)
            df = pd.DataFrame(list(translations.items()), columns=['Frame', 'Amino Acid Sequence'])
            st.dataframe(df, use_container_width=True)

    with tab3:
        st.write("Generates the complement and reverse complement of a DNA sequence.")
        if st.button("Calculate Complements"):
            rev_comp = stoo.reverse_complement(input_sequence)
            st.markdown(f"**Reverse Complement:**")
            st.code(rev_comp, language="text")

    st.divider()

    # --- ORF Finder ---
    st.subheader("Open Reading Frame (ORF) Finder")
    min_len = st.number_input("Minimum Protein Length (Amino Acids)", min_value=10, value=50, step=10)

    if st.button("Find ORFs"):
        orfs = stoo.find_orfs(input_sequence, min_len)
        if orfs:
            st.write(f"Found {len(orfs)} ORFs with a minimum length of {min_len} amino acids.")
            df_orfs = pd.DataFrame(orfs)
            st.dataframe(df_orfs, use_container_width=True)
        else:
            st.warning("No ORFs found with the specified minimum length.")
st.divider()

# --- GC Sliding Window ---
st.subheader("GC Content Sliding Window")
col1, col2 = st.columns(2)
window_size = col1.number_input("Window Size (bp)", min_value=1, value=100, step=10)
step_size = col2.number_input("Step Size (bp)", min_value=1, value=10, step=1)

if st.button("Calculate GC Window"):
    positions, gc_values = stoo.sliding_window_gc(input_sequence, window_size, step_size)
    if positions:
        fig = viz.plot_gc_content(positions, gc_values)
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("Sequence is too short for the selected window size.")
st.divider()

# --- Codon Usage ---
st.subheader("Codon Usage Analysis")
if st.button("Calculate Codon Usage"):
    codon_counts = stoo.codon_usage(input_sequence)
    if codon_counts:
        fig = viz.plot_codon_heatmap(codon_counts)
        st.plotly_chart(fig, use_container_width=True)

        # Show the raw data in an expandable section
        with st.expander("Show Raw Data"):
            st.dataframe(pd.DataFrame(list(codon_counts.items()), columns=['Codon', 'Count']), use_container_width=True)
    else:
        st.warning("Could not calculate codon usage. Ensure the sequence is a valid DNA sequence.")
else:
    st.info("Please paste a sequence or upload a FASTA file to begin analysis.")
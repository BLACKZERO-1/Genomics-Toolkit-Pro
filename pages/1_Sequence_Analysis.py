import streamlit as st
from modules import sequence_tools as stoo
from modules import visualization as viz
from Bio.Seq import Seq
import pandas as pd

st.set_page_config(page_title="Sequence Analysis", layout="wide")

st.title("🔬 Sequence Analysis Tools")

# --- Input Section ---
st.header("1. Input Sequence")

input_sequence = st.text_area("Paste your sequence here (DNA, RNA, or Protein)", height=200, placeholder=">Example Sequence\nACGTACGTACGT...")
uploaded_file = st.file_uploader("Or upload a FASTA file", type=["fasta", "fa"])

# Logic to handle file upload
if uploaded_file:
    try:
        fasta_string = uploaded_file.getvalue().decode("utf-8")
        sequence_lines = fasta_string.splitlines()[1:]
        input_sequence = "".join(sequence_lines)
    except Exception as e:
        st.error(f"Error parsing FASTA file: {e}")

# --- Sequence Sanitization ---
# Clean the sequence to remove non-alphabetic characters before any analysis
cleaned_sequence = ""
if input_sequence:
    cleaned_sequence = ''.join(filter(str.isalpha, input_sequence)).upper()

st.divider()

# --- Analysis Results Header ---
st.header("2. Analysis Results")

# --- Basic Analysis ---
st.subheader("Basic Analysis")
if cleaned_sequence:
    seq_len = len(cleaned_sequence)
    gc_content = stoo.calculate_gc_content(cleaned_sequence)
    col1, col2 = st.columns(2)
    col1.metric("Sequence Length", f"{seq_len} bp")
    col2.metric("GC Content", f"{gc_content:.2f} %")
else:
    st.warning("Provide a sequence to see basic analysis.")
st.divider()

# --- Transcription, Translation, and Complements ---
st.subheader("Transcription, Translation & Complements")
tab1, tab2, tab3 = st.tabs(["Transcription (DNA -> RNA)", "6-Frame Translation", "Complements"])
with tab1:
    if st.button("Transcribe DNA to RNA"):
        if cleaned_sequence:
            rna_seq = stoo.transcribe(cleaned_sequence)
            st.code(rna_seq, language="text")
        else:
            st.error("Please provide a sequence first.")
with tab2:
    if st.button("Translate in 6 Frames"):
        if cleaned_sequence:
            translations = stoo.translate_six_frames(cleaned_sequence)
            df = pd.DataFrame(list(translations.items()), columns=['Frame', 'Amino Acid Sequence'])
            st.dataframe(df, use_container_width=True)
        else:
            st.error("Please provide a sequence first.")
with tab3:
    if st.button("Calculate Complements"):
        if cleaned_sequence:
            rev_comp = stoo.reverse_complement(cleaned_sequence)
            st.markdown(f"**Reverse Complement:**")
            st.code(rev_comp, language="text")
        else:
            st.error("Please provide a sequence first.")
st.divider()

# --- ORF Finder ---
st.subheader("Open Reading Frame (ORF) Finder")
min_len = st.number_input("Minimum Protein Length (Amino Acids)", min_value=10, value=50, step=10)
if st.button("Find ORFs"):
    if cleaned_sequence:
        orfs = stoo.find_orfs(cleaned_sequence, min_len)
        if orfs:
            st.write(f"Found {len(orfs)} ORFs with a minimum length of {min_len} amino acids.")
            df_orfs = pd.DataFrame(orfs)
            st.dataframe(df_orfs, use_container_width=True)
        else:
            st.warning("No ORFs found with the specified minimum length.")
    else:
        st.error("Please provide a sequence first.")
st.divider()

# --- GC Sliding Window ---
st.subheader("GC Content Sliding Window")
col1, col2 = st.columns(2)
window_size = col1.number_input("Window Size (bp)", min_value=1, value=100, step=10)
step_size = col2.number_input("Step Size (bp)", min_value=1, value=10, step=1)
if st.button("Calculate GC Window"):
    if cleaned_sequence:
        positions, gc_values = stoo.sliding_window_gc(cleaned_sequence, window_size, step_size)
        if positions:
            fig = viz.plot_gc_content(positions, gc_values)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("Sequence is too short for the selected window size.")
    else:
        st.error("Please provide a sequence first.")
st.divider()

# --- Codon Usage ---
st.subheader("Codon Usage Analysis")
if st.button("Calculate Codon Usage"):
    if cleaned_sequence:
        codon_counts = stoo.codon_usage(cleaned_sequence)
        if codon_counts:
            fig = viz.plot_codon_heatmap(codon_counts)
            st.plotly_chart(fig, use_container_width=True)
            with st.expander("Show Raw Data"):
                st.dataframe(pd.DataFrame(list(codon_counts.items()), columns=['Codon', 'Count']), use_container_width=True)
        else:
            st.warning("Could not calculate codon usage.")
    else:
        st.error("Please provide a sequence first.")
st.divider()

# --- Restriction Sites ---
st.subheader("Restriction Site Analysis")
enzyme_list = ["EcoRI", "BamHI", "HindIII", "NotI", "SacI", "SpeI", "XbaI"]
selected_enzyme = st.selectbox("Select an Enzyme", options=enzyme_list)
if st.button(f"Find Sites for {selected_enzyme}"):
    if cleaned_sequence:
        cut_sites = stoo.find_restriction_sites(cleaned_sequence, selected_enzyme)
        if cut_sites:
            st.success(f"Found {len(cut_sites)} cut site(s) for {selected_enzyme} at positions: {', '.join(map(str, cut_sites))}")
            seq_display = list(cleaned_sequence)
            for site in reversed(cut_sites):
                seq_display.insert(site, "|")
            st.code("".join(seq_display))
        else:
            st.warning(f"No cut sites found for {selected_enzyme} in the sequence.")
    else:
        st.error("Please provide a sequence first.")
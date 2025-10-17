import streamlit as st
import pandas as pd
from Bio import SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from modules import sequence_tools as stoo
from modules import alignment_tools as atools
from modules import phylo_tools as ptools
from modules import visualization as viz
import io
import streamlit as st
import os
import streamlit as st
# (your other imports)

# --- Main Page Configuration ---
st.set_page_config(
    page_title="Genomics Toolkit Pro",
    page_icon="🧬",
    layout="wide"
)
# --- Load Custom CSS ---
import os # Add this import at the top with your others

def local_css(file_name):
    # Get the absolute path of the directory this script is in
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the full path to the CSS file
    css_file_path = os.path.join(_this_dir, file_name)

    try:
        with open(css_file_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        st.error(f"Could not find CSS file at: {css_file_path}")

# Load the CSS
local_css("style.css")

# (The rest of your app.py code...)
# --- Hide the default "app" page ---
# This part is now empty, ensuring that Streamlit defaults to the first page
# in the 'pages' directory ('0_🏠_Home.py') and does not show an "app" page.

# --- Initialize Session State ---
if 'sequence_input' not in st.session_state: st.session_state.sequence_input = ""
if 'alignment_input' not in st.session_state: st.session_state.alignment_input = ""
if 'phylo_input' not in st.session_state: st.session_state.phylo_input = ""
if 'sequences' not in st.session_state: st.session_state.sequences = []
if 'alignment' not in st.session_state: st.session_state.alignment = None

# --- Sidebar Navigation ---
page = st.sidebar.radio("", ["🏠 Home", "🔬 Sequence Analysis", "⛓️ Sequence Alignment", "🌳 Phylogenetic Analysis"])

# =====================================================================================
# --- HOME PAGE ---
# =====================================================================================
if page == "🏠 Home":
    st.title("🧬 Genomics Toolkit Pro")
    st.markdown("Welcome to your one-stop solution for advanced genomics analysis.")
    st.divider()
    st.header("Toolkit Modules")
    col1, col2, col3 = st.columns(3, gap="large")
    with col1:
        with st.container(border=True):
            st.subheader("🔬 Sequence Analysis")
            st.write("Perform detailed analysis of DNA, RNA, or protein sequences.")
    with col2:
        with st.container(border=True):
            st.subheader("⛓️ Sequence Alignment")
            st.write("Align multiple sequences using ClustalW and MUSCLE.")
    with col3:
        with st.container(border=True):
            st.subheader("🌳 Phylogenetic Analysis")
            st.write("Construct and visualize phylogenetic trees from alignments.")
    st.divider()
    st.header("Project Dashboard")
    d1, d2, d3 = st.columns(3)
    d1.metric("Tools Implemented", "20+", "Sequence, Alignment, Phylo")
    d2.metric("Analyses Performed", "1,245", "7%")
    d3.metric("Sample Datasets", "3", "FASTA format")

# =====================================================================================
# --- SEQUENCE ANALYSIS PAGE ---
# =====================================================================================
elif page == "🔬 Sequence Analysis":
    st.title("🔬 Sequence Analysis Tools")
    st.header("1. Input Sequence")
    st.subheader("Load Sample Data")
    if st.button("Load Sample DNA"):
        try:
            with open("data/sample_dna.fasta", "r") as f:
                st.session_state.sequence_input = "".join(f.read().splitlines()[1:])
            st.success("Sample DNA sequence loaded!")
            st.rerun()
        except FileNotFoundError:
            st.error("Error: sample_dna.fasta not found.")
    st.subheader("Enter Your Data")
    input_sequence = st.text_area("Paste sequence", value=st.session_state.sequence_input, height=150)
    uploaded_file = st.file_uploader("Or upload a FASTA file", type=["fasta", "fa"], key="seq_analysis_uploader")
    if uploaded_file:
        try:
            input_sequence = "".join(uploaded_file.getvalue().decode("utf-8").splitlines()[1:])
            st.session_state.sequence_input = input_sequence
        except Exception as e:
            st.error(f"Error parsing file: {e}")
    cleaned_sequence = ''.join(filter(str.isalpha, input_sequence)).upper().replace("U", "T")
    st.divider()
    st.header("2. Analysis Results")
    # All your sequence analysis tools follow...
    # (This part of your code was correct and is preserved here)
    st.subheader("Basic Analysis")
    if cleaned_sequence:
        col1, col2 = st.columns(2)
        col1.metric("Length", f"{len(cleaned_sequence)} bp")
        col2.metric("GC Content", f"{stoo.calculate_gc_content(cleaned_sequence):.2f} %")
    else: st.warning("Provide a sequence to see basic analysis.")
    st.divider()
    st.subheader("Transcription, Translation & Complements")
    tab1, tab2, tab3 = st.tabs(["Transcription", "Translation", "Complements"])
    with tab1:
        if st.button("Transcribe DNA to RNA"):
            if cleaned_sequence: st.code(stoo.transcribe(cleaned_sequence))
            else: st.error("Provide a sequence first.")
    with tab2:
        if st.button("Translate in 6 Frames"):
            if cleaned_sequence: st.dataframe(pd.DataFrame(list(stoo.translate_six_frames(cleaned_sequence).items()), columns=['Frame', 'Sequence']))
            else: st.error("Provide a sequence first.")
    with tab3:
        if st.button("Calculate Complements"):
            if cleaned_sequence: st.code(stoo.reverse_complement(cleaned_sequence))
            else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("ORF Finder")
    min_len = st.number_input("Min Protein Length (AA)", 10, 100, 50, 10)
    if st.button("Find ORFs"):
        if cleaned_sequence:
            orfs = stoo.find_orfs(cleaned_sequence, min_len)
            if orfs: st.dataframe(pd.DataFrame(orfs))
            else: st.warning("No ORFs found.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("GC Sliding Window")
    win, step = st.columns(2)
    win_size = win.number_input("Window Size (bp)", 1, 200, 100, 10)
    step_size = step.number_input("Step Size (bp)", 1, 100, 10, 1)
    if st.button("Calculate GC Window"):
        if cleaned_sequence:
            pos, gc = stoo.sliding_window_gc(cleaned_sequence, win_size, step_size)
            if pos: st.plotly_chart(viz.plot_gc_content(pos, gc))
            else: st.warning("Sequence too short.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("Codon Usage Analysis")
    if st.button("Calculate Codon Usage"):
        if cleaned_sequence:
            counts = stoo.codon_usage(cleaned_sequence)
            if counts: st.plotly_chart(viz.plot_codon_heatmap(counts))
            else: st.warning("Could not calculate.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("Restriction Site Analysis")
    enzyme = st.selectbox("Enzyme", ["EcoRI", "BamHI", "HindIII", "NotI", "SacI", "SpeI", "XbaI"])
    if st.button(f"Find {enzyme} Sites"):
        if cleaned_sequence:
            sites = stoo.find_restriction_sites(cleaned_sequence, enzyme)
            if sites: st.success(f"Found {len(sites)} site(s) at: {', '.join(map(str, sites))}")
            else: st.warning("No sites found.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("Primer Design")
    s, e = st.columns(2)
    start = s.number_input("Target Start", 0, 500, 100)
    end = e.number_input("Target End", 1, 1000, 300)
    if st.button("Design Primers"):
        if cleaned_sequence:
            if end > start and end <= len(cleaned_sequence): st.dataframe(pd.DataFrame([stoo.design_primers(cleaned_sequence, start, end)]))
            else: st.error("Invalid region.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("Motif Search")
    pattern = st.text_input("Motif (Regex)", "A[ATGC]G")
    if st.button("Search Motif"):
        if cleaned_sequence:
            matches = stoo.find_motif(cleaned_sequence, pattern)
            if matches: st.dataframe(pd.DataFrame(matches, columns=['Start', 'End']))
            else: st.warning("No matches found.")
        else: st.error("Provide a sequence first.")
    st.divider()
    st.subheader("k-mer Analysis")
    k = st.slider("k-mer size (k)", 2, 6, 3)
    if st.button("Analyze k-mers"):
        if cleaned_sequence:
            counts = stoo.kmer_analysis(cleaned_sequence, k)
            if counts: st.plotly_chart(viz.plot_kmer_distribution(counts))
            else: st.warning("Could not perform analysis.")
        else: st.error("Provide a sequence first.")

# =====================================================================================
# --- SEQUENCE ALIGNMENT PAGE ---
# =====================================================================================
elif page == "⛓️ Sequence Alignment":
    st.title("⛓️ Sequence Alignment Tools")
    st.header("1. Input Sequences")
    if st.button("Load Sample Alignment"):
        try:
            with open("data/sample_multi.fasta", "r") as f:
                st.session_state.alignment_input = f.read()
            st.success("Sample alignment data loaded!")
            st.rerun()
        except FileNotFoundError:
            st.error("Error: sample_multi.fasta not found.")
    input_text = st.text_area("Paste multi-FASTA", value=st.session_state.alignment_input, height=250)
    uploaded_file = st.file_uploader("Or upload multi-FASTA file", type=["fasta", "fa"], key="align_uploader")
    fasta_content = input_text or (uploaded_file.getvalue().decode("utf-8") if uploaded_file else None)
    if fasta_content:
        try:
            st.session_state.sequences = list(SeqIO.parse(io.StringIO(fasta_content), "fasta"))
            st.session_state.alignment = None
            st.success(f"Loaded {len(st.session_state.sequences)} sequences.")
        except Exception as e:
            st.error(f"Error parsing FASTA: {e}")
    st.divider()
    st.header("2. Pairwise Alignment")
    if st.session_state.sequences:
        records = {rec.id: rec for rec in st.session_state.sequences}
        ids = list(records.keys())
        id1 = st.selectbox("Seq 1", ids)
        id2 = st.selectbox("Seq 2", ids, index=min(1, len(ids)-1))
    else:
        st.selectbox("Seq 1", ["No sequences loaded"])
        st.selectbox("Seq 2", ["No sequences loaded"])
    st.subheader("Parameters")
    align_type = st.radio("Type", ("Global", "Local"))
    m, mm, g = st.columns(3)
    match = m.number_input("Match", value=1.0)
    mismatch = mm.number_input("Mismatch", value=-1.0)
    gap = g.number_input("Gap", value=-0.5)
    if st.button("Run Pairwise"):
        if st.session_state.sequences:
            res = atools.pairwise_align(str(records[id1].seq), str(records[id2].seq), match, mismatch, gap, align_type)
            st.code(res)
        else: st.error("Upload sequences first.")
    st.divider()
    st.header("3. Multiple Sequence Alignment (MSA)")
    msa_method = st.radio("Tool", ("ClustalW", "MUSCLE"))
    if st.button("Run MSA"):
        if st.session_state.sequences and len(st.session_state.sequences) > 2:
            with st.spinner(f"Running {msa_method}..."):
                st.session_state.alignment = atools.multiple_align(st.session_state.sequences, msa_method)
            if st.session_state.alignment: st.success("MSA complete.")
            else: st.error("MSA failed.")
        else: st.warning("Requires at least 3 sequences.")
    st.divider()
    st.header("4. Alignment Analysis")
    if st.session_state.alignment:
        aln = st.session_state.alignment
        st.text_area("Alignment", str(aln), height=300)
        st.download_button("Download", str(aln), "alignment.aln")
        st.subheader("Consensus")
        st.code(atools.generate_consensus(aln))
        st.subheader("Conservation")
        st.plotly_chart(viz.plot_alignment_conservation(atools.calculate_conservation(aln)))
        st.subheader("Identity Matrix")
        st.plotly_chart(viz.plot_identity_heatmap(atools.identity_matrix(aln), [rec.id for rec in aln]))
    else:
        st.info("Run MSA to see analysis.")

# =====================================================================================
# --- PHYLOGENETICS PAGE ---
# =====================================================================================
elif page == "🌳 Phylogenetic Analysis":
    st.title("🌳 Phylogenetic Analysis")
    st.header("1. Input Alignment")
    if st.button("Load Sample Alignment Data"):
        try:
            with open("data/sample_multi.fasta", "r") as f:
                st.session_state.phylo_input = f.read()
            st.success("Sample alignment data loaded!")
            st.rerun()
        except FileNotFoundError:
            st.error("Error: sample_multi.fasta not found.")
    input_text = st.text_area("Paste alignment", value=st.session_state.phylo_input, height=200)
    uploaded_file = st.file_uploader("Or upload alignment file", type=["fasta", "fa", "aln"], key="phylo_uploader")
    aln_content = input_text or (uploaded_file.getvalue().decode("utf-8") if uploaded_file else None)
    if aln_content:
        stringio = io.StringIO(aln_content)
        try:
            st.session_state.alignment = AlignIO.read(stringio, "clustal")
            st.success("Clustal alignment loaded.")
        except ValueError:
            stringio.seek(0)
            try:
                st.session_state.alignment = AlignIO.read(stringio, "fasta")
                st.success("FASTA alignment loaded.")
            except Exception as e:
                st.error(f"Error parsing alignment: {e}")
    st.divider()
    st.header("2. Tree Construction")
    if 'alignment' in st.session_state and st.session_state.alignment:
        method = st.radio("Method", ["Neighbor-Joining (NJ)", "UPGMA"])
        bootstrap_replicates = st.slider("Bootstrap Replicates", 10, 1000, 100, 10)
        if st.button("Build Tree"):
            with st.spinner("Building tree and running bootstrap... This may be slow."):
                tree_method_short = 'nj' if method == "Neighbor-Joining (NJ)" else 'upgma'
                main_tree = ptools.build_nj_tree(st.session_state.alignment) if tree_method_short == 'nj' else ptools.build_upgma_tree(st.session_state.alignment)
                if isinstance(main_tree, str):
                    st.error(f"Failed: {main_tree}")
                else:
                    tree_with_support = ptools.apply_bootstrap(main_tree, st.session_state.alignment, bootstrap_replicates, tree_method_short)
                    img_path = ptools.draw_tree(tree_with_support)
                    st.image(img_path)
                    st.subheader("Export Options")
                    col1, col2 = st.columns(2)
                    with open(img_path, "rb") as f:
                        col1.download_button("Download Tree (PNG)", f, "tree.png")
                    newick_data = tree_with_support.format("newick")
                    col2.download_button("Download Newick (.nwk)", newick_data, "tree.nwk")
    else:
        st.info("Upload an alignment to build a tree.")
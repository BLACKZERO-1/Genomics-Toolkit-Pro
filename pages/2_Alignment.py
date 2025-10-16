import streamlit as st
from Bio import SeqIO
from modules import alignment_tools as atools
from modules import visualization as viz
import io

st.set_page_config(page_title="Sequence Alignment", layout="wide")

st.title("⛓️ Sequence Alignment Tools")

# --- Input Section ---
st.header("1. Input Sequences")

# --- NEW CODE ADDED HERE ---
# Initialize session state for the text input
if 'alignment_input' not in st.session_state:
    st.session_state.alignment_input = ""

# Add a section for loading sample data
st.subheader("Load Sample Data")
if st.button("Load Sample Alignment"):
    try:
        with open("data/sample_multi.fasta", "r") as f:
            # Read the file content and store it in session state
            st.session_state.alignment_input = f.read()
        st.success("Sample alignment data loaded!")
    except FileNotFoundError:
        st.error("Error: sample_multi.fasta not found in the 'data' folder.")

st.subheader("Enter Your Data")
# --- END OF NEW CODE ---

# The text_area now uses session_state to hold its value
input_text = st.text_area(
    "Paste your multi-FASTA sequences here",
    value=st.session_state.alignment_input, # This links the text area to the button
    height=250,
    placeholder=">Seq1\nACGT...\n>Seq2\nTGCA..."
)
st.write("--- OR ---")
uploaded_file = st.file_uploader(
    "Upload a multi-FASTA file",
    type=["fasta", "fa"]
)

# Parse sequences and store them in session state for persistence
if 'sequences' not in st.session_state:
    st.session_state.sequences = []
if 'alignment' not in st.session_state:
    st.session_state.alignment = None

fasta_content = None
if input_text:
    fasta_content = input_text
elif uploaded_file is not None:
    fasta_content = uploaded_file.getvalue().decode("utf-8")

if fasta_content:
    try:
        stringio = io.StringIO(fasta_content)
        st.session_state.sequences = list(SeqIO.parse(stringio, "fasta"))
        st.session_state.alignment = None # Reset alignment when new sequences are loaded
        st.success(f"Successfully loaded {len(st.session_state.sequences)} sequences.")
    except Exception as e:
        st.error(f"Error parsing FASTA content: {e}")

st.divider()

# --- Pairwise Alignment Section ---
st.header("2. Pairwise Alignment")
if st.session_state.sequences:
    seq_records = {rec.id: rec for rec in st.session_state.sequences}
    seq_ids = list(seq_records.keys())
    selected_id_1 = st.selectbox("Select the first sequence", seq_ids)
    selected_id_2 = st.selectbox("Select the second sequence", seq_ids, index=min(1, len(seq_ids)-1))
else:
    st.selectbox("Select the first sequence", ["No sequences loaded"])
    st.selectbox("Select the second sequence", ["No sequences loaded"])

st.subheader("Alignment Parameters")
align_type = st.radio("Select Alignment Type", ("Global", "Local"))
col1, col2, col3 = st.columns(3)
match_score = col1.number_input("Match Score", value=1.0)
mismatch_penalty = col2.number_input("Mismatch Penalty", value=-1.0)
gap_penalty = col3.number_input("Gap Penalty", value=-0.5)

if st.button("Perform Pairwise Alignment"):
    if st.session_state.sequences:
        seq1_record = seq_records[selected_id_1]
        seq2_record = seq_records[selected_id_2]
        
        alignment_result = atools.pairwise_align(
            str(seq1_record.seq), str(seq2_record.seq), match_score, mismatch_penalty, gap_penalty, align_type
        )
        st.subheader("Pairwise Alignment Result")
        st.code(alignment_result, language="text")
    else:
        st.error("Please paste or upload sequences first to perform pairwise alignment.")

st.divider()

# --- Multiple Sequence Alignment Section ---
st.header("3. Multiple Sequence Alignment (MSA)")
msa_method = st.radio("Select MSA Tool", ("ClustalW", "MUSCLE"))

if st.button("Perform Multiple Sequence Alignment"):
    if st.session_state.sequences:
        if len(st.session_state.sequences) > 2:
            with st.spinner(f"Running {msa_method}..."):
                st.session_state.alignment = atools.multiple_align(st.session_state.sequences, msa_method)
            if st.session_state.alignment:
                st.success("Multiple Sequence Alignment complete.")
            else:
                st.error("MSA failed. Check the error message above.")
        else:
            st.warning("Please provide at least three sequences for multiple alignment.")
    else:
        st.error("Please paste or upload sequences first to perform MSA.")

st.divider()

# --- MSA Analysis Section ---
st.header("4. Alignment Analysis")
st.write("Results from the last successful Multiple Sequence Alignment will be shown here.")

if st.session_state.alignment:
    alignment = st.session_state.alignment
    
    st.subheader("MSA Result")
    st.text_area("Alignment", str(alignment), height=400)
    st.download_button(
        label="Download Alignment (Clustal format)",
        data=str(alignment),
        file_name="alignment.aln",
        mime="text/plain"
    )
    
    st.subheader("Consensus Sequence")
    consensus_seq = atools.generate_consensus(alignment)
    st.code(consensus_seq)
    
    st.subheader("Conservation Score Plot")
    conservation_scores = atools.calculate_conservation(alignment)
    fig_conservation = viz.plot_alignment_conservation(conservation_scores)
    st.plotly_chart(fig_conservation, use_container_width=True)
    
    st.subheader("Pairwise Identity Matrix")
    seq_ids_alignment = [rec.id for rec in alignment]
    identity_mat = atools.identity_matrix(alignment)
    fig_heatmap = viz.plot_identity_heatmap(identity_mat, seq_ids_alignment)
    st.plotly_chart(fig_heatmap, use_container_width=True)
else:
    st.info("Run a Multiple Sequence Alignment to see the analysis results.")
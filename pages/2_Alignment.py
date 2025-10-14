import streamlit as st
from Bio import SeqIO
from modules import alignment_tools as atools
import io

st.set_page_config(page_title="Sequence Alignment", layout="wide")

st.title("⛓️ Sequence Alignment Tools")

# --- Input Section ---
st.header("1. Input Sequences")

input_text = st.text_area(
    "Paste your multi-FASTA sequences here",
    height=250,
    placeholder=">Seq1\nACGT...\n>Seq2\nTGCA..."
)

st.write("--- OR ---")

uploaded_file = st.file_uploader(
    "Upload a multi-FASTA file",
    type=["fasta", "fa"]
)

sequences = []
fasta_content = None

if input_text:
    fasta_content = input_text
elif uploaded_file is not None:
    fasta_content = uploaded_file.getvalue().decode("utf-8")

if fasta_content:
    try:
        stringio = io.StringIO(fasta_content)
        sequences = list(SeqIO.parse(stringio, "fasta"))
        st.success(f"Successfully loaded {len(sequences)} sequences.")
    except Exception as e:
        st.error(f"Error parsing FASTA content: {e}")

st.divider()

if sequences:
    # --- Pairwise Alignment Section ---
    st.header("2. Pairwise Alignment")

    # ** THE FIX: Create a dictionary of records and a list of IDs for the selectbox **
    seq_records = {rec.id: rec for rec in sequences}
    seq_ids = list(seq_records.keys())

    # Use the simple list of IDs in the selectbox
    selected_id_1 = st.selectbox("Select the first sequence", seq_ids)
    selected_id_2 = st.selectbox("Select the second sequence", seq_ids, index=min(1, len(seq_ids)-1))

    # Retrieve the full SeqRecord object using the selected ID
    seq1_record = seq_records[selected_id_1]
    seq2_record = seq_records[selected_id_2]

    # Alignment parameters
    st.subheader("Alignment Parameters")
    align_type = st.radio("Select Alignment Type", ("Global", "Local"))
    col1, col2, col3 = st.columns(3)
    match_score = col1.number_input("Match Score", value=1.0)
    mismatch_penalty = col2.number_input("Mismatch Penalty", value=-1.0)
    gap_penalty = col3.number_input("Gap Penalty", value=-0.5)

    if st.button("Perform Pairwise Alignment"):
        if seq1_record and seq2_record:
            seq1 = str(seq1_record.seq)
            seq2 = str(seq2_record.seq)
            
            alignment_result = atools.pairwise_align(
                seq1, seq2, match_score, mismatch_penalty, gap_penalty, align_type
            )
            
            st.subheader("Alignment Result")
            st.code(alignment_result, language="text")
        else:
            st.warning("Please select two sequences to align.")
else:
    st.info("Paste sequences or upload a multi-FASTA file to begin.")
import streamlit as st
from Bio import AlignIO
from modules import phylo_tools as ptools
import io

st.set_page_config(page_title="Phylogenetics", layout="wide")

st.title("🌳 Phylogenetic Analysis")

# --- Input Section ---
st.header("1. Input Alignment")

# --- NEW CODE ADDED HERE ---
# Initialize session state to hold the text input
if 'phylo_input' not in st.session_state:
    st.session_state.phylo_input = ""

# Add a section for loading sample data
st.subheader("Load Sample Data")
if st.button("Load Sample Alignment Data"):
    try:
        with open("data/sample_multi.fasta", "r") as f:
            # Read the file content and store it in session state
            st.session_state.phylo_input = f.read()
        st.success("Sample alignment data loaded!")
    except FileNotFoundError:
        st.error("Error: sample_multi.fasta not found in the 'data' folder.")

st.subheader("Enter Your Data")
# --- END OF NEW CODE ---

# The text_area now uses the value from session_state
input_text = st.text_area(
    "Paste your pre-aligned sequences here (Clustal or FASTA format)",
    value=st.session_state.phylo_input, # This links the text area to the button
    height=300,
    placeholder=">Seq1\nACGT-G\n>Seq2\nAC-TTG\n..."
)
st.write("--- OR ---")
uploaded_file = st.file_uploader(
    "Upload a pre-aligned FASTA or Clustal (.aln) file",
    type=["fasta", "fa", "aln"]
)

alignment = None
alignment_content = None
if input_text:
    alignment_content = input_text
elif uploaded_file is not None:
    alignment_content = uploaded_file.getvalue().decode("utf-8")

if alignment_content:
    stringio = io.StringIO(alignment_content)
    try:
        alignment = AlignIO.read(stringio, "clustal")
        st.success("Clustal alignment loaded successfully.")
    except ValueError:
        stringio.seek(0)
        try:
            alignment = AlignIO.read(stringio, "fasta")
            st.success("FASTA alignment loaded successfully.")
        except Exception as e:
            st.error(f"Error parsing alignment file: {e}")

st.divider()

# --- Tree Construction Section ---
if alignment:
    st.header("2. Tree Construction")

    method = st.radio("Select Tree Construction Method", ["Neighbor-Joining (NJ)", "UPGMA"])

    bootstrap_replicates = st.slider("Number of Bootstrap Replicates", min_value=10, max_value=1000, value=100, step=10)

    if st.button("Construct Phylogenetic Tree"):
        with st.spinner("Building tree and running bootstrap... This may be slow."):
            main_tree = None
            tree_method_short = 'nj'
            if method == "Neighbor-Joining (NJ)":
                main_tree = ptools.build_nj_tree(alignment)
                tree_method_short = 'nj'
            elif method == "UPGMA":
                main_tree = ptools.build_upgma_tree(alignment)
                tree_method_short = 'upgma'

            if isinstance(main_tree, str):
                st.error(f"Failed to build tree: {main_tree}")
            else:
                tree_with_support = ptools.apply_bootstrap(main_tree, alignment, bootstrap_replicates, tree_method_short)

                image_path = ptools.draw_tree(tree_with_support)

                st.subheader("Phylogenetic Tree")
                st.success(f"Bootstrap analysis with {bootstrap_replicates} replicates completed. Support values are shown on branches.")
                st.image(image_path)

                st.divider()
                st.subheader("Export Options")
                col1, col2 = st.columns(2)

                # Column 1: Download Tree Image (PNG)
                with open(image_path, "rb") as file:
                    col1.download_button(
                        label="Download Tree as PNG",
                        data=file,
                        file_name="phylogenetic_tree.png",
                        mime="image/png"
                    )

                # Column 2: Download Tree Data (Newick format)
                newick_data = tree_with_support.format("newick")
                col2.download_button(
                    label="Download Tree as Newick (.nwk)",
                    data=newick_data,
                    file_name="phylogenetic_tree.nwk",
                    mime="text/plain"
                )
else:
    st.info("Paste or upload an alignment file to begin.")
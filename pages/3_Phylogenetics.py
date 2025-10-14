import streamlit as st
from Bio import AlignIO
from modules import phylo_tools as ptools
import io

st.set_page_config(page_title="Phylogenetics", layout="wide")

st.title("🌳 Phylogenetic Analysis")

# --- Input Section ---
st.header("1. Input Alignment")
input_text = st.text_area(
    "Paste your pre-aligned sequences here (Clastal or FASTA format)",
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

    if st.button("Construct Phylogenetic Tree"):
        with st.spinner("Building and visualizing tree..."):
            tree = ptools.build_nj_tree(alignment)

            if isinstance(tree, str):
                st.error(f"Failed to build tree: {tree}")
            else:
                # The new draw function is simpler and doesn't need a style
                image_path = ptools.draw_tree(tree)

                st.subheader("Phylogenetic Tree")
                st.image(image_path)

                with open(image_path, "rb") as file:
                    st.download_button(
                        label="Download Tree Image (PNG)",
                        data=file,
                        file_name="phylogenetic_tree.png",
                        mime="image/png"
                    )
else:
    st.info("Paste or upload an alignment file to begin.")
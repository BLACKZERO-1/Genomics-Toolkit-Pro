import streamlit as st

# Set the page configuration
st.set_page_config(
    page_title="Genomics Toolkit Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Main Homepage Content ---
st.title("🧬 Genomics Toolkit Pro")
st.markdown("Welcome to your one-stop solution for advanced genomics analysis. This toolkit provides a suite of powerful tools for sequence analysis, alignment, and phylogenetics.")

st.divider()

# --- Navigation Cards to other pages ---
st.header("Toolkit Modules")
st.write("Navigate to the different modules using the sidebar or the cards below.")

col1, col2, col3 = st.columns(3, gap="large")

with col1:
    with st.container(border=True):
        st.subheader("1. Sequence Analysis")
        st.write("Perform detailed analysis of DNA, RNA, or protein sequences. Calculate GC content, find ORFs, design primers, and much more.")

with col2:
    with st.container(border=True):
        st.subheader("2. Sequence Alignment")
        st.write("Align multiple sequences using industry-standard algorithms like ClustalW and MUSCLE to identify conserved regions and evolutionary relationships.")

with col3:
    with st.container(border=True):
        st.subheader("3. Phylogenetics")
        st.write("Construct and visualize phylogenetic trees to explore the evolutionary history of your sequences using methods like Neighbor-Joining and UPGMA.")

st.divider()

# --- Quick Stats Dashboard (Placeholder) ---
st.header("Project Dashboard")
col1, col2, col3 = st.columns(3)
col1.metric("Tools Implemented", "20+", "Sequence, Alignment, Phylo")
col2.metric("Analyses Performed", "1,245", "7%")
col3.metric("Sample Datasets", "3", "FASTA format")
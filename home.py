import streamlit as st

# --- Load Custom CSS ---
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
local_css("style.css")

# --- Sidebar Message ---
st.sidebar.success("Select a tool from the pages above.")

# --- Main Homepage Content ---
st.title("🧬 Genomics Toolkit Pro")
st.markdown("Welcome to your one-stop solution for advanced genomics analysis.")
st.divider()
st.header("Toolkit Modules")
# ... (Your homepage content from before) ...
col1, col2, col3 = st.columns(3, gap="large")
with col1:
    with st.container(border=True):
        st.subheader("🔬 Sequence Analysis")
        st.write("Perform detailed analysis of DNA, RNA, or protein sequences.")
with col2:
    with st.container(border=True):
        st.subheader("⛓️ Sequence Alignment")
        st.write("Align multiple sequences using algorithms like ClustalW and MUSCLE.")
with col3:
    with st.container(border=True):
        st.subheader("🌳 Phylogenetic Analysis")
        st.write("Construct and visualize phylogenetic trees from alignments.")
# 🧬 Genomics Toolkit Pro

A web-based toolkit for performing common genomics and bioinformatics tasks, built with Python and Streamlit. This application allows users to analyze, align, and build phylogenetic trees from biological sequences.

## Features

The toolkit is organized into three main modules:

1.  **🔬 Sequence Analysis:**
    * Load sequence data by pasting text or uploading a FASTA file.
    * Calculate basic stats like sequence length and GC content.
    * Transcribe DNA to RNA and translate DNA in all six reading frames.
    * Find Open Reading Frames (ORFs) with a user-defined minimum length.
    * Analyze GC content distribution with a sliding window plot.
    * Visualize codon usage frequency with an interactive heatmap.
    * Identify restriction enzyme cut sites.
    * Design basic forward and reverse primers for a target region.
    * Search for sequence motifs using regular expressions.
    * Analyze k-mer frequency distribution.

2.  **⛓️ Sequence Alignment:**
    * Load multiple sequences by pasting a multi-FASTA or uploading a file.
    * Perform pairwise alignments (Global and Local) between any two sequences.
    * Run multiple sequence alignments using integrated command-line tools (ClustalW and MUSCLE).
    * Analyze alignment results with a consensus sequence, a conservation score plot, and a pairwise identity matrix heatmap.

3.  **🌳 Phylogenetic Analysis:**
    * Load a pre-aligned set of sequences.
    * Construct phylogenetic trees using Neighbor-Joining (NJ) and UPGMA methods.
    * Perform bootstrap analysis to assess the statistical support for tree branches.
    * Visualize the final tree with bootstrap values displayed on the branches.
    * Export the tree as a PNG image or in the standard Newick (.nwk) format.

## Setup and Installation (WSL - Ubuntu)

Follow these steps to run the application locally.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/BLACKZERO-1/Genomics-Toolkit-Pro.git](https://github.com/BLACKZERO-1/Genomics-Toolkit-Pro.git)
    cd Genomics-Toolkit-Pro
    ```

2.  [cite_start]**Create and activate the Conda environment:** [cite: 33-34]
    ```bash
    conda create -n genomics python=3.10 -y
    conda activate genomics
    ```

3.  [cite_start]**Install dependencies:** [cite: 36-37]
    ```bash
    conda install -c conda-forge -c bioconda biopython clustalw muscle -y
    pip install streamlit pandas numpy plotly scikit-bio ete3 matplotlib PyQt5
    ```

4.  [cite_start]**Run the Streamlit app:** [cite: 40]
    ```bash
    streamlit run app.py
    ```
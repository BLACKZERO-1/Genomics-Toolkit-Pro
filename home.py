import streamlit as st
import base64
from io import BytesIO

# Base64 encoded microscope icon (64x64 px PNG)
microscope_base64 = """
iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAABQ0lEQVR4Xu3aS07DMBBF0U6CsTXiBhZZvkMwTezf8aQIIbMZUAGG1QoQHmAJnnQloEmWCN2s9f98ofZFcD6d3altLUuLNvR+A4vfAOzw7Zv7rqnAwjircoTiRWeWITRWvREmswSYMwkfPQ+GwWg8HM+vNYW8Ea7sp14Pf2hAjIQTzQTCfwk6lbRY6HuHmUUZWMD8diUzxlF2/BwL3o8ptJlsGh4pJlJVGGLY6Ngul++dRha6rkZQDhsKwTxDwKAyIOQ0Qu+fUvpDsSuJU7dZhIxkVrClDwNAYyApAgNG+T7GONlOAKcT6sGFhUlHt3u9a0nlzTziCJQxAQYHWaKt/Bh7ZFHZb213gBvIv87dePUFa5c6r8M/9BjoCfjYTDnqLf+8AAAAASUVORK5CYII=
"""

# Base64 encoded chain link icon (64x64 px PNG)
chain_link_base64 = """
iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAABl0lEQVR4Xu3asU4DMRCA4bcgO+wlPcEkg7KxuKOex7u27fGicyEzSmeQpm0CkVfqNwdEYnk8lkiSJEmSZ7DWfS+kAH18BnITxRfAjgJPZLUXOC54zBTI3gsAeU4PATPBbgCqO3jaAdqeiIt7RyEbMWYFgUQvRa3+LB3B0kxbBf1DcAVUI9GPC8bof4Hp4cVesn9N3v7gcwJ31LdoKYV7SmFA9gfHR7zDyhHcwcTNQw5GfZHNJ0pg/kAi0Tvn0VVwIjwrOT3ZTwosHJ8hWxLn37SEA9y+QhD8cP+tbNSR/7spaSO3hY4XxiYznsUUQe+w5LhjsUM07MdFfYB+iNCXeHsGIY+ZxjuMgByl8ZPeXh8Z84ibNyCbX3c3TpL5a4KM+3rQZMqzA5ufodb5NUMrIhoTnxfQvDbEaGXPNL5JiaDwRmJ2sacfTb6Evzn+lm9EmSZIkSZLko/8Am0uVk+7PnY9wAAAABJRU5ErkJggg==
"""

# Base64 encoded tree icon (64x64 px PNG)
tree_base64 = """
iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAABNklEQVR4Xu3aMQ6CQBQE0Mv0pNnCbOBEoGzBRxTbhrMPbDJLxAxPztAXZW9OhRpjFUjPLOJkjkfi+lf2zUdOSvvA7gAPhHC8BU+eA0wCF4HfYO7AWOPYMUnp8ZjdGAyxoXCWyKbgaTqSNocUcfMOKjdWMU8MDmkf5wAe7SPDjNXCHtZPlhpYsL8I5UqJIboqybGi2COVEo5nRhNOlnkYcgxdhEwlyMi5H9jHC8Q6e0O8HIZGeKro/kPUaUnvni6LSjMlLvmJmFtRFf6Dj3SrbF+jW/dPY9ql5jq6syMzM37R1N2ZchI57gPbyd6HZ7P7KxZvtf+BnuCxbLv53gAAAMAfz4LXqdjDl9qGp4UAAAAASUVORK5CYII=
"""

# Using the icons again for section images as placeholders — swap for better images if you have
seq_analysis_img_base64 = microscope_base64
seq_alignment_img_base64 = chain_link_base64
phylo_analysis_img_base64 = tree_base64

def base64_to_bytesio(base64_str):
    return BytesIO(base64.b64decode(base64_str))

def home_page():
    st.markdown(
        """
        <style>
        .main-title {
            font-family: 'Helvetica Neue', sans-serif;
            font-size: 48px;
            background: linear-gradient(90deg, #00c6ff, #0072ff);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            font-weight: bold;
            margin-bottom: 0;
        }
        .subtitle {
            font-family: 'Helvetica Neue', sans-serif;
            font-size: 18px;
            color: #666666;
            margin-top: 0;
            margin-bottom: 40px;
        }
        .feature-title {
            font-weight: 600;
            font-size: 22px;
            color: #0072ff;
            margin-bottom: 10px;
        }
        .feature-description {
            font-size: 16px;
            color: #333333;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown('<h1 class="main-title">🧬 Genomics Toolkit Pro</h1>', unsafe_allow_html=True)
    st.markdown(
        '<p class="subtitle">A web-based toolkit for performing common genomics and bioinformatics tasks, built with Python and Streamlit.</p>',
        unsafe_allow_html=True,
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.image(base64_to_bytesio(microscope_base64), width=60)
        st.markdown('<div class="feature-title">Sequence Analysis</div>', unsafe_allow_html=True)
        st.markdown(
            """
            <div class="feature-description">
            Load, analyze and visualize sequence data including ORFs, GC content, codon usage, and primer design.
            </div>""",
            unsafe_allow_html=True,
        )
        st.image(base64_to_bytesio(seq_analysis_img_base64), use_column_width=True)

    with col2:
        st.image(base64_to_bytesio(chain_link_base64), width=60)
        st.markdown('<div class="feature-title">Sequence Alignment</div>', unsafe_allow_html=True)
        st.markdown(
            """
            <div class="feature-description">
            Perform global/local pairwise and multiple sequence alignments with tools like ClustalW and MUSCLE.
            </div>""",
            unsafe_allow_html=True,
        )
        st.image(base64_to_bytesio(seq_alignment_img_base64), use_column_width=True)

    with col3:
        st.image(base64_to_bytesio(tree_base64), width=60)
        st.markdown('<div class="feature-title">Phylogenetic Analysis</div>', unsafe_allow_html=True)
        st.markdown(
            """
            <div class="feature-description">
            Build and visualize phylogenetic trees with bootstrap analysis and export options.
            </div>""",
            unsafe_allow_html=True,
        )
        st.image(base64_to_bytesio(phylo_analysis_img_base64), use_column_width=True)

    # Optional: Add 3D model or animation video here if you have one
    # st.video('assets/3d_model_animation.mp4')
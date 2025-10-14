from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt

def build_nj_tree(alignment):
    """
    Builds a Neighbor-Joining (NJ) phylogenetic tree from a multiple alignment.
    """
    try:
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dist_matrix)
        tree.root_at_midpoint()
        return tree
    except Exception as e:
        return str(e)

def draw_tree(tree):
    """
    Renders a phylogenetic tree using Biopython's draw function with matplotlib.
    Returns the path to the image file.
    """
    image_path = "tree.png"

    # Create a matplotlib figure
    fig = plt.figure(figsize=(10, 15), dpi=100)
    ax = fig.add_subplot(1, 1, 1)

    # Draw the tree on the axes
    Phylo.draw(tree, axes=ax, do_show=False)

    # Save the figure to a file
    plt.savefig(image_path)
    plt.close(fig) # Close the figure to free up memory

    return image_path
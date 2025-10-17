from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo
from Bio.Phylo.Consensus import get_support
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import random

def build_nj_tree(alignment):
    """Builds a Neighbor-Joining (NJ) phylogenetic tree."""
    try:
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(alignment)
        tree.root_at_midpoint()
        return tree
    except Exception as e:
        return str(e)

def build_upgma_tree(alignment):
    """Builds a UPGMA phylogenetic tree."""
    try:
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'upgma')
        tree = constructor.build_tree(alignment)
        tree.root_at_midpoint()
        return tree
    except Exception as e:
        return str(e)

def apply_bootstrap(main_tree, alignment, replicates, method):
    """
    Manually performs bootstrapping and adds support values to the main tree.
    """
    bootstrap_alignments = []
    aln_len = alignment.get_alignment_length()

    for i in range(replicates):
        indices = [random.randrange(aln_len) for _ in range(aln_len)]
        new_records = []
        for record in alignment:
            new_seq_str = "".join([record.seq[i] for i in indices])
            new_records.append(SeqRecord(Seq(new_seq_str), id=record.id, description=record.description))
        bootstrap_alignments.append(MultipleSeqAlignment(new_records))

    trees = []
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, method)
    for aln in bootstrap_alignments:
        tree = constructor.build_tree(aln)
        trees.append(tree)

    support_tree = get_support(main_tree, trees)

    return support_tree

def draw_tree(tree):
    """Renders a phylogenetic tree with bootstrap values using matplotlib."""
    image_path = "tree.png"
    fig = plt.figure(figsize=(10, 15), dpi=100)
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)

    plt.savefig(image_path)
    plt.close(fig)
    return image_path
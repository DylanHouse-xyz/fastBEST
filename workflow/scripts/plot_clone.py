import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import pandas as pd


def get_clone(cluster_file):
    """
    Reads input csv file and returns a pandas series 
    """
    df = pd.read_csv(cluster_file, header = None)
    return df.set_index(0)[1]


def build_and_draw_clone_tree(adjacency_list_file, mutation_to_clone, ax=None):
    """
    Parses mutation relationships, collapses them into clones, and plots the tree.

        Args:
            adjacency_list_file: File that contains the adjacency list output from fastbe search
            mutation_to_clone: CSV file containing cluster to mutation data
        Returns:
            A clonal tree
                
    """
    clone_tree = nx.DiGraph()

    with open(adjacency_list_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            mutations = line.split()
            # The first mutation is the ancestral/parent variant
            ancestor_clone = mutation_to_clone.get(mutations[0])

            # Subsequent mutations on the line are descendent variants
            for descendent_mut in mutations[1:]:
                descendent_clone = mutation_to_clone.get(descendent_mut)

                # Connect only if the mutation meets a distinct sub-clone
                if ancestor_clone != descendent_clone:
                    clone_tree.add_edge(ancestor_clone, descendent_clone)
                    parents = list(clone_tree.predecessors(descendent_clone))
                    print(f"Parents of {descendent_clone}: {parents}")
                    
    

                    
    
    
    # Draw the phylogenetic tree
    pos = graphviz_layout(clone_tree, prog="dot")
    labels = {node: node for node in clone_tree.nodes()}
    nx.draw(clone_tree, pos, with_labels=False, arrows=True, ax=ax, node_color="lightblue", node_size=800)
    nx.draw_networkx_labels(clone_tree, pos, labels, font_size=10, ax=ax)


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct and plot clonal evolutionary trees from mutation adjacency lists."
    )
    parser.add_argument("--cluster",required=True, help="CSV file mapping individual mutations to their designated clone/cluster assignment.")
    parser.add_argument("adjacency_lists",help="One or more text files containing mutation adjacency hierarchies.",nargs="+")
    args = parser.parse_args()

    # Map unique clones to mutations
    mutation_to_clone = get_clone(args.cluster)

    fig, axes = plt.subplots(figsize=(8,6))
    if len(args.adjacency_lists) == 1:
        axes = [axes]

    for i, adjacency_list_file in enumerate(args.adjacency_lists):
        build_and_draw_clone_tree(adjacency_list_file, mutation_to_clone, axes[i])

    output_png = "clonal_lineage_trees.png"
    plt.savefig(output_png, dpi=300)


if __name__ == "__main__":
    main()
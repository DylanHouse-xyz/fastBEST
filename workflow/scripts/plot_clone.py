import argparse
import os
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from networkx.drawing.nx_agraph import graphviz_layout


def get_clone(cluster_file):
    """
    Reads input csv file and returns a pandas series
    """
    df = pd.read_csv(cluster_file, header = None)
    return df.set_index(0)[1]


def build_and_draw_clone_tree(adjacency_list_file, mutation_to_clone, output_dir, ax=None):
    """
    Parses mutation relationships, collapses them into clones, and plots the tree. Outputs a clonal tree .png, and a list containing each clones parents in a .txt format.

        Args:
            adjacency_list_file: File that contains the adjacency list output fron fastbe search

            mutation_to_clone: CSV file containing cluster to mutation data

            output_dir = The directory where both outputs will be saved.
        Returns:
            A clonal tree & a .txt file containing a list of parents (0-indexing)

    """
    clone_tree = nx.DiGraph()
    parents_file = []

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

                    if parents:
                        parents_file.append(int(parents[0]))
    



    # Draw the phylogenetic tree
    pos = graphviz_layout(clone_tree, prog="dot")
    labels = {node: node for node in clone_tree.nodes()}
    nx.draw(clone_tree, pos, with_labels=False, arrows=True, ax=ax, node_color="lightblue", node_size=800)
    nx.draw_networkx_labels(clone_tree, pos, labels, font_size=10, ax=ax)

    # Write list of parents to parents.txt.
    output_file_path = os.path.join(output_dir, 'parents.txt')
    with open(output_file_path, 'w') as p:
        p.write(str(parents_file))

def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct and plot clonal evolutionary trees from mutation adjacency lists."
    )
    parser.add_argument("--cluster",required=True, help="CSV file mapping individual mutations to their designated clone/cluster assignment.")
    parser.add_argument("adjacency_list",help="The adjacency list output from fastBE search")
    parser.add_argument('-o', '--output',help="The directory in which the clonal tree will be stored." )
    args = parser.parse_args()

    # Map unique clones to mutations
    mutation_cluster = get_clone(args.cluster)
    build_and_draw_clone_tree(args.adjacency_list, mutation_cluster, args.output)
    output_png = os.path.join(args.output, 'clone_tree.png')
    plt.savefig(output_png, dpi=300)


if __name__ == "__main__":
    main()

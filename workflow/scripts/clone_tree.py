import argparse
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd  
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph(adjacency_list_file, cluster_file, ax=None):
    """
Reads the adjacency list and the clustering results to draw a directed mutation tree. The nodes are colored according to their assigned clone, and the edges represent the parent-child relationships between mutations.

        Args:
            adjacency_list_file (str): The path to the .txt file containing the adjacency list of the mutation tree.
            mutation_map_file (str): The path to the .csv file containing the mutation and clone mapping.
    """
    T = nx.DiGraph()

    # Reads the cluster .csv file that contains both mutations & assigned clones
    df = pd.read_csv(cluster_file)
    mutation_to_clone = dict(zip(df['mutation'], df['clone']))

    with open(adjacency_list_file, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            
            nodes = list(map(int, line.split()))
            parent_node = nodes[0]
            for node in nodes[1:]:
                T.add_edge(parent_node, node)

    # Assigning a unique color to a clone.
    unique_clones = sorted(df['clone'].unique())
    color_map_plt = plt.get_cmap('tab20')
    clone_colors = {clone: color_map_plt(i / max(1, len(unique_clones)-1)) 
                    for i, clone in enumerate(unique_clones)}

    node_colors = [clone_colors.get(mutation_to_clone.get(node), 'lightgrey') 
                   for node in T.nodes()]

    pos = graphviz_layout(T, prog='dot')
    
    nx.draw(
        T, pos, 
        with_labels=False, 
        arrows=True,
        arrowsize=1,
        alpha=0.7, 
        ax=ax, 
        node_color=node_colors,  
        node_size=5,
        edge_color='gray'
    )

def main():
    parser = argparse.ArgumentParser(description="Draw a colored clone tree")
    parser.add_argument("adjacency_list", help="The .txt file containing the adjacency list")
    parser.add_argument("cluster_file", help="The .csv file containing mutation and clone columns")
    args = parser.parse_args()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title(f"Tree: {args.adjacency_list}")
    
    draw_graph(args.adjacency_list, args.cluster_file, ax)

    plt.tight_layout()

    output_filename = args.adjacency_list.replace(".txt", "_tree.png")
    plt.savefig(output_filename, dpi=300)
    print(f"Tree visualization saved to {output_filename}")

if __name__ == "__main__":
    main()
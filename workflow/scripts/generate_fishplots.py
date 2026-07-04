from collections import defaultdict, deque

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fishplotpy import FishPlotData, fishplot
import copy

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Parses fastBEST intermediate files and produces tumour series fish plots and final data results.
        """
    )
    parser.add_argument(
        "--ccf-file",
        help="Raw CCF file calculated from frequency matrix and clusters.",
        required=True
    )
    parser.add_argument(
        "--parent-file",
        help="Two-column file connecting child clones to their parents.", #A required parameter. Script will refuse to execute without it.
        required=True
    )
    parser.add_argument(
        "--clusters-file",
        help="(Optional) File from fastBE associating mutations with clusters. This is useful to ensure that indices match, as clones will undergo reindexing during the plotting process.",
        required=False,
        default=None
    )
    parser.add_argument(
        "--output",
        help="(Optional) Path to output for fish plot.",
        default="./fish.png"
    )
    parser.add_argument(
        "--sample-order",
        help="(Optional) Comma-separated list of index values for the final output order. (i.e., 3,2,1,5,6) Omitted options will not be plotted.",
    )

    return parser.parse_args()

def recursive_scale_descendants(clone_dict:dict, clone_key:str, sample_index:int, scale:float):
    clone_dict[clone_key]['sample_freq'][sample_index] *= scale
    #print(f'Scaled {clone_key} by {scale} ({clone_dict[clone_key]["sample_freq"][sample_index]})')
    for d in clone_dict[clone_key]['descendants']:
        recursive_scale_descendants(clone_dict, d, sample_index, scale)

def adjust_clone_frequencies(clone_dict:dict):
    n_samples = len(clone_dict[list(clone_dict.keys())[0]]['sample_freq'])
    for sample in range(n_samples):
        # Increase parents that are too low.
        for i in reversed(sorted(list(clone_dict.keys()))):
            descendants = clone_dict[i]['descendants']
            own_freq = clone_dict[i]['sample_freq'][sample]
            sum_children = 0
            for d in descendants:
                sum_children += clone_dict[d]['sample_freq'][sample]
            #print(f'{i}->{descendants} : {sum_children}')
            if sum_children > 0 and sum_children >= own_freq:
                print(f'Sum of children {",".join([str(s) for s in descendants])} ({sum_children:.3f}) exceeds parent {i} ({own_freq:.3f}) in sample {sample}. Adjusting...')
                clone_dict[i]['sample_freq'][sample] = sum_children * 1.01
                #print(clone_dict[i]['sample_freq'][sample])
        
        # Decrease all that are too high.
        for i in sorted(list(clone_dict.keys())):
            if clone_dict[i]['descent_level'] == 0:
                descendants = clone_dict[i]['descendants']
                own_freq = clone_dict[i]['sample_freq'][sample]
                if own_freq > 1.0:
                    print(f'Frequency of clone {i} in sample {sample} ({own_freq:.2f}) > 1.0. Reducing...')
                    scale_factor = 1/own_freq
                    recursive_scale_descendants(clone_dict, i, sample, scale_factor)

def plot_clone_network(children:list, parents:list):
    ### Plots a network of clones. Used for debugging; currently not very good. ###
    import networkx as nx
    import matplotlib.pyplot as plt
    
    # 1. Your existing lists
    nodes = list(children)
    edges = []

    for i in range(len(children)):
        if parents[i] != -1:
            edges.append( (children[i],parents[i]) )
    
    # 2. Initialize an empty Graph object
    G = nx.Graph()
    
    # 3. Add your nodes and edges
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    # 4. Draw the network
    plt.figure(figsize=(8, 6))
    
    # networkx provides various layouts (spring, circular, shell, etc.)
    pos = nx.spring_layout(G, seed=42) 
    
    nx.draw(
        G, pos, 
        with_labels=True, 
        node_color='lightblue', 
        node_size=1200, 
        font_size=10, 
        font_weight='bold', 
        edge_color='gray'
    )


def generate_fishplots(ccf_file:str, parents_file:str, clusters_file:str=None, output_file:str=None, sample_order:list=None):
    
    if output_file is None:
        output_file = "./fish.png"

    # Initial load-in of data
    clone_dict = {}

    root_clone = 0
    tree_depth = 0

    founder_clones = []

    lineage_data = pd.read_csv(parents_file).sort_values("Clone_Index")
    #print(lineage_data)

    #parents,children,reindex_dict = reindex_tumour_tree(lineage_data["Parent_Clone"],lineage_data["Clone_Index"])
    parents = lineage_data["Parent_Clone"]
    children = lineage_data["Clone_Index"]

    #print(parents)
    #print(children)

    clust_dict = {}

    # If the optional clusters file is provided, associate the mutations with the clones before they're reindexed.ß
    if clusters_file is None:
        clust_dict = {}
    else:
        var_clusters = pd.read_csv(clusters_file,sep=",")
        for i in range(len(var_clusters["clone"])):
            row = var_clusters.iloc[i]
            clust_dict[int(row["mutation"])] = int(row["clone"])

    ## Re-enable this when we actually make the network useful.
    #plot_clone_network(children,parents)

    ancestor_dict = {}

    #parents = [0, 1, 1, 3]
    parents_clean = []
    children_clean = []

    for i in range(len(parents)):
        clone = int(children[i])
        parent = int(parents[i])
        ancestor_dict[clone] = parent

        #print(f"{parent}->{clone}")
        
        if parent == -1:
            root_clone = clone
        else:
            if parent == root_clone:
                founder_clones.append(clone)
            parents_clean.append(parent)
            children_clean.append(clone)

    for i in children:
        clone_dict[i] = {'parent':ancestor_dict[i],'full_ancestry':[],'descendants':[]}
        next_ancestor = clone_dict[i]['parent']
        while next_ancestor != -1:
            clone_dict[i]['full_ancestry'].append(next_ancestor)
            next_ancestor = ancestor_dict[next_ancestor]
        clone_dict[i]['descent_level'] = len(clone_dict[i]['full_ancestry'])
        tree_depth = max(clone_dict[i]['descent_level'],tree_depth)

    for i in clone_dict:
        parent = clone_dict[i]['parent']
        if parent in clone_dict:
            clone_dict[parent]['descendants'].append(i)

    frac_data = pd.read_csv(ccf_file) # .drop(str(root_clone),axis=1)

    n_samples = len(frac_data["Sample_Index"])

    # If an explicit list wasn't given, go with what was in the input file.
    if sample_order is None:
        sample_order = list(frac_data["Sample_Index"])
    else:
        sample_order = [int(x) for x in sample_order.split(",")]

    n_clones = len(children_clean)

    for i in clone_dict:
        clone_dict[i]["sample_freq"]={}
        for j in range(n_samples):
            clone_dict[i]["sample_freq"][j] = float(frac_data[str(i)][j])
        clone_dict[i]['variants'] = []
        for k in clust_dict:
            if clust_dict[k] == i:
                clone_dict[i]['variants'].append(k)

    adj_clone_dict = copy.copy(clone_dict)

    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    convert_id_dict = {-1:"@"}

    count = 0
    for i in letters:
        convert_id_dict[count]=i
        count += 1

    #print(convert_id_dict)

    letter_clone_dict = {}
    letter_level_dict = {}

    for i in adj_clone_dict:
        letter_id = convert_id_dict[i]
        clone_level = adj_clone_dict[i]['descent_level']
        letter_level_dict[letter_id]=clone_level
        parent = adj_clone_dict[i]['parent']
        letter_clone_dict[letter_id] = copy.copy(adj_clone_dict[i])
        letter_clone_dict[letter_id]['parent'] = convert_id_dict[parent]
        ancestry = adj_clone_dict[i]['full_ancestry']
        letter_ancestry = [convert_id_dict[x] for x in ancestry]
        letter_clone_dict[letter_id]['full_ancestry'] = letter_ancestry
        descendants = adj_clone_dict[i]['descendants']
        letter_descendants = [convert_id_dict[x] for x in descendants]
        letter_clone_dict[letter_id]['descendants'] = letter_descendants

    #for i in letter_clone_dict:
    #    print(i)
    #    print(letter_clone_dict[i])
    #    print("\n")

    renumber_dict = {}

    root_letter = '$'
    for i in letter_level_dict:
        if letter_level_dict[i] == 0:
            root_letter = i
            #null = letter_clone_dict.pop(i)
    
    renumber_dict[root_letter]=0
    renumber_dict['@']=-1

    n_clones = len(letter_clone_dict)

    current_clone = 1
    current_level = 1

    while current_clone < n_clones:
        for i in letter_level_dict:
            if letter_level_dict[i] == current_level:
                renumber_dict[i] = current_clone
                #print(f"{i} -> {renumber_dict[i]}")
                current_clone += 1
        current_level += 1

    #print(renumber_dict)

    final_clone_dict = {}

    #print(letter_clone_dict)

    for i in letter_clone_dict:
        
        final_id = renumber_dict[i]
        if final_id == -1:
            continue
        parent = letter_clone_dict[i]['parent']
        final_clone_dict[final_id] = copy.copy(letter_clone_dict[i])
        final_clone_dict[final_id]['parent'] = renumber_dict[parent]
        ancestry = letter_clone_dict[i]['full_ancestry']
        number_ancestry = [renumber_dict[x] for x in ancestry]
        final_clone_dict[final_id]['full_ancestry'] = number_ancestry
        descendants = letter_clone_dict[i]['descendants']
        number_descendants = [renumber_dict[x] for x in descendants]
        final_clone_dict[final_id]['descendants'] = number_descendants

    adjust_clone_frequencies(final_clone_dict)

    null = final_clone_dict.pop(renumber_dict[root_letter])

    # Dump adj_clone_dict back to the necessary fishplot data.

    freq_list_all = []
    parents_final = []

    for i in sorted(list(final_clone_dict.keys())):
        #print(i)
        sample_freq_list = []
        for j in sample_order:
            sample_freq_list.append(final_clone_dict[i]['sample_freq'][j])
        freq_list_all.append(sample_freq_list)
        parents_final.append(final_clone_dict[i]['parent'])

    timepoints = list(range(0,len(sample_order)))

    # Convert to the formats that FishplotPy expects and generate plot.

    print(final_clone_dict)

    frac_table = np.array(freq_list_all) * 100 # Expects percentages.

    print("\nFishplots inputs:")
    print(frac_table)
    print(parents_final)
    print(timepoints)

    fp_data = FishPlotData(
    frac_table=frac_table,
    parents=parents_final,
    timepoints=timepoints,
        fix_missing_clones=True
    )

    fp_data.layout_clones()

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Generate the plot onto the axes
    fishplot(fp_data, ax=ax, shape="spline") # Using recommended spline shape

    # 5. Display the plot
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Wrote out file {output_file}")
    #plt.show()

def main():
    """Main execution logic."""
    # 1. Grab the parsed arguments
    args = parse_args()
    #print(args)

    generate_fishplots(ccf_file=args.ccf_file, parents_file=args.parent_file, output_file=args.output, clusters_file=args.clusters_file, sample_order=args.sample_order)


if __name__ == "__main__":
    main()

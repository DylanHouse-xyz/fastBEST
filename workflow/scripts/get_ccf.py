import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_ccf(mutation_to_clone, vaf_matrix, filepath):
    """
    Calculates the cancer cell fraction (CCF) & the clonal composition in
    cancer samples from outputs of fastBE cluster.

        Args:
            mutation_to_clone: CSV file containing cluster to mutation data
            vaf_matrix: variant allele frequency matrix in a .txt format according to fastBE's specifications
            filepath: The path to output.
        Returns:
            A CSV file containing the CCF data per clone and sample.
    """
    threshold = 0.05

    clusters_df = pd.read_csv(mutation_to_clone)
    clusters_df = clusters_df.sort_values("mutation", ascending = True).reset_index(drop = True)
    unique_clusters = sorted(clusters_df["clone"].unique())

    cluster_to_idx = {
        cluster: idx for idx, cluster in enumerate(unique_clusters)
    }

    root_mutation = clusters_df.iloc[0,1]
    clusters_df = clusters_df[clusters_df.iloc[:, 1] != root_mutation].reset_index(drop=True)


    variants_per_cluster = clusters_df["clone"].value_counts()
    cluster_counts_vector = np.array(
        [variants_per_cluster[clone] for clone in unique_clusters]
    )

    vaf_matrix = pd.read_csv(vaf_matrix, sep = " ", header = None).values
    num_samples, num_variants = vaf_matrix.shape


    num_clusters = len(unique_clusters)
    sum_matrix = np.zeros((num_samples, num_clusters))

    # Iterates over matrix of zeros though variant and sample index and
    # replaces value with that of the VAF matrix value. We then find the
    # clone it belongs to.
    for sample_idx in range(num_samples):
        for variant_idx in range(num_variants):

            vaf_val = vaf_matrix[sample_idx, variant_idx]

            # Find which cluster this variant position belongs to
            cluster_name = clusters_df.loc[variant_idx, "clone"]
            cluster_idx = cluster_to_idx[cluster_name]


            sum_matrix[sample_idx, cluster_idx] += vaf_val

    avg_vaf_matrix = sum_matrix / cluster_counts_vector

    # Multiply by 2 to get CCF (assume diploid)
    ccf_matrix = avg_vaf_matrix * 2
    ccf_matrix = np.where(ccf_matrix < threshold, 0.0, ccf_matrix)


    ccf_df = pd.DataFrame(ccf_matrix, columns=unique_clusters)
    ccf_df.index.name = "Sample_Index"
    ccf_df.to_csv(filepath)
    print("CCF file produced in specified directory")
    return ccf_df

def plot_ccf(df):
    df.plot(kind = 'bar', stacked = 'True', title = 'Cancer Cell Fraction Proportion', use_index = True)
    output_png = "ccf.png"
    plt.savefig(output_png, dpi=300)
    
def main():
    parser = argparse.ArgumentParser("Calculates the cancer cell fraction from fastBE's cluster output and VAF matrix")
    parser.add_argument('Cluster_Variant', help = "The cluster-variant .csv file")
    parser.add_argument('VAF_Matrix', help = "The variant allele frequency .txt file")
    parser.add_argument('output', help = "Filepath to output")
    args = parser.parse_args()

    if args.Cluster_Variant and args.VAF_Matrix:
        df = get_ccf(args.Cluster_Variant, args.VAF_Matrix, args.output)
        plot_ccf(df)

    else:
        raise ValueError("Ensure to input both cluster-variant csv file from fastBE output & the VAF .txt matrix in fastBE format")

if __name__ == "__main__":
    main()

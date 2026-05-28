# Get CCF, plot proportion, plot phylogenetic tree, and plot fishplot will
# be here.

rule get_ccf:
    input:
        matrix = rules.vcf_converter.output.matrix
        clone_map = rules.optimized_fastbe_cluster.output.cluster
    output:
        ccf = "results/{tumors}/{tumors}_ccf.csv"
    conda:
        "envs/scripts.yaml"
    shell:
        "python3 scripts/get_ccf.py {input.clone_map} {input.matrix} {output.ccf}"


# Plot phylogenetic tree :)

rule plot_tree:
    input:
        clone_map = rules.optimized_fastbe_cluster.output.cluster
        adjacency_list = rules.fastbe_search.output.tree
    output:
        phylogenetic_tree = "clonal_lineage_trees.png"
    conda:
        "envs/scripts.yaml"
    shell:
        "mkdir -p results{tumors}/tree && "
        "python3 scripts/plot_clone.py --cluster {input.clone_map} {input.adjacency_list}"

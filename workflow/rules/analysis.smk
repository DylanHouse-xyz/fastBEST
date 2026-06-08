# Get CCF, plot proportion, plot phylogenetic tree, and plot fishplot will
# be here.

rule get_ccf:
    input:
        matrix = rules.append_column.output.af_matrix,
        clone_map = rules.optimized_fastbe_cluster.output.cluster
    output:
        ccf = "results/{tumors}/{tumors}_ccf.csv"
    resources:
        mem_mb = 3000,
        runtime = "2m",
        slurm_partition = "normal"
    conda:
        "../envs/scripts.yaml"
    message:
        "Calculating the cancer cell fraction (CCF) of each clone in each sample given a variant allele frequency (VAF) matrix and a clone-to-mutation map from fastBE." 
    shell:
        "python3 scripts/get_ccf.py {input.clone_map} {input.matrix} {output.ccf}"


# Plot phylogenetic tree :)

rule plot_tree:
    input:
        clone_map = rules.optimized_fastbe_cluster.output.cluster,
        adjacency_list = rules.fastbe_search.output.tree
    output:
        phylogenetic_tree = "results/{tumors}/tree/clonal_evolution.png"
    resources:
        mem_mb = 3000,
        runtime = "2m",
        slurm_partition = "normal"
    conda:
        "../envs/scripts.yaml"
    message:
        "Plotting a clonal tree."
    shell:
        "mkdir -p results/{wildcards.tumors}/tree && "
        "python3 scripts/plot_clone.py --cluster {input.clone_map} {input.adjacency_list} -o {output.phylogenetic_tree}"

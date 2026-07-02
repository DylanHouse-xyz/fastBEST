
# Plot phylogenetic tree :)

rule plot_tree:
    input:
        clone_map = rules.optimized_fastbe_cluster.output.cluster,
        adjacency_list = rules.fastbe_search.output.tree
    output:
        phylogenetic_tree = directory("results/{tumors}/tree/"),
        history = "results/{tumors}/tree/parents.csv"
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
        "python3 scripts/plot_clone.py {input.adjacency_list}--cluster {input.clone_map} -o {output.phylogenetic_tree}"


# Get CCF, plot proportion, plot phylogenetic tree, and plot fishplot will
# be here.

rule get_ccf:
    input:
        matrix = rules.append_column.output.af_matrix,
        labels = rules.vcf_converter.output.labeled_matrix,
        clone_map = rules.optimized_fastbe_cluster.output.cluster,
        parent_history = rules.plot_tree.output.history,
        manifest = config["manifest"]
    output:
        ccf = "results/{tumors}/{tumors}_ccf.csv",
        raw_ccf = "results/{tumors}/{tumors}_ccf-raw.csv"
    resources:
        mem_mb = 3000,
        runtime = "2m",
        slurm_partition = "normal"
    conda:
        "../envs/scripts.yaml"
    message:
        "Calculating the cancer cell fraction (CCF) of each clone in each sample given a variant allele frequency (VAF) matrix and a clone-to-mutation map from fastBE." 
    shell:
        "python3 scripts/get_ccf.py --clusters_file {input.clone_map} --vaf {input.matrix} --labels {input.labels }  --parents {input.parent_history} --output {output.ccf} --manifest {input.manifest}"



rule fishplot:
    input:
        parents = rules.plot_tree.output.history,
        raw_ccf = rules.get_ccf.output.raw_ccf,
        clusters = rules.optimized_fastbe_cluster.output.cluster,
    output:
        fishplot = "results/{tumors}/{tumors}_fishplot.png"
    conda:
        "../envs/fishplot.yaml"
    message:
        "Generating a fishplot to map out clonal evolutionary dynamics"
    shell:
        "python3 scripts/generate_fishplots.py --ccf-file {input.ccf} --parent-file {input.parents} --clusters-file {input.clusters} --output {output.fishplot}"

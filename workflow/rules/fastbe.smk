# Helper function to read the optimal number of clones after the find_kneedle_point checkpoint
def get_optimal_clones(wildcards):
    checkpoint_output = checkpoints.find_kneedle_point.get(**wildcards).output.optimal_clones
    with open(checkpoint_output, "r") as f:
        k_value = f.read().strip()
    return k_value

# Searches for the most optimal adjacency list for a given frequency matrix
rule fastbe_search:
    input:
        af_matrix = rules.vcf_converter.output.matrix
    output:
        tree = "results/{tumors}/fastbe/{tumors}_tree.txt",
        meta = "results/{tumors}/fastbe/{tumors}_results.json"
    params:
        name = "results/{tumors}/fastbe/{tumors}"
    log:
        "logs/{tumors}/fastbe_search.log"
    conda:
        "envs/fastbe.yaml"
    threads: 32
    resources:
        mem_mb = 24000
    shell:
        "fastbe search {input.tree} -o {params.name} -f 1"

# Initial cluster of variants with arbitrary clones
rule initial_fastbe_cluster:
    input:
        af_matrix = rules.vcf_converter.output.matrix,
        tree = rules.fastbe_search.output.tree
    output:
        meta_file = "results/{tumors}/fastbe/initial_cluster/initial_clustering_results.json"
    params:
        meta_dir = directory("results/{tumors}/fastbe/initial_cluster")
    conda:
        "envs/fastbe.yaml"
    shell:
        "mkdir -p results/{tumors}/fastbe/initial_cluster && "
        "fastbe cluster -k 6 -o {params.meta_dir}/initial {input.tree} {input.af_matrix}"

# Find optimal number of clones using the kneedle algorithm
checkpoint find_kneedle_point:
    input:
        meta = rules.initial_fastbe_cluster.output.meta_file
    output:
        optimal_clones = "results/{tumors}/fastbe/optimal_clones.txt"
    conda:
        "envs/scripts.yaml"
    shell:
        "python3 find_mutation.py -i {input.meta} -o {output.optimal_clones}"

rule optimized_fastbe_cluster:
    input:
        matrix = rules.vcf_converter.output.matrix,
        tree = rules.fastbe_search.output.tree,
        k_file = lambda wildcards: checkpoints.find_kneedle_point.get(**wildcards).output.optimal_clones
    output:
        cluster = "results/{tumors}/fastbe/fastbe_optimized_k_clustering.csv",
        meta = "results/{tumors}/fastbe/fastbe_optimized_k_clustering_results.json"
    params:
        k = lambda wildcards: get_optimal_clones(wildcards)
    conda:
        "envs/fastbe.yaml"
    shell:
        "fastbe cluster -k {params.k} -o results/{tumors}/fastbe/fastbe_optimized_k {input.tree} {input.matrix}"

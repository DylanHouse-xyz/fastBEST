# Entire workflow from splitting intervals to Mutect2 variant calling and filtering: followed by
# phylogenetic reconstruction with fastBE.

# Ensure to run snakemake with snakemake -s Reconstruction-workflow


configfile: "config/config.yaml"
configfile: "config/samples.yaml"


INTERVAL_SHARD_COUNT = 24
INTERVAL_SHARD_IDS = [str(i).zfill(4) for i in range(INTERVAL_SHARD_COUNT)]

include: "rules/mutect2.smk"
include: "rules/filter.smk"
include: "rules/get_variants.smk"
include: "rules/fastbe.smk"
include: "rules/analysis.smk"

rule all:
    input:
        expand("results/{tumors}/mutect_merged.stats", tumors = config["samples"]),
        expand("results/{tumors}/read_orientation_model.tar.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtered_all.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtering_stats.tsv", tumors = config["samples"]),
        expand("results/{tumors}/pass_variants.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/af_matrix.csv", tumors = config["samples"]),
        expand("results/{tumors}/fastbe/fastbe_optimized_k_clustering.csv", tumors = config["samples"]),
        expand("results/{tumors}/fastbe/fastbe_optimized_k_clustering_results.json", tumors = config["samples"]),
        expand("results/{tumors}/{tumors}_ccf.csv"),
        expand("scripts/clonal_lineage_trees.png")

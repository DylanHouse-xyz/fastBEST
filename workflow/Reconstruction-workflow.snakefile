# Entire workflow from splitting intervals to Mutect2 variant calling and filtering: followed by
# phylogenetic reconstruction with fastBE.

# Mutect2 workflow modified from https://github.com/GavinHaLab/mutect2_snakemake

# Ensure to run snakemake with snakemake -s Reconstruction-workflow


configfile: "config/config.yaml"
configfile: "config/samples.yaml"


CHROMOSOMES = config["chromosomes"]


rule all:
    input:
        expand("results/{tumors}/mutect_merged.stats", tumors = config["samples"]),
        expand("results/{tumors}/read_orientation_model.tar.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtered_all.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtering_stats.tsv", tumors = config["samples"]),
        expand("results/{tumors}/pass_variants.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/af_matrix.csv", tumors = config["samples"]),



# Mutect2 somatic variant call rule to output vcf, stats, and artifacts split per chromosome
rule Mutect2:
    input:
        
        tumor_bam = lambda wildcards: config["samples"][wildcards.tumor][0],
        normal_bam = lambda wildcards: config["samples"][wildcards.tumor][2],
    output:
        vcf = temp("results/{tumor}/unfiltered_{chromosomes}.vcf.gz"),
        tar = temp("results/{tumor}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = temp("results/{tumor}/unfiltered_{chromosomes}.vcf.gz.stats")
    params:
        ref = config["ref_genome"],
        germ = config["germline_resource"],
    log:
        "logs/mutect2/{tumor}_{chromosomes}_mutect2.txt",
    shell:
        "(gatk Mutect2 "
        "-R {params.ref} "
        "-I {input.tumor_bam} "
        "-I {input.normal_bam} "
        "--germline-resource {params.germ} "
        "--intervals {wildcards.chromosomes} "
        "--f1r2-tar-gz {output.tar} "
        "--output {output.vcf}) 2> {log}"


# Merge all vcf.gz.stats
rule merge_mutect_stats:
    input:
        stats = expand("results/{{tumors}}/unfiltered_{chromosome}.vcf.gz.stats", chromosome=CHROMOSOMES),
    output:
        protected("results/{tumors}/mutect_merged.stats"),
    log:
        "logs/merge_mutect_stats/{tumors}_merge_mutect_stats.txt",
    params:
        stats_flags = lambda wildcards, input: " ".join([f"-stats {f}" for f in input.stats]),
    shell:
        "(gatk MergeMutectStats "
        "{params.stats_flags} "
        "-O {output}) 2> {log}"

rule learn_read_orientation_model:
    input:
        tar = expand("results/{{tumors}}/unfiltered_{chromosome}_f1r2.tar.gz", chromosome=CHROMOSOMES),
    output:
        protected("results/{tumors}/read_orientation_model.tar.gz"),
    log:
        "logs/learn_read_orientation_model/{tumors}_learn_read_orientation_model.txt",
    params:
        tar_flags = lambda wildcards, input: " ".join([f"-I {f}" for f in input.tar]),
    shell:
        "(gatk LearnReadOrientationModel "
        "{params.tar_flags} "
        "-O {output}) 2> {log}"

rule Merge_Results:
    input:
        vcfs = expand("results/{{tumors}}/unfiltered_{chromosome}.vcf.gz", chromosome=CHROMOSOMES),
    output:
        vcf = temp("results/{tumors}/final_unfiltered_merged.vcf.gz"),
    params:
        vcf_list = lambda wildcards, input: " ".join([f"-I {v}" for v in input.vcfs]),
    shell:
        "gatk MergeVcfs {params.vcf_list} -O {output.vcf}"

rule filter_mutect_calls:
    input:
        unfiltered_vcf = "results/{tumors}/final_unfiltered_merged.vcf.gz",
        read_orientation_model = "results/{tumors}/read_orientation_model.tar.gz",
        mutect_stats = "results/{tumors}/mutect_merged.stats",
    output:
        filtered_vcf = protected("results/{tumors}/filtered_all.vcf.gz"),
        filtering_stats = protected("results/{tumors}/filtering_stats.tsv"),
    params:
        reference_genome = config["ref_genome"],
    log:
        "logs/filter_mutect_calls/{tumors}_filter_mutect_calls.txt",
    shell:
        "(gatk FilterMutectCalls "
        "-R {params.reference_genome} "
        "-V {input.unfiltered_vcf} "
        "--ob-priors {input.read_orientation_model} "
        "--stats {input.mutect_stats} "
        "--filtering-stats {output.filtering_stats} "
        "-O {output.filtered_vcf}) 2> {log}"

rule extract_pass_variants:
    input:
        filtered_vcf = "results/{tumors}/filtered_all.vcf.gz",
    output:
        pass_variants_vcf = protected("results/{tumors}/pass_variants.vcf.gz"),
    params:
        reference_genome = config["ref_genome"],
    log:
        "logs/extract_pass_variants/{tumors}_extract_pass_variants.txt",
    shell:
        "(gatk SelectVariants "
        "-R {params.reference_genome} "
        "-V {input.filtered_vcf} "
        "--exclude-filtered "
        "--create-output-variant-index "
        "-O {output.pass_variants_vcf}) 2> {log}"


rule vcf_converter:
        input:
                vcf="results/{tumors}/pass_variants.vcf.gz",
        output:
                protected("results/{tumors}/af_matrix,csv"),
        log:
                "logs/{tumors}/vcf_converter.log",
        script:
                "python scripts/vcfconverter.py -i {{input}} -o {{output}}"



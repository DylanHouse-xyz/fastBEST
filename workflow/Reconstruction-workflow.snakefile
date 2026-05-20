# Entire workflow from splitting intervals to Mutect2 variant calling and filtering: followed by
# phylogenetic reconstruction with fastBE.

# Ensure to run snakemake with snakemake -s Reconstruction-workflow


configfile: "config/config.yaml"
configfile: "config/samples.yaml"


INTERVAL_SHARD_COUNT = 24
INTERVAL_SHARD_IDS = [str(i).zfill(4) for i in range(INTERVAL_SHARD_COUNT)]


rule all:
    input:
        expand("results/{tumors}/mutect_merged.stats", tumors = config["samples"]),
        expand("results/{tumors}/read_orientation_model.tar.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtered_all.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtering_stats.tsv", tumors = config["samples"]),
        expand("results/{tumors}/pass_variants.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/af_matrix.csv", tumors = config["samples"]),




rule split_intervals:
    output:
        expand("results/intervals/{scatter}-scattered.interval_list", scatter=INTERVAL_SHARD_IDS),
    params:
        intervals = config["intervals"],
        reference = config["ref_genome"],
        scatter_count = INTERVAL_SHARD_COUNT
    shell:
        "gatk SplitIntervals "
        "-R {params.reference} "
        "-L {params.intervals} "
        "--scatter-count {params.scatter_count} "
        "-O results/intervals"


# Mutect2 somatic variant call rule to output vcf, stats, and artifacts split per interval shard
rule Mutect2:
    input:
        tumor_bam = lambda wildcards: config["samples"][wildcards.tumor][0],
        normal_bam = lambda wildcards: config["samples"][wildcards.tumor][2],
        intervals = "results/intervals/{scatter}-scattered.interval_list"
    output:
        vcf = temp("results/{tumor}/unfiltered_{scatter}.vcf.gz"),
        tar = temp("results/{tumor}/unfiltered_{scatter}_f1r2.tar.gz"),
        stats = temp("results/{tumor}/unfiltered_{scatter}.vcf.gz.stats")
    threads: 4
    resources: 
        mem_mb = 24000
    params:
        ref = config["ref_genome"],
        germ = config["germline_resource"],
        tumor_input = lambda wildcards, input: " ".join([f"-I {b}" for b in input.tumor_bam]),
        normal_name = lambda wildcards, input: config["samples"][wildcards.tumor][1]
    log:
        "logs/mutect2/{tumor}_{scatter}_mutect2.txt",
    shell:
        "(gatk Mutect2 "
        "-R {params.ref} "
        "{params.tumor_input} "
        "-I {input.normal_bam} "
        "-normal {params.normal_name} "
        "--germline-resource {params.germ} "
        "-L {input.intervals} "
        "--f1r2-tar-gz {output.tar} "
        "--native-pair-hmm-threads {threads} "
        "--java-options '-Xmx16g -XX:+UseParallelGC' "
        "--output {output.vcf}) 2> {log}"


# Merge all vcf.gz.stats
rule merge_mutect_stats:
    input:
        stats = expand("results/{{tumors}}/unfiltered_{scatter}.vcf.gz.stats", scatter=INTERVAL_SHARD_IDS),
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
        tar = expand("results/{{tumors}}/unfiltered_{scatter}_f1r2.tar.gz", scatter=INTERVAL_SHARD_IDS),
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
        vcfs = expand("results/{{tumors}}/unfiltered_{scatter}.vcf.gz", scatter=INTERVAL_SHARD_IDS),
    output:
        vcf = temp("results/{tumors}/final_unfiltered_merged.vcf.gz"),
    params:
        vcf_list = lambda wildcards, input: " ".join([f"-I {v}" for v in input.vcfs]),
    shell:
        "gatk MergeVcfs {params.vcf_list} -O {output.vcf}"

rule get_pileup_summaries:
    input:
        tumor_bam = lambda wildcards: config["samples"][wildcards.tumors][0]
    output:
        protected("results/{tumors}/pileup_summaries.table")
    params:
        known_polymorphic_sites = config["known_polymorphic_sites"],
        tumors = lambda wildcards, input: " ".join([f"-I {b}" for b in input.tumor_bam]),
    log:
        "logs/get_pileup_summaries/{tumors}_get_pileup_summaries.txt"
    shell:
        "(gatk GetPileupSummaries "
        "{params.tumors} "
        "-V {params.known_polymprhic_sites} "
        "-L {params.known_polymorphic_sites} "
        "-O {output}) 2> {log}"

rule get_pileup_summary_normal:
    input:
        lambda wildcards: config["samples"][wildcards.tumors][2],
    output:
        protected("results/{tumors}/normal_pileup_summaries.table"),
    params:
        known_polymorphic_sites = config["known_polymorphic_sites"],
        tumor = lambda wildcards, input: config["samples"][wildcards.tumors][2] 
    log:
        "log/get_pileup_summary_normal/{tumors}_get_pileup_summaries.txt"
    shell:
        "(gatk GetPileupSummaries "
        "-I {params.tumor} "
        "-V {params.known_polymorphic_sites} "
        "-L {params.known_polymorphic_sites} "
        "-O {output}) 2> {log}"

rule calculate_contamination:
    input:
        summary_table = "results/{tumors}/pileup_summaries.table",
        normal_summary = "resuts/{tumors}/normal_pileup_summaries.table",
    output:
        segments_table = protected("results/{tumors}/segments.table"),
        contamination_table = protected("results/{tumors}/contamination.table")
    log:
        "logs/calculate_contamination/{tumors}_calculate_contamination.txt",
    shell:
        "(gatk CalculateContamination "
        "-I {input.summary_table} "
        "-I {input.normal_summary} "
        "-tumor-segmentation {otput.segments_table} "
        "-matched {input.normal_summary} "
        "-O {output.contamination_table}) 2> {log}"

rule filter_mutect_calls:
    input:
        unfiltered_vcf = "results/{tumors}/final_unfiltered_merged.vcf.gz",
        read_orientation_model = "results/{tumors}/read_orientation_model.tar.gz",
        mutect_stats = "results/{tumors}/mutect_merged.stats",
        segments_table = "results/{tumors}/segments.table",
        contamination_table = "results/{tumors}/contamination.table",
        intervals = config["intervals"],
    output:
        filtered_vcf = protected("results/{tumors}/filtered_all.vcf.gz"),
        filtering_stats = protected("results/{tumors}/filtering_stats.tsv"),
    params:
        reference_genome = config["ref_genome"],
        max_events_in_region = 1,
    log:
        "logs/filter_mutect_calls/{tumors}_filter_mutect_calls.txt",
    shell:
        "(gatk FilterMutectCalls "
        "-R {params.reference_genome} "
        "-V {input.unfiltered_vcf} "
        "--tumor-segmentation {input.segments_table} "
        "--contamination-table {input.contamination_table} "
        "--ob-priors {input.read_orientation_model} "
        "--max-events-in-region {params.max_events_in_region} "
        "--stats {input.mutect_stats} "
        "--filtering-stats {output.filtering_stats} "
        "-O {output.filtered_vcf}) 2> {log}"

rule filter_ffpe_artifacts:
    input:
        filtered_vcf = "results/{tumors}/filtered_all.vcf.gz",
    output:
        ffpe_filtered_vcf = protected("results/{tumors}/filtered_no_ffpe_artifacts.vcf.gz"),
    log:
        "logs/filter_ffpe_artifacts/{tumors}_filter_ffpe_artifacts.txt",
    shell:
        """ (bcftools view -e '((REF="C" && ALT="T") || (REF="G" && ALT="A"))' {input.filtered_vcf} -O z -o {output.ffpe_filtered_vcf} && bcftools index -t {output.ffpe_filtered_vcf}) > {log} 2>&1"""

rule extract_pass_variants:
    input:
        filtered_vcf = "results/{tumors}/filtered_no_ffpe_artifacts.vcf.gz",
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
        "--select-type-to-include SNP "
        "--create-output-variant-index "
        "-O {output.pass_variants_vcf}) 2> {log}"


rule vcf_converter:
        input:
                vcf="results/{tumors}/pass_variants.vcf.gz",
        output:
                matrix=protected("results/{tumors}/af_matrix.csv"),
        log:
                "logs/{tumors}/vcf_converter.log",
        shell:
                "python scripts/vcfconverter.py -i {input.vcf} -o {output.matrix}"


rule fastbe_search:
    input:
        af_matrix = "results/{tumors}/af_matrix.csv"
    output:
        tree = protected("results/{tumors}/fastbe/{tumors}")
        metadata = protected("results/{tumors}/fastbe/{tumors}.json")
    log:
        "logs/{tumors}/fastbe_search.log"
    threads = 32
    
    shell:
        "fastbe search {input.af_matrix} -o {output.tree} -t 32 -f1" 




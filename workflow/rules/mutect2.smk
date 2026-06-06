

rule split_intervals:
    output:
        expand("results/intervals/{scatter}-scattered.interval_list", scatter=INTERVAL_SHARD_IDS),
    params:
        intervals = config["intervals"],
        reference = config["ref_genome"],
        scatter_count = INTERVAL_SHARD_COUNT
    resources:
        mem_mb = 2048,
        runtime="1m",
        slurm_partition="MSC"
    conda:
        "../envs/mutect2.yaml"
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
        mem_mb = 24000,
        runtime="24h",
        slurm_partition="highmem"
    params:
        ref = config["ref_genome"],
        germ = config["germline_resource"],
        tumor_input = lambda wildcards, input: " ".join([f"-I {b}" for b in input.tumor_bam]),
        normal_name = lambda wildcards, input: config["samples"][wildcards.tumor][1]
    conda:
        "../envs/mutect2.yaml"
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
    resources:
        mem_mb = 2048,
        runtime = "3m",
        slurm_partition="MSC"
    conda:
        "../envs/mutect2.yaml"
    shell:
        "(gatk MergeMutectStats "
        "{params.stats_flags} "
        "-O {output}) 2> {log}"

rule Merge_Results:
    input:
        vcfs = expand("results/{{tumors}}/unfiltered_{scatter}.vcf.gz", scatter=INTERVAL_SHARD_IDS),
    output:
        vcf = temp("results/{tumors}/final_unfiltered_merged.vcf.gz"),
    params:
        vcf_list = lambda wildcards, input: " ".join([f"-I {v}" for v in input.vcfs]),
    resources:
        mem_mb = 2048,
        runtime = "3m",
        slurm_partition="MSC"
    conda:
        "../envs/mutect2.yaml"
    message:
        "Merging split VCF files."
    shell:
        "gatk MergeVcfs {params.vcf_list} -O {output.vcf}"



rule learn_read_orientation_model:
    input:
        tar = expand("results/{{tumors}}/unfiltered_{scatter}_f1r2.tar.gz", scatter=INTERVAL_SHARD_IDS),
    output:
        protected("results/{tumors}/read_orientation_model.tar.gz"),
    log:
        "logs/learn_read_orientation_model/{tumors}_learn_read_orientation_model.txt",
    params:
        tar_flags = lambda wildcards, input: " ".join([f"-I {f}" for f in input.tar]),
    conda:
        "envs/mutect2.yaml"
    shell:
        "(gatk LearnReadOrientationModel "
        "{params.tar_flags} "
        "-O {output}) 2> {log}"

rule get_pileup_summaries:
    input:
        tumor_bam = lambda wildcards: config["samples"][wildcards.tumors][0]
    output:
        protected("results/{tumors}/pileup_summaries.table")
    params:
        known_polymorphic_sites = config["known_polymorphic_sites"],
        tumors = lambda wildcards, input: " ".join([f"-I {b}" for b in input.tumor_bam]),
    conda:
        "envs/mutect2.yaml"
    log:
        "logs/get_pileup_summaries/{tumors}_get_pileup_summaries.txt"
    shell:
        "(gatk GetPileupSummaries "
        "{params.tumors} "
        "-V {params.known_polymorphic_sites} "
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
    conda:
       "envs/mutect2.yaml" 
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
        normal_summary = "results/{tumors}/normal_pileup_summaries.table",
    output:
        segments_table = protected("results/{tumors}/segments.table"),
        contamination_table = protected("results/{tumors}/contamination.table")
    conda:
        "envs/mutect2.yaml"
    log:
        "logs/calculate_contamination/{tumors}_calculate_contamination.txt",
    shell:
        "(gatk CalculateContamination "
        "-I {input.summary_table} "
        "-I {input.normal_summary} "
        "-tumor-segmentation {output.segments_table} "
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
    conda:
        "envs/mutect2.yaml"
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
    conda:
        "envs/bcftools.yaml"
    shell:
        """ (bcftools view -e '((REF="C" && ALT="T") || (REF="G" && ALT="A"))' {input.filtered_vcf} -O z -o {output.ffpe_filtered_vcf} && bcftools index -t {output.ffpe_filtered_vcf}) > {log} 2>&1"""

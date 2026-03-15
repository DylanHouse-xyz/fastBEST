rule mutect2:
    input:
        bam=expand("{bam_dir}/{tumor}",
            bam_dir=config['bam_dir'],
            tumor=config['samples']['tumors']),
        fasta=config['ref']['fasta'],
    output:
        vcf="results/{batch}.vcf.gz".format(batch=config['batch_name']),
    params:
        bams=lambda wc, input: " ".join(["-I " + b for b in input.bam]),
        normal=config['samples']['normal'],
        germline=config['ref']['germline'],
        intervals=config['ref']['intervals'],
    threads: 8
    resources:
        mem_mb=53248,
    log:
        "logs/mutect2_{batch}.log".format(batch=config['batch_name']),
   shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' Mutect2 "
        "-R {input.fasta} "
        "{params.bams} "
        "-normal {params.normal} "
        "--germline-resource {params.germline} "
        "-L {params.intervals} "
        "-O {output.vcf} "
        "--native-pair-hmm-threads {threads} "
        "&> {log}"

rule filter_mutect:
    input:
        vcf="results/{batch}.vcf.gz".format(batch=config['batch_name']),
        ref=config['ref']['fasta'],
    output:
        vcf="results/{batch}_filtered.vcf.gz".format(batch=config['batch_name']),
    log:
        "logs/filter_mutect_{batch}.log".format(batch=config['batch_name']),
    shell:
        "gatk FilterMutectCalls "
        "-R {input.ref} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        "&> {log}"

rule mutect2:
    input:
        map=os.path.join(config["bam_dir"], "{sample}.bam"),
        normal=os.path.join(config["bam_dir"], config["samples"]["normal"] + ".bam"),
        fasta=config['ref']['fasta'],
        intervals=config['ref']['intervals'],
        germline=config['ref']['germline']
    output:
        vcf="results/{sample}.vcf.gz"
    params:
        normal_name=config['samples']['normal']
    threads: 8
    resources:
        mem_mb=9216
    log:
        log_file="logs/mutect2_{sample}.log"
    wrapper:
        "v7.6.0/bio/gatk/mutect"

rule filter_mutect:
    input:
        vcf="results/{sample}.vcf.gz",
        ref=config['ref']['fasta']
    output:
        vcf="results/{sample}_filtered.vcf.gz"
    log:
        log_file="logs/filter_mutect_{sample}.log"
    wrapper:
        "v7.6.0/bio/gatk/filtermutectcalls"

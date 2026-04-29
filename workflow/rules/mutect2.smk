


# Mutect2 variant calling to output a single VCF file and artifact priors.
rule mutect2:
    input:
        bam=expand("{bam_dir}/{tumor}",
            bam_dir=config['bam_dir'],
            tumor=config['samples']['tumors']),
        fasta=config['ref']['fasta'],
    output:
        vcf="results/{batch}.vcf.gz".format(batch=config['batch_name']),
        f1r2="results/{batch}.f1r2.tar.gz".format(batch=config['batch_name'])
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
    wrapper:
        "v7.6.0/bio/gatk/mutect"


# LearnReadOrientationModel to generate maximum likelihood estimates from f1r2 file
rule OrientationModel:
    input:
        f1r2="results/{batch}.f1r2.tar.gz".format(batch=config['batch_name']),
    output:
        "results/{batch}_artifacts_prior.tar.gz".format(batch=config['batch_name']),
    resources:
        mem_mb=1024
    log:
        "logs/{batch}_learnreadorientationbias".format(batch=config['batch_name']),
    wrapper:
        "v7.6.0/bio/gatk/learnreadorientationmodel"

# Filter Mutect2
rule filter_mutect:
    input:
        vcf="results/{batch}.vcf.gz".format(batch=config['batch_name']),
        ref=config['ref']['fasta'],
        f1r2="results/{batch}_artifacts_prior.tar.gz".format(batch=config['batch_name']),

    output:
        vcf="results/{batch}_filtered.vcf.gz".format(batch=config['batch_name']),
    log:
        "logs/filter_mutect_{batch}.log".format(batch=config['batch_name']),
    wrapper:
        "v7.6.0/bio/gatk/filtermutectcalls"


# A wildcard constraint to limit num to single digit numbers. This would avoid other characters such as _ from raising syntax errors.
wildcard_constraints:
    num = r"\d+"


# A rule to run fastqc on all raw reads and generates a html report of the quality of said reads.
rule fastqc:
    input:
        files = lambda wildcards: config["reads"][wildcards.read][0][int(wildcards.num) - 1], # We use wildcard.num - 1 here because of python indexing.
    output:
        html ="results/fastqc/{read}_r{num}.html",
        zip = "results/fastqc/{read}_r{num}_fastqc.zip",
    threads: 1
    resources:
        mem_mb = 2024,
        slurm_partition = "MSC",
        runtime = "3m"
    message:
        "Running fastqc to check quality of sequence and generate a report"
    conda:
        "envs/alignment.yaml"
    wrapper:
        "v7.6.0/bio/fastqc"

# Note: As how Snakemake rules amd outputs work, this rule wil only run if these index files aren't present in the specified directory. Usually, for large genomes this could take ~30 minutes or less depending on the number of threads you soecify. 
rule bwa_index:
    input:
        genome = config["ref_genome"],
    output:
        idx = multiext("results/bwa/idx/hs38DH", ".0123", ".amb" ".ann", ".bwt.2bit.64", ".pac")
    threads: 32
    resources:
        mem_mb = 140000,
        slurm_partition = "highmem",
        runtime = "1h"
    message:
        "Generating index files for reference Genome"
    conda:
        "envs/alignment.yaml"
    log:
        "logs/bwa_index/genome.log"
    shell:
        "bwa-mem2 index -p hs38DH {input.genome}"



# Bwa-mem2 alignment of reads using samtools sorting. Input reads are in config/samples.yaml under reads. If bwa-index was run, the index files directly point to them.
rule bwa_mem2:
    input:
        reads= lambda wildcards: config["reads"][wildcards.read][0],
        idx = rules.bwa_index.output.idx,
    output:
        mapped = temp("results/mapped/{read}.bam"),
    params:
        extra=r"-R '@RG\tID:{read}\tSM:{read}'"
        sort='samtools',
        sort_order="coordinate"
    threads: 16,
    conda:
        "envs/alignment.yaml"
    message:
        "alignment with bwa mem2"
    wrapper:
        "v9.4.1/bio/bwa-mem2/mem"

# For reasons, this step where you use samtools collate, fixate, and sort again of the aligned reads are necessary for marking PCR duplicates prior to variant calling. 
rule samtools_fixmate:
    input:
        aln = rules.bwa_mem2.output.mapped,
    output:
        mc = temp("results/mapped/{read}_fix.bam"),
        sort = temp("results/mapped/{read}_sort.bam"),
        collate = temp("results/mapped/{read}_collate.bam")
    conda:
        "envs/alignment.yaml",
    shell:
        "samtools collate -o {output.collate} {input.aln} && "
        "samtools fixmate -m {output.collate} {output.mc} && "
        "samtools sort -o {output.sort} {output.mc}"




# Marking PCR duplicates to remove the faux statistical power these duplicates could infer. As stated in the documenation, this doesn't remove these duplicates, but rather marks them so they wont be used downstream.
rule samtools_markdup:
    input:
        aln = rules.samtools_fixmate.output.sort
    output:
        bam = "results/mapped/{read}.markdup.bam",
        idx = "results/mapped/{read}.markdup.bam.csi",
    params:
        extra="-c --no-PG -r",
    threads: 2,
    conda:
        "envs/alignment.yaml",
    message:
        "Marking PCR duplicates prior to analysis",
    wrapper:
        "v9.4.1/bio/samtools/markdup"




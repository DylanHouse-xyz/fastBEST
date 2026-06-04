
# A wildcard constraint to limit num to single digit numbers. This would avoid other characters such as _ from raising syntax errors.
wildcard_constraints:
    num = r"\d+"


# A rule to run fastqc on all raw reads and generates a html report of the quality of said reads.
rule fastqc:
    input:
        lambda wildcards: config["reads"][wildcards.read][0][int(wildcards.num) - 1], # We use wildcard.num - 1 here because of python indexing.
    output:
        html ="results/fastqc/{read}_r{num}.html",
        zip = "results/fastqc/{read}_r{num}_fastqc.zip",
    threads: 1
    resources:
        mem_mb = 2024
    message:
        "Running fastqc to check quality of sequence and generate a report"
    conda:
        "envs/alignment.yaml" 
    wrapper:
        "v7.6.0/bio/fastqc"


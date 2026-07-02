rule extract_pass_variants:
    input:
        filtered_vcf = "results/{tumors}/filtered_no_ffpe_artifacts.vcf.gz",
    output:
        pass_variants_vcf = protected("results/{tumors}/pass_variants.vcf.gz"),
    params:
        reference_genome = config["ref_genome"],
    resources:
        mem_mb = 1024,
        runtime = "1m",
        slurm_partition = "MSC"
    conda:
        "../envs/mutect2.yaml",
    log:
        "logs/extract_pass_variants/{tumors}_extract_pass_variants.txt",
    message:
        "Extracting all variants that passed the filter. Outputs a pass_variants vcf."
    shell:
        "(gatk SelectVariants "
        "-R {params.reference_genome} "
        "-V {input.filtered_vcf} "
        "--restrict-alleles-to BIALLELIC "
        "--exclude-filtered "
        "--select-type-to-include SNP "
        "--create-output-variant-index "
        "-O {output.pass_variants_vcf}) 2> {log}"


rule vcf_converter:
    input:
        vcf=rules.extract_pass_variants.output.pass_variants_vcf,
    output:
        matrix=protected("results/{tumors}/af_matrix.txt"),
        labeled_matrix = "results/{tumors}/labeled_matrix.txt,
    resources:
        mem_mb = 1024,
        runtime = "1m",
        slurm_partition = "MSC"
    log:
        "logs/{tumors}/vcf_converter.log",
    conda:
        "../envs/scripts.yaml",
    message:
        "Converting passed variants into a fastBE-ready format. Rows are samples and columns are distinct mutations."
    shell:
        "python3 scripts/vcfconverter.py -i {input.vcf} -o {output.matrix} -label {output.labeled_matrix}"

rule append_column:
    input:
        matrix = rules.vcf_converter.output.matrix,
    output:
        af_matrix = protected("results/{tumors}/af_matrix_root.txt"),
        labelled_root = "results/{tumors}/labelled_root_matrix.txt",
    resources:
        mem_mb = 1024,
        runtime = "1m",
        slurm_partition = "MSC"
    log:
        "logs/{tumors}/append_column.log"
    message:
        "Appending a column of 1 to the variant allele frequency matrix. This rule is used when the root is unknown."
    shell:
        """
        awk '{{print "1.0", $0}}' {input.matrix} > {output.af_matrix}
        """

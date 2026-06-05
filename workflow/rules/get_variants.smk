rule extract_pass_variants:
    input:
        filtered_vcf = "results/{tumors}/filtered_no_ffpe_artifacts.vcf.gz",
    output:
        pass_variants_vcf = protected("results/{tumors}/pass_variants.vcf.gz"),
    params:
        reference_genome = config["ref_genome"],
    conda:
        "envs/mutect2.yaml",
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
        matrix=protected("results/{tumors}/af_matrix.txt"),
    log:
        "logs/{tumors}/vcf_converter.log",
    conda:
        "envs/scripts.yaml",
    shell:
        "python3 scripts/vcfconverter.py -i {input.vcf} -o {output.matrix}"

rule append_column:
    input:
        matrix = rules.vcf_converter.input.vcf
    output:
        af_matrix = protected("results/{tumors}/af_matrix_root.txt")
    log:
        "logs/{tumors}/append_column.log"
    message:
        "Appending a column of 1 to the variant allele frequency matrix. This rule is used when the root is unknown."
    shell:
        """
        awk '{print "1.0", $0}' {input.matrix} > {output.af_matrix}
        """

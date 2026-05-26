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

rule vcf_converter:
	input:
		vcf="filtered/{sample}.vcf.gz",
	output:
		"results/af_matrix.csv",
	log:
		"logs/vcf_converter.log",
	script:
		"scripts/vcfconverter.py"

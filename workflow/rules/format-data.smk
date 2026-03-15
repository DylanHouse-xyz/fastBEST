rule vcf_converter:
        input:
                vcf="results/{batch}_filtered.vcf.gz".format(batch=config['batch_name']),
        output:
                "results/{batch}_af_matrix.csv".format(batch=config['batch_name']),
        log:
                "logs/{batch}_vcf_converter.log".format(batch=config['batch_name']),
        script:
                "scripts/vcfconverter.py"

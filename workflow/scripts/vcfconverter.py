
from cyvcf2 import VCF
import pandas as pd
import argparse


def extract_vcf_data(vcf_file):
    """
This function takes the input VCF file and creates a Pandas dataframe containing the allele frequences of each sample. The commented out variables are there for those who are interested in the layout of the dataframe. Additionally, this function changes any allele frequency less than 0.05 to 0.
    """
    vcf = VCF(vcf_file)



    allele_frequencies = []
    # variant_ids = []
    # sample_names = vcf.samples

    for v in vcf:
        if v.format("AF") is not None:
            allele_frequencies.append(v.format("AF").flatten())
            #variant_ids.append(f"{v.CHROM}:{v.POS}:{v.REF}>{v.ALT[0]}")

        else:
            allele_frequencies.append(v.gt_alt_depths / (v.gt_ref_depths + v.gt_alt_depths))
           # variant_ids.append(f"{v.CHROM}:{v.POS}:{v.REF}>{v.ALT[0]}")




    allele_frequencies = pd.DataFrame(allele_frequencies).transpose()
    #allele_frequencies.columns = variant_ids
    #allele_frequencies.index = sample_names


    allele_frequencies = allele_frequencies.where(allele_frequencies > 0.05, 0.0)


    return allele_frequencies


def fastbe_matrix(allele_frequencies, matrix):
    """
    Simple function to convert a Pandas dataframe into a format that fastBE can use.
    """

    allele_frequencies.to_csv(path_or_buf =matrix, sep = ' ', header = False, index = False,float_format='%.2f')







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts the allele frequencies from the  input Variant Calling Format (VCF) file to the required frequency matrix format per fastBE specifications. Mutations correspond to columns and rows correspond to unique samples. ')
    parser.add_argument('-i', help="input VCF file")
    parser.add_argument('-o', help="Your output file name")
    args = parser.parse_args()
    if args.i:
        vcf_matrix = extract_vcf_data(args.i)
        print("Converting vcf to matrix...")
    if args.o:
        fastbe_matrix(vcf_matrix, args.o)
        print("Your output file has been produced")

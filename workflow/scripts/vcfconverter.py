
from cyvcf2 import VCF
import pandas as pd
import argparse


def extract_vcf_data(vcf_file):
    """
This function takes the input VCF file and creates a Pandas dataframe containing the allele frequencies of each sample. The commented out variables are there for those who are interested in the layout of the dataframe. Additionally, this function changes any allele frequency less than 0.05 to 0.
    """
    vcf = VCF(vcf_file)



    allele_frequencies = []
    for v in vcf:
        if v.FILTER is None:
            if v.format("AF") is not None:
                allele_frequencies.append(v.format("AF").flatten())
            else:
                allele_frequencies.append(v.gt_alt_depths / (v.gt_ref_depths + v.gt_alt_depths))

    allele_frequencies = pd.DataFrame(allele_frequencies, columns = vcf.samples).transpose()

    allele_frequencies = allele_frequencies.where(allele_frequencies > 0.02, 0.0)
    #allele_frequencies = allele_frequencies.where(allele_frequencies < 0.50, 0.0)

    mutations_present = (allele_frequencies > 0).sum(axis=0)
    max_mutation = allele_frequencies.max(axis=0)
    filter_mutation = (mutations_present > 1) & (max_mutation < 0.50) | ((mutations_present == 1) & (max_mutation > 0.05) & (max_mutation < 0.50))
    allele_frequencies = allele_frequencies.loc[:, filter_mutation]

    junk = [col for col in allele_frequencies.columns if (allele_frequencies[col] > 0).all()]
    allele_frequencies = allele_frequencies.drop(columns = junk)
    #print(allele_frequencies)
    return allele_frequencies


def fastbe_matrix(allele_frequencies, matrix):
    """
    Simple function to convert a Pandas dataframe into a format that fastBE can use.
    """

    allele_frequencies.to_csv(path_or_buf =matrix, sep = ' ', header = False, index = False,float_format='%.2f')

def labeled_matrix(allele_frequencies, matrix):
    """
    Saves a copy of the matrix keeping both the sample row labels and mutation column headers.
    """
    allele_frequencies.to_csv(path_or_buf=matrix, sep=' ', header=True, index=True, float_format='%.2f')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts the allele frequencies from the  input Variant Calling Format (VCF) file to the required frequency matrix format per fastBE specifications. Mutations correspond to columns and rows correspond to unique samples. ')
    parser.add_argument('-i', help="input VCF file")
    parser.add_argument('-o', help="Your output file name")
    parser.add_argument('-label', help='A labelled version of the matrix')
    args = parser.parse_args()
    if args.i:
        vcf_matrix = extract_vcf_data(args.i)
        print("Converting vcf to matrix...")
    if args.o:
        fastbe_matrix(vcf_matrix, args.o)
        print("Your output file has been produced")

    if args.label:
        labeled_matrix(vcf_matrix, args.label)
        print(f"Labeled output file '{args.label}' has been produced.")

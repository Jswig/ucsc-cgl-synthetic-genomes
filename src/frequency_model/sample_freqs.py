import allel
import argparse
import numpy as np 
import pandas as pd
import os
import sys
import json

parser = argparse.ArgumentParser(
    description='Compute sample frequencies for each SNP from a VCF file'
)
parser.add_argument('input', help='Input VCF file')
parser.add_argument('output', help='Output JSON file')

def compute_sample_freqs(
    variants_vcf: str, 
    output: str,
):
    variants = (allel
        .vcf_to_dataframe(variants_vcf, fields=['POS', 'REF', 'ALT'])
        .drop(['ALT_2', 'ALT_3'], axis=1) # ALT_2, ALT_3 are always empty
    )
    genotypes = allel.read_vcf(variants_vcf, fields=['calldata/GT'])
    genotypes = genotypes['calldata/GT']
    haplo_1 = pd.DataFrame(genotypes[:,:,0])
    haplo_2 = pd.DataFrame(genotypes[:,:,1])
    
    num_samples = haplo_1.shape[1] * 2

    haplos = pd.concat([variants, haplo_1, haplo_2], axis=1)
    haplos = haplos[(haplos['REF'].str.len() == 1) & (haplos['ALT_1'].str.len() == 1)] # subset to only keep SNPs
    var_freqs = {}

    for pos in haplos['POS'].unique():
        sum_alt = 0
        pos_variants = haplos[haplos['POS'] == pos]
        pos_freqs = {}
        # frequencies of each variant at this position
        for row in pos_variants.itertuples(index=False):
            # each tuple is of form (index, POS, REF, ALT, ...)
            num_alt = np.sum(row[3:])
            sum_alt += num_alt 
            pos_freqs[row[2]] = num_alt / num_samples
        # frequency of the reference 
        pos_freqs[str(row[1])] =  float((num_samples - sum_alt) / num_samples)

        # add result to main dictionary
        var_freqs[int(pos)] = pos_freqs

    with open('reduced_sample.json', 'w') as fp:
        json.dump(var_freqs, fp)

    return 

if __name__ == '__main__':
    args=parser.parse_args()
    compute_sample_freqs(args.input, args.output)
import allel
import argparse
import numpy as np 
import pandas as pd
import json
from tqdm import tqdm

parser = argparse.ArgumentParser(
    description='Compute sample frequencies for each SNP from a VCF file'
)
parser.add_argument('input', help='Input VCF file')
parser.add_argument('output', help='Output JSON file')
parser.add_argument('--snp_only', type=bool, default=True, help='Keep only SNPs')

def compute_sample_freqs(
    variants_vcf: str, 
    output: str,
    snp_only: bool,
):
    variants = (allel
        .vcf_to_dataframe(variants_vcf, fields=['POS', 'REF', 'ALT'])
        .drop(['ALT_2', 'ALT_3'], axis=1) # ALT_2, ALT_3 are always empty
    )
    genotypes = allel.read_vcf(variants_vcf, fields=['calldata/GT'])
    genotypes = genotypes['calldata/GT']
    # scikit-allel reads missing values as -1
    haplo_1 = pd.DataFrame(genotypes[:,:,0]).replace(-1, 0)
    haplo_2 = pd.DataFrame(genotypes[:,:,1]).replace(-1, 0)
    
    num_samples = haplo_1.shape[1] * 2

    haplos = pd.concat([variants, haplo_1, haplo_2], axis=1)
    if snp_only:
        haplos = haplos[(haplos['REF'].str.len() == 1) & (haplos['ALT_1'].str.len() == 1)] # subset to only keep SNPs
    
    var_freqs = {}

    for pos in tqdm(haplos['POS'].unique()):
        sum_alt = 0
        pos_variants = haplos[haplos['POS'] == pos]
        freqs = []
        variants = []
        # frequencies of each variant at this position
        for row in pos_variants.itertuples(index=False):
            # each tuple is of form (POS, REF, ALT, ...)
            num_alt = np.sum(row[3:])
            sum_alt += num_alt 
            variants.append(row[2])
            freqs.append(num_alt / num_samples)
        # frequency of the reference 
        variants.append(row[1])
        freqs.append((num_samples - sum_alt) / num_samples)
        # add result to main dictionary
        var_freqs[int(pos)] = (variants, freqs)

    with open(output, 'w') as fp:
        json.dump(var_freqs, fp)

    return 

if __name__ == '__main__': 
    args=parser.parse_args()
    compute_sample_freqs(args.input, args.output, args.snp_only)
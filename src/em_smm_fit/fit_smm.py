import allel
import argparse
import numpy as np
import tqdm
from numba import jit 

parser = argparse.ArgumentParser(
    description='fit symmetric mixture model'
)
parser.add_argument('input', help='Input VCF file')

def fit_em_smm(variants_vcf: str):
    variants = (allel
		.vcf_to_dataframe(variants_vcf, fields=['POS', 'REF', 'ALT'])
		.drop(['ALT_2', 'ALT_3'], axis=1) # ALT_2, ALT_3 are always empty
	)
	genotypes = allel.read_vcf(variants_vcf, fields=['calldata/GT'])
	genotypes = genotypes['calldata/GT']


if __name__ == '__main__':
    args = parser.parse_args()

    fit_em_smm(args.input)
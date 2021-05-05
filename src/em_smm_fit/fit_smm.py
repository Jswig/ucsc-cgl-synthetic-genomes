import allel
import argparse
import multiprocessing
import numpy as np
import tqdm
from numba import jit 

parser = argparse.ArgumentParser(
	description='fit symmetric mixture model'
)
parser.add_argument('input', help='Input VCF file')
parser.add_argument(
	'--n_iterations', type=int, default=20, help="Number of iterations of EM algorithm"
)
parser.add_argument(
	'-K', type=int, default=10, help="Number of components of mixture model"
)

rng = np.random.default_rng(42)

def fit_em_smm(variants_vcf: str, n_iterations: int, K: int):
	variants = (allel
		.vcf_to_dataframe(variants_vcf, fields=['POS', 'REF', 'ALT'])
		.drop(['ALT_2', 'ALT_3'], axis=1) # ALT_2, ALT_3 are always empty
	)
	genotypes = allel.read_vcf(variants_vcf, fields=['calldata/GT'])
	genotypes = genotypes['calldata/GT']	
	# scikit-allel reads missing values as -1
	genotypes = np.where(genotypes == -1, 0, genotypes)
	haplo_1 = genotypes[:,:,0]
	haplo_2 = genotypes[:,:,1]

	# find locus with largest number of variants in sample
	max_n_variants = (variants
		.groupby('POS')
		.count()
		.sort_values(by='REF')
	)[-1]
	n_loci = len(genotypes)

	# em initialization
	groups = rng.integers(0,5, size=n_samples)
	group_ini = rng.random(size=6)
	group_probs = group_ini / np.sum(group_ini) # make these probability vectors
	variant_ini = rng.random(size=(K, n_loci, max_n_variants)) 
	variant_probs = variant_ini / np.sum(variant_ini, axis=2, keepdims=1)

	# TODO: EM loop
	for i in range(n_iterations):
		pass

if __name__ == '__main__':
	args = parser.parse_args()

	fit_em_smm(
		variants_vcf=args.input,
		n_iterations=args.n_iterations,
		K=args.K,
	)
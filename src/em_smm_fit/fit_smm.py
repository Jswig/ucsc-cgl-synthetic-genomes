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

@numba.jit
def _em_loop():
	pass

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
	
	# TODO better representation of haplotypes using ordinal encoding

	# find locus with largest number of variants in sample
	# add 1 to account for fact that we always have a reference
	max_n_variants = (variants
		.groupby('POS')
		.count()
		.sort_values(by='REF')
	)['REF'][-1] + 1 
	n_loci = genotypes.shape[0]
	n_samples = genotypes.shape[1] * 2 # we treat each chromosome independantly

	# em initialization
	groups_e_ini = rng.random(size=(K, n_samples))
	groups_e = group_e_ini / np.sum(groups_e_ini, axis=1, keepdims=1)
	group_ini = rng.random(size=6)
	group_probs = group_ini / np.sum(group_ini) # make this a probability vector
	variant_ini = rng.random(size=(K, n_loci, max_n_variants)) 
	# TODO: add step filtering this to correct number of variants
	variant_probs = variant_ini / np.sum(variant_ini, axis=2, keepdims=1) # make these  probability vectors

	# TODO: EM loop
	for i in range(n_iterations):
		# group belonging update
		for r in range(n_samples):
			# TODO group belonging update
			pass
		# group probability update
		# TODO 
		
		# variant probabilities update
		# TODO

if __name__ == '__main__':
	args = parser.parse_args()

	fit_em_smm(
		variants_vcf=args.input,
		n_iterations=args.n_iterations,
		K=args.K,
	)
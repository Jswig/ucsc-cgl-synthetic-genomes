import allel
import argparse
import multiprocessing
import numpy as np
import os

from numba import jit, prange
from tqdm import tqdm

parser = argparse.ArgumentParser(
	description='fit symmetric mixture model'
)
parser.add_argument('input', help='Input VCF file')
parser.add_argument('ouput', help='Output path')
parser.add_argument(
	'--n_iterations', type=int, default=20, help="Number of iterations of EM algorithm"
)
parser.add_argument(
	'-K', type=int, default=10, help="Number of components of mixture model"
)
rng = np.random.default_rng(42)

# @jit(nopython=False, parallel=True, fastmath=True)
def _em_loop(
	n_iterations: int,
	K: int,
	n_samples: int,
	n_loci: int,
	max_n_variants: int,
	group_e: np.ndarray,
	group_probs: np.ndarray, 
	variant_probs: np.ndarray, 
	haplotypese: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:	

	all_variant_cols = np.arange(n_loci)

	for i in range(n_iterations):
		# group expectation by sample update
		for r in range(n_samples): # this can be parallelized
			probs_alpha = [
				group_probs[alpha] * np.prod(
					variant_probs[alpha, all_variant_cols, haplotypes[r]]
				)
				for alpha in range(K)
			]
			group_e[r,:] = [
				probs_alpha[alpha] / np.sum(probs_alpha) 
				for alpha in range(K)
			]
		# NOTE if paralellizing, might need to be careful about assigning to group_e,
		# might trigger a race condition

		# group probability update
		group_probs = np.sum(groups_e, axis=0) / n_samples
		
		# variant probabilities update
		for alpha in range(K): # this can also be parallelized
			norm_ct = group_probs[alpha]*n_samples
			for k in range(n_loci):
				for i in range(max_n_variants):
					variant_samples = np.where(haplotypes[:,k] == i, 1, 0) # find samples with given variant
					group_probs[alpha, k, i] = np.sum(group_e[variant_samples, alpha]) / norm_ct
	return (group_probs, groups_e, variant_probs)

def fit_em_smm(
	variants_vcf: str, n_iterations: int, K: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

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
	
	# TODO better representation of haplotypes using ordinal encoding?
	# not sure if this is the best option

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
	group_e_ini = rng.random(size=(K, n_samples))
	group_e = group_e_ini / np.sum(groups_e_ini, axis=1, keepdims=1)
	group_ini = rng.random(size=6)
	group_probs = group_ini / np.sum(group_ini) # make this a probability vector
	variant_ini = rng.random(size=(K, n_loci, max_n_variants)) 
	# TODO: add step filtering this to correct number of variants
	variant_probs = variant_ini / np.sum(variant_ini, axis=2, keepdims=1) # make these  probability vectors

	# EM Loop
	return _em_loop(
		n_iterations,
		K,
		n_samples,
		n_loci,
		max_n_variants,
		group_e,
		group_probs,
		variant_probs,
		haplotypes,
	)
	
if __name__ == '__main__':
	args = parser.parse_args()

	(group_probs, groups_e, variant_probs) = fit_em_smm(
		variants_vcf=args.input,
		n_iterations=args.n_iterations,
		K=args.K,
	)
	np.save(os.path.join(args.output, 'group_probs.npy'), group_probs, allow_pickle=False)
	np.save(os.path.join(args.output, 'group_e.npy'), group_e, allow_pickle=False)
	np.save(os.path.join(args.output, 'variant_probs.npy'), variant_probs, allow_pickle=False)
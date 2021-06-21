import allel
import argparse
import numpy as np
import os
from numba import jit
from typing import Tuple

from src.lib.math import logsum_approx

parser = argparse.ArgumentParser(description='fit symmetric mixture model')
parser.add_argument('input', help='Input VCF file')
parser.add_argument('output', help='Output path')
parser.add_argument(
	'--n_iterations', type=int, default=20, help='Number of iterations of EM algorithm'
)
parser.add_argument('-K', type=int, default=10, help='Number of components of mixture model')
parser.add_argument(
	'--use_approx', dest='logsum_approx', action='store_true', help='Approximate log of sum with convexity lower bound in high-dimensional cases'
)
parser.add_argument(
	'--no_approx', dest='logsum_approx', action='store_false', help='Do no approximate log of a sum wiht convexity lower bound'
)
parser.set_defaults(logsum_approx=False)
parser.add_argument('--seed', type=int, default=42, help='Random seed')

@jit(nopython=True, fastmath=True)
def _em_loop(
	n_iterations: int,
	K: int,
	n_samples: int,
	n_loci: int,
	n_variants_pos: np.ndarray,
	group_e: np.ndarray,
	group_probs: np.ndarray, 
	variant_probs: np.ndarray, 
	haplotypes: np.ndarray,
	logsum_approx: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:	

	flattened_offset = np.arange(0, n_loci*2, 2)
 
	for i in range(n_iterations):
		print('Starting iteration ', i, ' /', n_iterations, '\n')
		for r in range(n_samples):
			log_variant_probs = np.log(variant_probs) # NOTE might cause problems for variant probs that are 0?
			log_probs_alpha = np.array([
				np.sum(log_variant_probs[alpha].flatten()[flattened_offset + haplotypes[r]]) 
				# workaround for not being able to use more than one advanced index in numba
				# normally would use variant_probs[alpha, np.arange(n_loci), haplos[r]]
				for alpha in range(K)
			]) 
			denominator = 0
			if logsum_approx: # use approximation for log of a sum
				denominator = logsum_approx(group_probs * log_probs_alpha)
			else:
				denominator = np.log(np.sum(group_probs * np.exp(log_probs_alpha)))     
			new_group_e = np.log(group_probs) + log_probs_alpha - denominator
			new_group_e = np.exp(new_group_e)
			if logsum_approx: # rescale to sum up to 1 if using convex approx
				new_group_e = new_group_e / np.sum(new_group_e)
			if np.isinf(np.max(new_group_e)):
				raise ValueError('Infinite expectation')
			else:
				group_e[r,:] = new_group_e

		group_probs = np.sum(group_e, axis=0) / n_samples
		print('Group probabilities: ', group_probs)
		# variant probabilities update
		for alpha in range(K):
			norm_ct = group_probs[alpha] * n_samples
			for k in range(n_loci):
				for i in range(n_variants_pos[k]): # ignore placeholders in array
					variant_samples = np.where(haplotypes[:,k] == i, True, False) # find samples with given variant
					new_prob = np.sum(group_e[variant_samples, alpha]) / norm_ct
					variant_probs[alpha, k, i] = max(new_prob, 10e-12) # prevent issues from taking log(0)
		print('Variant probabilities: ', variant_probs[:2, :10])    
	return (group_probs, group_e, variant_probs)

@jit(nopython=True)
def _encode_haplotypes(
 variants_pos: np.ndarray, haplos: np.ndarray, n_samples: int, n_loci: int,
) -> np.ndarray:
	haplos_encoded = np.full((n_samples, n_loci), 0)
	haplos_idx = 0 
	for k, pos in enumerate(np.unique(variants_pos)):
		n_variants = np.sum(variants_pos == pos)
		for idx in range(1, n_variants+1):
			res = np.where(
				haplos[:, haplos_idx] == 1, idx, haplos_encoded[:,k]
			)
			haplos_encoded[:,k] = res
			haplos_idx += 1
	return haplos_encoded

def fit_em_smm(
	variants_vcf: str, n_iterations: int, K: int, seed: int, logsum_approx: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
	haplos = np.hstack((haplo_1, haplo_2)).T

	n_variants_pos = (variants  # find number of variants by position
		.groupby('POS')			# add 1 to account for fact that we always have a reference
		.count()['REF']
		.values) + 1
	max_n_variants = np.sort(n_variants_pos)[-1]
	n_loci = len(variants['POS'].unique())
	n_samples = haplos.shape[0] 
	haplos = _encode_haplotypes(variants['POS'].values, haplos, n_samples, n_loci)

   	# em initialization
	rng = np.random.default_rng(seed)
	group_e_ini = rng.random(size=(n_samples, K))
	group_e = group_e_ini / np.sum(group_e_ini, axis=1, keepdims=1)
	group_probs = np.full(6, 1/K) # make this a probability vector
	variant_ini = rng.random(size=(K, n_loci, max_n_variants)) 
	variant_probs = variant_ini / np.sum(variant_ini, axis=2, keepdims=1)
	# TODO: add step filtering this to correct number of variants

	return _em_loop(
		n_iterations, K, n_samples,
		n_loci,
		n_variants_pos,
		group_e,
		group_probs,
		variant_probs,
		haplos,
		logsum_approx,
	)
	
if __name__ == '__main__':
	args = parser.parse_args()
	(group_probs, group_e, variant_probs) = fit_em_smm(
		variants_vcf=args.input,
		n_iterations=args.n_iterations,
		K=args.K,
		seed=args.seed,
		logsum_approx=args.logsum_approx,
	)
	np.save(os.path.join(args.output, 'group_probs.npy'), group_probs, allow_pickle=False)
	np.save(os.path.join(args.output, 'group_e.npy'), group_e, allow_pickle=False)
	np.save(os.path.join(args.output, 'variant_probs.npy'), variant_probs, allow_pickle=False)
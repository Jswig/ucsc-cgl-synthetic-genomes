import argparse 
import numpy as np 
import pandas as pd 

from os.path import join
from numba import jit 
from typing import Dict
from src.lib.vcf import samples_to_vcf

parser = argparse.ArgumentParser(
    description='Generate samples from mixture model'
)
parser.add_argument(
    'input', 
    help='Input .npy files parent directory (group_probsm group_e, variant_probs)'
)
parser.add_argument('output', help='Output VCF file')
parser.add_argument(
    '-n', '--n_samples', default=100, type=int, help='Number of sample genotypes to generate'
)
parser.add_argument('s', '--seed', default=42, type=int, help='Random seed')

def _generate_samples(
    n_samples: int,
    group_probs: np.ndarray,
    variant_probs: np.ndarray,
    seed: int,
) -> Dict[int, np.array]:
    samples = []
    
    rng = np.random.default_rng(seed)
    
    n_loci = variant_probs.shape[1]
    max_n_variants = variant_probs.shape[2]
    variant_codes = np.arange(max_n_variants)
    n_groups = len(group_probs)

    haplotypes = np.full((n_samples*2, n_loci), 0)

    groups = rng.choice(
        np.arange(0, len(group_probs)),
        size=n_samples*2,
        p=group_probs,
    )
    group_counts = [
        np.count_nonzero(groups == i)
        for i in range(K)
    ] 

    # iterate over groups, generating appropriate number of samples from 
    # each group
    for alpha in range(n_groups):
        for k in range(n_loci):
            pass

    # NOTE: this returns an array of samples. Maybe convert to a dict at a 
    # later point for compatibility with existing VCF writing code 

    # shuffle columns of resulting array (not really necessary, just prevents it 
    # from looking weird)


if __name__ == '__main__':
    args = parser.parse_args()

    group_e = np.load(join(args.input, 'group_e.npy'))
    group_probs = np.load(join(args.input, 'group_e.npy'))
    variant_probs = np.load(join(args.input, 'variant_probs.npy')) 
    
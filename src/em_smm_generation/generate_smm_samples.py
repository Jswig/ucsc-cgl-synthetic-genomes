import argparse 
import numpy as np 
import pandas as pd 

from numba import jit 
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
    '-n', '--n_samples', default=100, type=int, help='Number of samples to generate'
)

def generate_samples():
    pass

if __name__ == '__main__':
    args = parser.parse_args()
    pass
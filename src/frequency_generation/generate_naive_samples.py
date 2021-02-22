import allel 
import argparse
import json
import numpy as np 
import pandas as pd
from typing import List, Tuple

parser = argparse.ArgumentParser(
    description='Generate artificial samples using naive frequency model'
)
parser.add_argument('input', help='Sample variant frequencies in JSON format')
parser.add_argument('output', help='Output VCF file')
parser.add_argument(
    '-n', '--n_samples', default=100, type=int, help='Number of samples to generate'
)

def generate_samples(freqs: str, n_samples: int) -> Tuple[List[int], np.ndarray]:
    with open(freqs, 'r')  as fp:
        freqs_dict = json.load(fp)
    n_variants = len(freqs_dict)
    positions = freqs_dict.keys()
    samples = np.zeros([positions, n_samples])
    return (positions, samples)

def samples_to_vcf(output: str, samples: np.ndarray):
    with open(output, 'w') as vcf:
        pass

if __name__ == '__main__':
    args = parser.parse_args()
    samples = generate_samples(args.input, args.n_samples)
    samples_to_vcf(args.output, samples)
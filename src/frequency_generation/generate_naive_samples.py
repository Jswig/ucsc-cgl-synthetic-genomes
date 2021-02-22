import allel 
import argparse
import json
import numpy as np 
import pandas as pd
from tqdm import tqdm
from typing import List, Tuple

parser = argparse.ArgumentParser(
    description='Generate artificial samples using naive frequency model'
)
parser.add_argument('input', help='Sample variant frequencies in JSON format')
parser.add_argument('output', help='Output VCF file')
parser.add_argument(
    '-n', '--n_samples', default=100, type=int, help='Number of samples to generate'
)

def generate_samples(freqs_dict: dict, n_samples: int) ->  np.ndarray:
    samples = []
    for v in tqdm(freqs_dict.values()):
        samples.append(
            np.random.choice(v[0], size=n_samples*2, p=v[1])
        ) # n_samples*2 as we need two samples to get the genotype
    return samples

def samples_to_vcf(
    freqs_dict: dict, output: str, n_samples: int, chrom: int, samples: np.ndarray
):
    vcf_header = (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    )   
    ids = [f'NF{str(k).zfill(5)}\t' for k in range(n_samples)]
    vcf_header = vcf_header + ''.join(ids) + '\n'

    with open(output, 'w') as vcf:
        vcf.write(vcf_header)

        for pos in tqdm(freqs_dict.keys()):
            bases = freqs_dict[pos][0]
            ref = bases.pop()
            for variant in bases:
                vcf.write(f'17\t{pos}\t.\t{ref}\t{variant}\t.\t.\t.\tGT\t')



if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.input, 'r')  as f:
        freqs_dict = json.load(f)
    samples = generate_samples(freqs_dict, args.n_samples)
    print(samples[1])
    samples_to_vcf(freqs_dict, args.output, samples, args.n_samples)

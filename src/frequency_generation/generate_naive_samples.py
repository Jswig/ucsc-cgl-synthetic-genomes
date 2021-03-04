import allel 
import argparse
import json
import numpy as np 
import pandas as pd
from collections import deque
from tqdm import tqdm
from typing import Deque, Tuple

parser = argparse.ArgumentParser(
    description='Generate artificial samples using naive frequency model'
)
parser.add_argument('input', help='Sample variant frequencies in JSON format')
parser.add_argument('output', help='Output VCF file')
parser.add_argument(
    '-n', '--n_samples', default=100, type=int, help='Number of samples to generate'
)

def generate_samples(freqs_dict: dict, n_samples: int) ->  Deque[np.ndarray]:
    samples = deque()
    for v in tqdm(freqs_dict.values()):
        samples.append(
            np.random.choice(v[0], size=n_samples*2, p=v[1])
        ) # n_samples*2 as we need two samples to get the genotype
    # with open('output/samples_log.json', 'w') as log:
    #     l_samples = [list(sample) for sample in samples]
    #     json.dump(l_samples, log)  
    return samples

def samples_to_vcf(
    freqs_dict: dict, output: str, n_samples: int, samples: np.ndarray
):
    vcf_header = (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    )   
    ids = [f'NF{str(k).zfill(5)}\t' for k in range(n_samples-1)]
    ids.append(f'NF{str(n_samples-1).zfill(5)}')
    vcf_header = vcf_header + ''.join(ids) + '\n'

    with open(output, 'w') as vcf:
        vcf.write(vcf_header)

        for pos in tqdm(list(freqs_dict.keys())[:-1]):
            bases = freqs_dict[pos][0]
            ref = bases.pop()
            haplos = samples.popleft()
            haplo_1, haplo_2 = np.split(haplos, 2) # haplos has length 2*n_samples
            for variant in bases:
                vcf.write(f'13\t{pos}\t.\t{ref}\t{variant}\t.\t.\t.\tGT\t')
                haplo_1_has_var = np.where(haplo_1 == variant, 1, 0)
                haplo_2_has_var = np.where(haplo_2 == variant, 1, 0)
                genotypes = [
                    f'{haplo_1_has_var[k]}|{haplo_2_has_var[k]}\t'
                    for k in range(len(haplo_1) - 1)
                ]
                genotypes.append(
                    f'{haplo_1_has_var[len(haplo_1)-1]}|{haplo_2_has_var[len(haplo_1)-1]}'
                )
                vcf.write(''.join(genotypes))
                vcf.write('\n')
        

if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.input, 'r')  as f:
        freqs_dict = json.load(f)
    samples = generate_samples(freqs_dict, args.n_samples)
    samples_to_vcf(freqs_dict, args.output, args.n_samples, samples)

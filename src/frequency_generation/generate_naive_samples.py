import allel 
import argparse
import json
import numpy as np 
import pandas as pd
from collections import deque
from tqdm import tqdm
from typing import Dict, Tuple

parser = argparse.ArgumentParser(
	description='Generate artificial samples using naive frequency model'
)
parser.add_argument('input', help='Sample variant frequencies in JSON format')
parser.add_argument('output', help='Output VCF file')
parser.add_argument(
	'-n', '--n_samples', default=100, type=int, help='Number of samples to generate'
)

def generate_samples(
	freqs_dict: dict, 
	n_samples: int, 
	log: bool = False,
	log_path: str = None,
) ->  Dict[int, np.ndarray]:
	samples = {}
	for pos, v in tqdm(freqs_dict.items()):
		print(len(np.arange(0, len(v['freq']))))
		print(len(v['freq']))
		
		samples[pos] = np.random.choice(
			np.arange(0, len(v['freq'])), 
			size=n_samples*2, 	
			p=v['freq']
		)
	if log:
		with open(log_path, 'w') as log:
			l_samples = {k: list(v) for (k,v) in samples}
			json.dump(l_samples, log)  
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

		for pos, v in tqdm(freqs_dict.items()):
			refs = v['REF']
			variants = v['ALT']
			haplos = samples[pos]
			haplo_1, haplo_2 = np.split(haplos, 2) # haplos has length 2*n_samples 
			for i, (ref, var) in enumerate(zip(refs, variants)): 
				vcf.write(f'13\t{pos}\t.\t{ref}\t{var}\t.\t.\t.\tGT\t')
				haplo_1_has_var = np.where(haplo_1 == i, 1, 0)
				haplo_2_has_var = np.where(haplo_2 == i, 1, 0)
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
		# note: ensure keys in JSON are sorted 
		freqs_dict = json.load(f)
	samples = generate_samples(freqs_dict, args.n_samples)
	samples_to_vcf(freqs_dict, args.output, args.n_samples, samples)

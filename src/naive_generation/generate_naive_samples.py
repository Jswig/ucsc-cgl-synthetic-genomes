import allel 
import argparse
import json
import numpy as np 
from tqdm import tqdm
from typing import Dict, Tuple

from src.lib.vcf import samples_to_vcf

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
		samples[pos] = np.random.choice(
			np.arange(0, len(v['freq'])), 
			size=n_samples*2, 	
			p=v['freq']
		)
	if log:
		with open(log_path, 'w') as log:
			l_samples = {k: list(v) for (k,v) in samples	}
			json.dump(l_samples, log)  
	return samples

if __name__ == '__main__':
	args = parser.parse_args()
	with open(args.input, 'r')  as f:
		# note: ensure keys in JSON are sorted 
		freqs_dict = json.load(f)
	samples = generate_samples(freqs_dict, args.n_samples)
	samples_to_vcf(freqs_dict, args.output, args.n_samples, samples, '13')

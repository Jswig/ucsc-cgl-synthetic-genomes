import tqdm
import numpy as np
import json

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

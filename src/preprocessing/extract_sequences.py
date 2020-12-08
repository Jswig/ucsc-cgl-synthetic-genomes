import numpy as np 
import pandas as pd 
import allel

import Bio.SeqIO
import os.path 

def extract_from_vcf(
    reference_fasta: str, 
    variants_vcf: str, 
    sample_ids: str, 
    output: str,
    start: int, 
    end: int
):
    reference = Bio.SeqIO.read(reference_fasta, 'fasta')
    reference_seq = str(reference.seq)
    reference_gene = reference_seq[start:end]

    ids = pd.read_csv(sample_ids, header=None)
    ids = ids[0].tolist()

    variants = (allel
        .vcf_to_dataframe(variants_vcf, fields=['POS', 'REF', 'ALT'])
        .drop(['ALT_2', 'ALT_3'], axis=1) # ALT_2, ALT_3 are always empty
    )
    # shift POS to start from 0 in the gene's sequence
    variants['POS'] = variants['POS'] - start

    genotypes = allel.read_vcf(variants_vcf, fields=['calldata/GT'], samples=ids)
    genotypes = genotypes['calldata/GT']
    haplo_1 = pd.DataFrame(genotypes[:,:,0])
    haplo_2 = pd.DataFrame(genotypes[:,:,1])

    haplos = pd.concat([variants, haplo_1, haplo_2], axis=1)
    haplos = haplos[(haplos['REF'].str.len() == 1) & (haplos['ALT_1'].str.len() == 1)]
    haplos = haplos.reset_index() # reset index after filtering 

    # iterate over all haplotype columns, building the fulls sequence for each
    seqs = []
    for j in range(3, haplos.shape[1]):
        seq = list(reference_gene) # cannot modify a string
        for i in range(len(haplos)):    
            if haplos.iloc[i,j] == 1:
                alt = haplos.loc[i, 'ALT_1']
                seq[haplos.loc[i, 'POS']] = alt

        seq = pd.Series(seq, dtype='category')
        seqs.append(seq)
        print('Completed sample {}\n'.format(j))

    seqs_df = pd.DataFrame(seqs) # as feather is columnar prefer taking transpose after loading
    seqs_df.to_feather(output)


if __name__ == '__main__':

    if not os.path.isfile('data/interim/brca1_seqs.feather'):
        extract_from_vcf(
            reference_fasta='data/Homo_sapiens.GRCh38.dna.chromosome.17.fa', 
            variants_vcf='data/raw/brca1.vcf',
            sample_ids='data/interim/sample_names.txt',
            output='data/interim/brca1_seqs.feather',
            start=43044295,
            end=43170246,
        )

    if not os.path.isfile('data/interim/brca2_seqs.feather'):   
        extract_from_vcf(
            reference_fasta='data/Homo_sapiens.GRCh38.dna.chromosome.13.fa', 
            variants_vcf='data/raw/brca2.vcf',
            sample_ids='data/interim/sample_names.txt',
            output='data/interim/brca2_seqs.feather',
            start=32315086,
            end=32400267,
        )
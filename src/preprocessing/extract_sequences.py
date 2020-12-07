import numpy as np 
import pandas as pd 
import allel
from Bio import SeqIO


def extract_from_vcf(reference_fasta: str, variants_vcf: str, start: int, end: int):
    
    reference = list(SeqIO.parse(reference_fasta, "fasta"))
    variants = allel.read_vcf(variants_vcf, fields=["POS", "REF", "ALT"])

    reference_seq = str(reference.seq)
    reference_gene = reference_seq[start:end]


if __name__ == "__main__":

    extract_from_vcf(
        reference_fasta="data/Homo_sapiens.GRCh38.dna.chromosome.17.fa", 
        variants_vcf="data/raw/bcra1.vcf",
        start=43044295,
        end=43170246,
    )
    extract_from_vcf(
        reference_fasta="data/Homo_sapiens.GRCh38.dna.chromosome.13.fa", 
        variants_vcf="data/raw/bcra2.vcf",
        start=32315086,
        end=32400267,
    )
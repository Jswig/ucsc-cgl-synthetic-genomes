#!bin/bash

# downloads reference genome chromosomes from ensembl.org

SOURCE="ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna"
CHROM_13="Homo_sapiens.GRCh38.dna.chromosome.13.fa"
CHROM_17="Homo_sapiens.GRCh38.dna.chromosome.17.fa"

cd data

if [ ! -f $CHROM_13 ]
then
	wget $SOURCE/$CHROM_13.gz
	gunzip $CHROM_13
fi

if [ ! -f $CHROM_17 ]
then 
	wget $SOURCE/$CHROM_17.gz
	gunzip $CHROM_17.gz
fi 


#!bin/bash

if [ ! -f data/interim/brca_1_sample_names.txt ]
then
    bcftools query -l data/raw/brca1.vcf > data/interim/brca1_sample_names.txt
fi 

if [ ! -f data/interim/brca_2_sample_names.txt ]
then
    bcftools query -l  data/raw/brca2.vcf > data/interim/brca2_sample_names.txt
fi 
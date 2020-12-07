#!bin/bash

if [ ! -f data/interim/sample_names.txt ]
then
    bcftools query -l data/raw/brca1.vcf > data/interim/sample_names.txt
fi 

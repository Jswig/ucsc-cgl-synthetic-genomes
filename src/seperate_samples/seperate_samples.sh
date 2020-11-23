#!bin/sh

extract_samples() {
    for sample in `bcftools query -l $1` 
    do
        bcftools view -c1 -s $sample  -Ov
}

extract_samples $1 $2
extract_samples $3 $4
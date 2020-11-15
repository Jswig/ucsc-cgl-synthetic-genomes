#!bin/sh

extract_sample_names() {
    bcftools view -h $1 | 
        cut -s -f 10- | 
        sed '/^#/d' > "{$2}{$1}.txt"
}

extract_sample_names $1 $2
extract_sample_names $1 $3
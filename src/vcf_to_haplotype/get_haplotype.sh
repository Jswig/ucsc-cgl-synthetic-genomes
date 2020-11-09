#!bin/sh

praw="data/raw/"
pinterim="data/interim/"

haplo --database -i "${praw}f8_chr13_brca2_copd_hmb.vcf" -o "${pinterim}chr13_brca2_haplo.tsv"
haplo --database -i "${praw}f8_chr17_brca1_copd_hmb.vcf" -o "${pinterim}chr17_brca1_haplo.tsv"
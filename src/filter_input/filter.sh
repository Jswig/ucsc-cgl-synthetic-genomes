#!/bin/bash

cd data/raw

COMP_DATA = ($1).gz
bcftools view ($1) -Oz -o COMP_DATA
bcftools index COMP_DATA
bcf filter bcftools filter -r 13,17 COMP_DATA
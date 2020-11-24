#!bin/sh

usage() {
    cat <<HELP_USAGE
        $0 <gene1.vcf> <output_folder_1> <gene_2.vcf> <output_folder_2>
HELP_USAGE
}

extract_samples() {
    for sample in `bcftools query -l $1` 
    do
        bcftools view -c1 -s $sample -Ov $1 > "$2/$sample.vcf" 
    done
}

if [ $# == 0 ]
then
    usage
fi

extract_samples $1 $2
extract_samples $3 $4
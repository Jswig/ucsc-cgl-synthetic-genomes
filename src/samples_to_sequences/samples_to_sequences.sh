#!bin/sh

usage() {
    cat <<HELP_USAGE
        $0 <sample_folder_1> <sample_positions_1> <output_folder_1> <sample_folder_2> <sample_positions_2> <output_folder_2> <reference.fasta>
HELP_USAGE
}

extract_sequences() {
    for sample in $1/*.vcf
    do
        sample_name=`basename $sample .vcf`
        gatk FastaAlternateReferenceMaker \
            -R $4 \ 
            -O "$3$sample_name.fasta" \ 
            -L $2 \
            -V $sample 
    done 
}

if [ $# == 0]
then
    usage
    exit 1
fi

gatk CreateSequenceDictionary -R $7
samtools faidx $7

extract_sequences $1 $2 $3 $7
extract_sequences $4 $5 $6 $7



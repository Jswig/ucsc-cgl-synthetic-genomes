# 19-10-2020

BRCA1
BRCA2

- model individually 

What is biologically plausible ?

- Use HG-38 reference genome 
- use frequency of each variant 
- use edit distance  
- distribution of variance depending on population 
- think about ethnicity data (data should inform the discriminator) (population frequencies) 
- pathogenicity of a variant BRCAexchange.org
    - `clinical_significance_ENAGMA` benign, likely benign, pathogenic, unknown 
- mis-sense vs non-sense vs synonymous variants 
    - *synonymous*: mutation that creates the same amino acid 
    - *mis-sense*: create different amino acide
    - *non-sense*" early STOP codon
- location of mutations:
    - functional domain or not
- cohort frequency: compare to frequency in the specific dataset


# 16-11-2020


- Try with just SNPs (easy to code)
    - makes model evaluation easier 

- Try FASTARefMaker, if does not work just modify reference with SNPs from files

- It's biologically plausible to model each of a person's haplotypes independantly
    - just draw twice from the generator for each person

# 12-1-2020

We treat the haplotypes as independant. As 


Hyperparameter tuning for the GAN:
- tune using prior knowledge of expected parameters

# 7-12-2020

Model evaluation:

- Does it preserve intron/exon structure with splice sites?
    - MaxEntScan to detect number of splice sites
        - needs GATK to get data back in variant calling format
    - are functional domains conserved (just compare regions for equality)
        - also should be in variant calling format
    - edit distance 
    - are both variants dysfunctional (future work!)
    - ethnicity (is the variant similar to that from the population from that ethnicity) (future work!)

How many unique genes does the network generate?

The code can incorporate simple subsnititions easily (1 -> n or n->1)

Other variant calling tool: [hail](https://hail.is/docs/0.2/genetics/hail.genetics.Call.html)
- use relatedness tool in hail [link](https://hail.is/docs/0.2/methods/relatedness.html)
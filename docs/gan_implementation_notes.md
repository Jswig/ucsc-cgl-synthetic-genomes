# GAN implementation notes

## Data pre-processing

### Gene locations

BRCA 1 location
```
chr17:43,044,295-43,170,245
```

BRCA 2 location
```
chr13:32,315,086-32,400,266
```

reference genomes
```
GRCh38_latest_genomic.fna
```

### Storing data

Using default strings data is about 900kb/sample
Using `categorical` dtypes reduces this to 100kb
Can probably do better with numpy arrays? --> not really since `categorical` is 
really an integer
Advantages to keeping everything in pandas for one-hot encoding?
Even the csv should be reasonable

## Model

cf *Generating and Designing DNA with Deep Generative Models*

- Latent space sampled from standard normal distribution
- perform 5 discriminator updates for every generator update



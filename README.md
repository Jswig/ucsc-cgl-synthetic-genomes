# ucsc-cgl-synthetic-genomes

*Anders Poirel*

## Dependencies

- [Popper](https://popper.readthedocs.io/en/latest/)
- Docker

## Running all steps in the workflow

```sh
popper run -f wf.yml 
```

## Generating synthetic samples

### GAN

TODO

### Null model

<!-- #### With Popper

Change the steps in `wf.yml` to reflect the input and output files you wish to use,
as well as the user id 

```sh
popper run -f wf.yml frequency_model
popper run -f wf.yml frequency_generation
```

#### Without Popper -->

Run
```sh
docker run -it --mount ./ ./containers/allel/
```

In the container, 
```sh
python src/frequency_model/sample_frequencies.py \
    data/raw/brca2.vcf \
    output/models/sample_freqs.json &&

python src/frequency_generation/generate_naive_samples.py \
    generate_naive_samples.py \
    output/models/sample_freqs.json \
    output/samples/naive_samples.vcf \
    -n 10000 \
```
Omit the first command if `sample_freqs.json` already 
exists


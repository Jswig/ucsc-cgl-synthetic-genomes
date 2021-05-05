# ucsc-cgl-synthetic-genomes

*Anders Poirel*

## Dependencies

- [Popper](https://popper.readthedocs.io/en/latest/)
- Docker

## Container enginer configuration

If using Docker, and on a host OS where the files you wish to use 
require permission to access, edit `config.yml` and change `user` to have the appropiate uid.
Use the `--conf config.yml` with every call to Popper. 

*Note:* depending on the machine, additional configuration may be required to use the GPU.

## Running all steps in the workflow

```sh
popper run -f wf.yml 
```

## Generating synthetic samples

### GAN

TODO

### Naive (Null) model

#### With Popper

Change the steps in `wf.yml` to reflect the input and output files you wish to use

```sh
popper run -f wf.yml naive_fit
popper run -f wf.yml naive_generation
```

#### Without Popper

Run
```sh
docker run -it --mount ./ ./containers/allel/
```

In the container, 
```sh
python src/frequency_model/sample_frequencies.py \
    data/raw/brca2.vcf \
    output/models/brca2_freqs.json &&

python src/frequency_generation/generate_naive_samples.py \
    generate_naive_samples.py \
    -n 10000 \
    output/models/brca2_freqs.json \
    output/samples/naive_samples.vcf 
```
Omit the first command if `sample_freqs.json` already 
exists

## References

My term papers with explanations of relevant theoretical aspects are 
[here](https://github.com/Jswig/ucsc-gi-papers).


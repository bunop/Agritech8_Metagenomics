
# Agritech Spoke 8 - MetaGenomics

This project want to track steps to perform analysis of metagenomics data.

## Create the metadata files

In the `data` folder there are some `.xlsx` files that contain the metadata of the samples.
This file has been processed to fill the missing information using two ipython
notebooks: `scripts/metadata.ipynb` and `scripts/samplesheet.ipynb`. Please
execute these notebooks to generate the metadata files.

## Create the configuration files for nextflow

Simply call:

```bash
nf-core launch nf-core/ampliseq -r 2.11.0
```

and follow the instructions to create the configuration files. Two configuration files
has been created: `config/nf-params-bacteria.json` and `config/nf-params-fungi.json`
for *bacteria* and *fungi* respectively.

## Launch the Nextflow pipeline

To launch the pipeline, simply call:

```bash
nextflow run nf-core/ampliseq -r 2.11.0 -params-file config/nf-params-bacteria.json \
    -profile singularity -resume -c config/custom.config
```

for bacteria and

```bash
nextflow run nf-core/ampliseq -r 2.11.0 -params-file config/nf-params-fungi.json -profile singularity
```

for fungi.

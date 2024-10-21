
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

## Determine the best trimming parameters using figaro

[figaro](https://github.com/Zymo-Research/figaro) can analyze error rates in files and
suggest best parameter to be used with `dada2`: download the singularity image with

```bash
cd /home/core/nxf_singularity_cache/
singularity pull docker://quay.io/biocontainers/figaro:1.1.2--hdfd78af_0
cd -
```

Figaro requires all the files in the same folder. Get all bacteria files in `data/bacteria`
folder by calling the `scripts/bacteria.ipynb` Ipython notebook.
Figaro complains about reads of different lengths: use `seqkit` to ged rid of reads shorter
than 251 bp:

```bash
cd data
mkdir bacteria-fixed-size
cd bacteria
for fastq in $(ls *.gz); do seqkit seq -m 251 $fastq > ../bacteria-fixed-size/$fastq ; done
cd -
```

Next call figaro on bacteria folder:

```bash
singularity run -B /home/ /home/core/nxf_singularity_cache/figaro_1.1.2--hdfd78af_0.sif figaro.py \
    --outputDirectory test-figaro --ampliconLength 423 --forwardPrimerLength 17 --reversePrimerLength 21 \
    --inputDirectory data/bacteria-fixed-size
```

Select the first two results using figaro:

```bash
jq '.[0:2] | map({trimPosition, maxExpectedError, readRetentionPercent, score})' test-figaro/trimParameters.json
```

here's the output:

```json
[
  {
    "trimPosition": [
      250,
      231
    ],
    "maxExpectedError": [
      3,
      2
    ],
    "readRetentionPercent": 82.31,
    "score": 77.31462467136672
  },
  {
    "trimPosition": [
      241,
      240
    ],
    "maxExpectedError": [
      2,
      2
    ],
    "readRetentionPercent": 76.45,
    "score": 74.45350738362131
  }
]
```

First results has and higher retention, however it doesn't truncate R1 and has an higher *expected error*.
Chosen parameters are from the second result, and were passed to `config/nf-params-bacteria.json`

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

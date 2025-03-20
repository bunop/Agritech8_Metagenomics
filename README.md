
# Agritech Spoke 8 - MetaGenomics

This project want to track steps to perform analysis of metagenomics data.

## Create the metadata files

In the `data` folder there are some `.xlsx` files that contain the metadata of the samples.
This file has been processed to fill the missing information using two ipython
notebooks: `scripts/metadata.ipynb` and `scripts/samplesheet.ipynb`. Please
execute these notebooks to generate the metadata files.

### Update metadata

Metadata were then fixed to describe better samples. Now the current metadata
come from `metadata_bacteria_fix.xlsx` and `metadata_fungi_fix.xlsx`, which are
a manual edited version of `metadata_bacteria.tsv` and `metadata_fungi.tsv` created
with the previous step. File were then exported as TSV files.

## Create the configuration files for nextflow

Simply call:

```bash
nf-core launch nf-core/ampliseq -r 2.11.0
```

and follow the instructions to create the configuration files. Two configuration files
has been created: `config/nf-params-bacteria.json` and `config/nf-params-fungi.json`
for *bacteria* and *fungi* respectively.

## Determine the best trimming parameters

The dada creator suggested to use [figaro](https://github.com/Zymo-Research/figaro)
to determine which error rates and trimming parameters are better. I've found
different problems in running figaro since all the sequences are supposed to be 
of the same length: this is not true and using seqtk to remove all the sequence
and keep only the ones with the same length is a good idea drops out a lot of data.
The suggested parameters don't work with this data, so I've decided to manually 
explore the dada parameters. The different parameters combination are described 
in `analysis/01-about-dada2-filtering.ipynb` notebook. The best combination
was selected and used in the `config/nf-params-bacteria.json` file.

## Launch the Nextflow pipeline

To launch the pipeline, simply call:

```bash
nextflow run nf-core/ampliseq -r 2.11.0 -params-file config/nf-params-bacteria.json \
    -profile singularity -resume -c config/custom.config
```

for bacteria and

```bash
nextflow run nf-core/ampliseq -r 2.11.0 -params-file config/nf-params-fungi.json \
    -profile singularity -resume -c config/custom.config
```

for fungi.

> note: for fungi the pipeline need to be called. See
> [DADA2 ITS Pipeline Workflow (1.8)](https://benjjneb.github.io/dada2/ITS_workflow.html)
> for more information

## Analyze the results

The results of the pipeline are stored in the `results-bacteria` and `results-fungi`
folders. The `analysis` folder contains some notebooks to analyze the results.
This project can be managed using [R 4.4.2](https://cran.r-project.org/), 
[RStudio](https://posit.co/downloads/), [renv](https://rstudio.github.io/renv/articles/renv.html) 
and [targets](https://books.ropensci.org/targets/). 
You require also [quarto](https://quarto.org/docs/download/)
to render the final reports (should be present in a RStudio installation).

### Using Rstudio

Opening the `Agritech8_Metagenomics.Rproj` with Rstudio should initialize the
`renv` environment. When ready, install the required packages with:

```r
renv::restore()
```

Next, since targets is managed by projects, you need to declare which project
you want to compile, for example

```r
Sys.setenv(TAR_PROJECT = "technical_replicates")
```

`technical_replicates` is the name of the project you want to compile (inspect
`_targets.yaml` for the list of projects already configured).
You can check for pipeline errors using `targets` (a dependency installed with
`renv`):

```r
targets::tar_manifest(fields = all_of("command"))
```

You can also see the dependency graph of the pipeline, with information on which
process is outdated and need to be called by `targets`:

```r
targets::tar_visnetwork()
```

To call the pipeline, simply call:

```r
targets::tar_make()
```

This will run the pipeline and generate the results. The final reports are managed
with quarto, you should find the *html* files in `analysis` folder. To load data
from the *targets* pipeline, you can 
[tar_load](https://docs.ropensci.org/targets/reference/tar_load.html) to read a
specific object or [tar_load_everything](https://docs.ropensci.org/targets/reference/tar_load_everything.html)
to load all data in your *global environment*.

### Using GNU Make

If you don't want to use Rstudio, you can use GNU Make to compile the pipeline.
First, you need to install the required packages with:

```bash
make restore
```

Next you can build the projects one-by-one: first get a list of the available
projects with:

```bash
make list
```

Then simply build a project using its name, for example:

```bash
make technical_replicates
```

### Creating Krona report

In order to create report with *Krona*, you can call the `scripts/run_krona.R`
(or execute the proper rule in the Makefile). The script will generate the
krona input files and the command to run to generate the report.

## Add another analysis

Open `_targets.yaml` at the project root directory and define a new *target* section
with the source script and the store folder, for example:

```yaml
november_2023:
  script: scripts/07-november_2023.R
  store: stores/november_2023
```

Then create the script `scripts/07-november_2023.R` and add stuff using targets:
you can also copy another script and modify it. If you add new r packages, remember
to update the `renv.lock` file with:

```r
renv::install(<package_name>)
renv::snapshot()
```

Then you need to declare the target to compile with `TAR_PROJECT` environment
variable:

```r
Sys.setenv(TAR_PROJECT = "november_2023")
```

Remember also to update the `Makefile` in order to be able to compile the new project
without Rstudio, you need *only* to add the project name in the `PROJECTS` variable:

```make
# Start by specifying the project names and their corresponding scripts and stores.
PROJECTS := technical_replicates plot_iNEXT reactor_vs_algae duckweed_and_box november_2023
```

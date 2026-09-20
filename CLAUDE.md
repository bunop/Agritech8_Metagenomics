# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Amplicon metagenomics (16S bacteria / ITS fungi) analyses for Agritech Spoke 8. The work
is two-stage: an upstream `nf-core/ampliseq` Nextflow run produces phyloseq objects and
trees, and a downstream multi-project R `{targets}` pipeline turns those into statistics
and Quarto reports. See `README.md` for the provenance of the metadata files.

## Commands

R environment (renv; `renv.lock` pins R 4.5.2 / Bioconductor 3.21 — the README's
"R 4.4.2" is stale):

```bash
make restore        # renv::restore()
make list           # list the available targets projects
make <project>      # build one project, e.g. make november_2023-bacteria
make projects       # build them all
make all            # projects + krona_bacteria + krona_fungi
make prune          # tar_prune() every project
```

Equivalent from an R session (RStudio opens `Agritech8_Metagenomics.Rproj` and activates
renv automatically):

```r
Sys.setenv(TAR_PROJECT = "november_2023-bacteria")   # required: selects the project
targets::tar_manifest(fields = all_of("command"))    # inspect the plan
targets::tar_visnetwork()                            # dependency graph / what is outdated
targets::tar_make()
targets::tar_make(names = nov_2023_bray_distance_matrix)  # build a single target
```

There is no test suite, linter, or CI in this repo — `tar_make(names = ...)` on one
target is the way to exercise a single piece of the pipeline.

On the cluster, each project has a SLURM wrapper that activates the `R-4.5` conda env,
sets `TAR_PROJECT`, and runs `tar_make()` then `tar_prune()`:

```bash
sbatch scripts/07-november_2023-bacteria.sh
```

Upstream Nextflow (run outside this pipeline, produces `results-bacteria/` and
`results-fungi/`):

```bash
nf-core launch nf-core/ampliseq -r 2.11.0            # regenerate config/nf-params-*.json
nextflow run nf-core/ampliseq -r 2.11.0 -params-file config/nf-params-bacteria.json \
    -profile singularity -resume -c config/custom.config
```

Python (Poetry, `pyproject.toml`, `package-mode = false`) is used only for the metadata
notebooks `scripts/metadata.ipynb` and `scripts/samplesheet.ipynb`. `nbstripout` is a
dependency — notebooks are committed without outputs.

## Architecture

### Multi-project targets setup

There is **no `_targets.R`**. `_targets.yaml` declares six independent projects, each
with its own plan script and its own store:

```yaml
november_2023-bacteria:
  script: scripts/07-november_2023-bacteria.R
  store: stores/november_2023-bacteria
```

The `TAR_PROJECT` environment variable selects the active project; nothing works without
it. `Makefile`'s `PROJECTS` variable must be kept in sync with `_targets.yaml` — it is
the only place the make rules learn project names.

### Shared helpers

Every plan script calls `tar_source()` with no arguments, which sources all of `R/` into
every pipeline. Helpers are therefore global: changing one affects all six projects.

- `R/load_stuff.R` — `load_metadata()` (also derives `technical_rep`, `sample_name`,
  `date_condition` and the `label` factor), `load_phyloseq()`, `add_tree_to_phyloseq()`,
  `prune_phyloseq()`, `sort_phyloseq()` (microViz `ps_reorder`),
  `sort_by_factor_column()`, `get_otu_table()`, `get_samples_data()`,
  `subset_phyloseq()`
- `R/alpha.R` — `rarefy_alpha()` (phyloseq `estimate_richness` + Faith's PD via
  picante), `summarize_rarefactions()`, `reshape_rarefaction_data()`,
  `calculate_rarecurve()`, `calculate_kruskal_wallis()`, `calculate_dunn_test()`
- `R/beta.R` — `calculate_distance_matrix()` (vegan `vegdist`),
  `calculate_permanova()` (`adonis2`), `calculate_pairwise_permanova()`,
  `calculate_beta_dispersion()` (`betadisper` + ANOVA), `calculate_pcoa()` (ape,
  Cailliez), `calculate_nmds()` (`metaMDS`), `get_distances()`
- `R/merge.R` — technical replicates: `calculate_technical_replicate_distances()`,
  `calculate_technical_replicate_correlations()`, `merge_metadata()`,
  `merge_technical_replicates()`
- `R/misc.R` — `rarefy_phyloseq_object()`, `agglomerate_by_taxa()`,
  `phyloseq_to_ampvis2()`, `calculate_coverage_stats()`
- `R/plots.R` — `custom_heatmap()`, `plot_alpha_diversity()`, `plot_distances()`
- `R/krona.R` (`get_krona_cmd()`), `R/iNEXT_helper.R` (`combine_iNEXT_results()`)

Check these before writing a new function; the plan scripts are almost entirely calls
into them.

### Data flow of a typical plan

`data/metadata_{bacteria,fungi}_fix.tsv` + `results-*/phyloseq/dada2_phyloseq.rds` +
`results-*/qiime2/phylogenetic_tree/tree.nwk` (all tracked as `format = "file"` targets)
→ phyloseq object with tree → prune (bacteria only: drop Cyanobacteria/chloroplast) →
subset and order samples → rarefy 1000× → alpha/beta diversity, PERMANOVA, betadisper,
ordinations, heatmaps → `tar_quarto()` report.

`results-*/`, `data/`, `*.rds`, `*.png`, `*.pdf`, `*.html` are all gitignored.

### Reports

Each plan ends with `tarchetypes::tar_quarto()` rendering `analysis/0X-<project>.qmd`,
so `tar_make()` regenerates the HTML. The qmd does **not** depend on `TAR_PROJECT`; it
reads results with an explicit store path:

```r
store_path <- here::here("stores", "november_2023-bacteria")
targets::tar_load(nov_2023_metadata, store = store_path)
```

Supplementary tables are written by the report itself (`openxlsx::write.xlsx`) into
`analysis/0X-<project>/`, which `analysis/.gitignore` excludes. `analysis/.Rprofile`
activates the parent renv so a Jupyter R kernel started in `analysis/` finds the library.

### Cloud storage

`scripts/07-november_2023-bacteria.R` is the **only** project using S3 so far (still an
"attempt" — the other five are local-store):

```r
library(paws.storage)
tar_option_set(
  repository = "aws",
  repository_meta = "aws",
  resources = tar_resources(aws = tar_resources_aws(
    bucket = "r-targets-hay4l",
    prefix = "agritech8-metagenomics/november_2023-bacteria"
  ))
)
```

Credentials come from `AWS_PROFILE = "personal-s3"`, set in `.Rprofile` — i.e. a named
profile in `~/.aws/`, no keys or region in the repo. Neither the `Makefile` nor the
SLURM wrappers propagate AWS settings, so R must start in the project root for
`.Rprofile` to apply.

### Parallelism

`crew::crew_controller_local(workers = parallelly::availableCores() - 1,
seconds_idle = 60, tasks_max = 50)` with `seed = 42`. The 1000 rarefaction branches are
the reason crew is there.

## Conventions

- **Target names** are project-prefixed snake_case in scripts 07/08 (`nov_2023_*`) and
  unprefixed in 02/05/06. The suffix grammar is regular and worth following:
  `*_distance_matrices` (branched) → `*_distance_matrix` (averaged) → `*_pcoa_object` /
  `*_nmds_object` / `*_beta_dispersion` / `*_permanova` / `*_pairwise_*_permanova`, with
  grouping suffixes `_by_sample_name`, `_by_labels`. A `_tmp` suffix marks the stage
  before a sorting/ordering step.
- **Rarefaction idiom**: a `thousand_iterations = seq_len(1000)` target drives
  `pattern = map(...)`, `iteration = "list"`; replicate distance matrices are averaged
  with `as.dist(Reduce("+", x) / length(x))`. Adding a metric means repeating this
  matrices → matrix pair. UniFrac/wUniFrac call `phyloseq::distance()` directly in the
  target rather than `calculate_distance_matrix()`.
- **Sample labels**: `T0` = default, `NA` = reactor, `NP` = box, `CS` = reactor+cs,
  `LM` = box+dw. `custom_order_labels <- c("T0", "NA", "NP", "CS", "LM")` is itself a
  target, threaded through sorting, heatmaps, betadisper levels and alpha plots so
  everything stays aligned. Recent work moved comparisons from `date_condition` to
  `label`.
- **Technical replicates** share a `sample_name` (sampleID minus its trailing letter).
  In project 07 the merged object (`nov_2023_phyloseq_obj_merged`) feeds **only the
  barplots**; heatmaps and every diversity/ordination statistic still use the unmerged
  object. Merging is not yet ported to `scripts/08-november_2023-fungi.R`.
- **Fungi (ITS) differs**: no chloroplast pruning, and metadata is filtered *after*
  subsetting to keep only sampleIDs present in the phyloseq object, because some samples
  failed sequencing.
- **Stores**: `stores/*/.gitignore` commits only `meta/meta`, never `objects/`. Never
  hand-edit it — and expect any `tar_make()` to produce a large, noisy diff on that file
  (thousands of lines). That diff is normal, not a mistake.
- **Adding a project**: add the section to `_targets.yaml`, create
  `scripts/0X-<name>.R` (copy the closest existing one), append the name to `PROJECTS`
  in `Makefile`, add `analysis/0X-<name>.qmd`, and optionally a `scripts/0X-<name>.sh`
  SLURM wrapper.
- **New R packages**: `renv::install(<pkg>)` then `renv::snapshot()`, and add the
  package to the `packages` vector in `tar_option_set()` of the plan that needs it —
  targets workers do not inherit the session's attached packages.

## Domain stack

phyloseq (the core object), vegan, ape, picante, microViz, ampvis2, FSA,
pairwiseAdonis, iNEXT, plus tidyverse. Note there is **no differential abundance
analysis** — no DESeq2, ANCOM or ALDEx2. Krona reports are generated out-of-band by
`scripts/run_krona.R` (`make krona_bacteria` / `make krona_fungi`), which writes the
input files and prints the `ktImportText` command to run.

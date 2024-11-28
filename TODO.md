---
title: "TODO"
author: "Paolo Cozzi"
editor: source
---

- Merge samples with [phyloseq](https://joey711.github.io/phyloseq/merge.html)
- Explore ampvis2 package
  - see [heatmap](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html#heatmap)
- Explore [metagMisc](https://github.com/vmikk/metagMisc)
  - perform multiple rarefaction like Qiime2 does:
  ```r
  # Assuming 'physeq' is your phyloseq object
  rarefied_diversity <- phyloseq_mult_raref_div(
    physeq = physeq,
    SampSize = min(sample_sums(physeq)), # Rarefaction depth
    iter = 1000,                         # Number of iterations
    divindex = c("Observed", "Shannon"), # Diversity indices
    parallel = FALSE,                    # Set to TRUE for parallel processing
    verbose = TRUE
  )
  ```

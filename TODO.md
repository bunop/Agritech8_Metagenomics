---
title: "TODO"
author: "Paolo Cozzi"
editor: source
---

- Merge samples with [phyloseq](https://joey711.github.io/phyloseq/merge.html)
  - create a new quarto document for merged samples
  - read merged data from `_target` folder
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
- Explore [qiime2R](https://github.com/jbisanz/qiime2R)
  - see [tutorial](https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121)
- See [Alpha Diversity Tutorial](https://rstudio-pubs-static.s3.amazonaws.com/1071936_6115f873acbc4dc4a30b1380cc3885fb.html)
  - Run the Kruskall-Wallis test

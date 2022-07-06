gscramble-paper-sims-and-analyses
================

This repo holds all the code for reproducing the simulations, analyses,
and figures from the paper, “Title still not determined” by Authors,
Authors, Authors, about the R package ‘gscramble’.

It is set up to be relatively reproducible. Running it all the way we
did requires Snakemake.

The R packages for this are maintained by the ‘renv’ package. So, run it
with R version 4.2.0 (2022-04-22), and do these steps:

1.  Open the RStudio project.
2.  Issue `renv::restore()` in the console.
3.  Activate your snakemake environment in your Unix terminal
4.  Do `snakemake --use-conda --cores 8` to use 8 cores. Set the number
    of cores as appropriate to your system.

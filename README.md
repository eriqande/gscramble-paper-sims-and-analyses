gscramble-paper-sims-and-analyses
================

[![DOI](https://zenodo.org/badge/504563380.svg)](https://doi.org/10.5281/zenodo.14270660)

This repo holds the code and data for reproducing the simulations,
analyses, and figures performed by Eric C. Anderson that accompany the
paper, “GSCRAMBLE: Simulation of admixed individuals without reuse of
genetic material” by Eric C. Anderson, Rachael M. Giglio, Matthew G.
DeSaix and Timothy J. Smyser, to be published in Molecular Ecology
Resources.

The R package described in the paper is [available on
CRAN](https://CRAN.R-project.org/package=gscramble)

This repository is set up to be easily reproducible. Running all of it
the way we did requires Snakemake.

The R packages for this are maintained by the ‘renv’ package. So, run it
with R version 4.2.0 (2022-04-22), and do these steps:

1.  Open the RStudio project.
2.  Issue `renv::restore()` in the console.
3.  Activate your snakemake environment in your Unix terminal
4.  Do `snakemake --use-conda --cores 8` to use 8 cores. Set the number
    of cores as appropriate to your system.

## Notes during development:

To run this on Sedna, I did this:

``` sh
module load R/4.0.3
R

# and when it started it bootstrapped renv together and then
# in R I did:
renv::restore()

# that took a while, because it built a lot of it from source
# on Linux, but eventually got done.

# then I gave a test run (2 reps for each condition) with:
snakemake --cores 20 --use-conda --use-envmodules results/compiled/Q-values-from-sims.tsv.gz
```

To do the actual simulations on SEDNA I did this:

``` sj
snakemake --profile hpcc-profiles/slurm/sedna results/compiled/Q-values-from-sims.tsv.gz
```

### Some notes on GSPs to use

To get large numbers of simulated admixed individuals of different
admixture fractions when using gscramble we note that:

- `create_GSP(pop1 = "A", pop2 = "B", T, T, T, T)` gives us:
  - 2 Bx2-A’s as s10
  - 1 BX1-A as s9
  - 2 F2s as s11
  - 1 F2 as s7

And that only consumes 4 As and 2 Bs. So, that covers most of our needs,
except for a BX1 x F1 or BX1 x F2, which will give us the Q = 0.375
category.

### While developing the sims…

I have a couple little test cases on the original sims to run on my
laptop like this:

``` sh
snakemake -p  'results/bias_sims/freq_A/Qs_1/L_10000/N_50/n_3/Rep_113/unsupervised_Q.tsv' 'results/bias_sims/freq_A/Qs_1/L_10000/N_50/n_3/Rep_114/supervised_Q.tsv' --cores 8 --use-conda

or 

snakemake -p  'results/bias_sims/freq_A/Qs_1/L_100/N_50/n_3/Rep_113/unsupervised_Q.tsv' 'results/bias_sims/freq_A/Qs_1/L_100/N_50/n_3/Rep_114/supervised_Q.tsv' --cores 8 --use-conda
```

And I can use those to test versions using gscramble.

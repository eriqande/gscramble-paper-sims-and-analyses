log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

source("renv/activate.R");

# This script makes a set of plink files to feed into admixture
library(tidyverse)
source("R/bias-sim-funcs.R")

# Example Parameter values.  Will ultimately get these from snakemake.
if(snakemake@params$freq == "A") {
  fp <- list(L = 10000, a = 1, b = 8, t = 0.01)
} else {
  stop("Not a known frequency spec!")
}


if(snakemake@params$Q == "1") {
  Qs <- c(0, 0.125, 0.25, 0.375, 0.5)
} else {
  stop("Unknown Qs spec!")
}

N <- as.integer(snakemake@params$N)
n <- as.integer(snakemake@params$n)

fp$L <- as.integer(snakemake@params$L)

prefix <- snakemake@params$prefix

# The seed will be a transformed hash of the prefix
Rseed <- rlang::hash(prefix) %>%
  str_replace_all("[^0-9]", "") %>%
  str_sub(end = 8) %>%
  as.integer()

cat(Rseed, file = file.path(dirname(prefix), "R.seed"))
set.seed(Rseed)

# simulate the reference individuals
G <- sim_inds(fp, N)

# simulate the admixed individuals and add them in there
A <- sim_admixed(G, N, Qs, n)

# write it out
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
write_to_plink(M = A, prefix = prefix)

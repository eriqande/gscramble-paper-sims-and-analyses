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


if (snakemake@params$Q == "1" || snakemake@params$Q == "1_sim_the_knowns") {
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

# here we do some special things for the case when
# snakemake@params$Q is "1_sim_the_knowns". Namely, we
# replace the N reference individuals from each population
# with N individuals simulated by drawing with replacement
# from the original N.  This should produce even more profound
# biases.
if (snakemake@params$Q == "1_sim_the_knowns") {
  # simulate the new ones
  G2 <- sim_admixed(G = G, N = N, Qs = 0, n = N)
  # pick out the new ones
  G2_new <- G2[-(1:(2 * N)),]
  # reorder these new ones so that pop1 indivs are first:
  G2_n2 <- G2_new[c((N+1):(2*N), 1:N), ]
  # then rename them to be s1 and s2
  rownames(G2_n2)[1:N] <- "s1"
  rownames(G2_n2)[(N + 1):(2 * N)] <- "s2"
  
  # then replace the top 2N individuals in A with the G2_n2 guys
  A2 <- rbind(
    G2_n2,
    A[-(1:(2 * N)), ]
  )
  
  # set that to A and move ahead
  A <- A2
}


# write it out
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
write_to_plink(M = A, prefix = prefix)

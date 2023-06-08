log <- file(snakemake@log$log, open="wt")
sink(log, type = "output")
sink(log, type = "message")
saveRDS(snakemake, file = snakemake@log$snake_obj)


source("renv/activate.R");

# This script makes a set of plink files to feed into admixture
library(tidyverse)
source("R/bias-sim-funcs.R")

#snakemake <- read_rds("results/bias_sims/freq_A/Qs_1/L_120/N_50/n_3/Rep_113/snake_obj.rds")

sim_type <- snakemake@params$sim_type
if (!(sim_type %in% c("bias_sims", "gscramble_unlinked", "gscramble_linked"))) {
  stop("The bias_sims wildcard must be either bias_sims or gscramble. ")
}

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


if (sim_type == "bias_sims") {
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
} else {
  # here we are doing it via gscramble
  
  # Make a two-column per locus format.  Let the alleles be A for 0 and T for 1
  first_gc <- G
  first_gc[G == 1 | G == 2] <- "T"
  first_gc[G == 0] <- "A"
  second_gc <- G
  second_gc[G == 0 | G == 1] <- "A"
  second_gc[G == 2] <- "T"
  
  Gc <- matrix(rbind(first_gc, second_gc), nrow = nrow(G))
  
  Ge <- Gc   # Ge is our "Genos"
  
  # now, make the Imeta.  Assume the first half are from one pop and the second
  # are from another.  Im is our I_meta
  Im <-  tibble(
    group = c(rep("A", nrow(G) / 2), rep("B", nrow(G) / 2)),
    indiv = paste0(c(rep("A", nrow(G) / 2), rep("B", nrow(G) / 2)), "_", c(1:(nrow(G) / 2), 1:(nrow(G) / 2)))
  )
  
  
  # Now, let us also make some marker meta data.  We will assume a pig-sized genome from
  # RecRates, and we will put just randomly broadcast the markers in order onto those
  # chromosomes.   RR is our RecRates.
  RR <- read_rds("input/RecRates_Big.rds")
  
  # get the length of each chromosome and then sample the number of markers
  # on each in a proportional fashion
  clengths <- RR %>%
    mutate(cf = factor(chrom, levels = unique(chrom)), .before = chrom) %>%
    group_by(chrom) %>%
    summarise(length = chrom_len[1], max_pos = max(end_pos))
  
  num_on_chr <- rmultinom(1, size = ncol(G), prob = clengths$length)[,1]
  
  # now, our markers are on different chromosomes and we broadcast them
  # uniformly within each. 
  CL_plus <- clengths[ rep(1:nrow(clengths), num_on_chr), ] %>%
    group_by(chrom) %>%
    mutate(pos = sort(floor(runif(n = n(), min = 2, max = length)))) %>%
    ungroup()
  
  # Finally we name the loci locus_1, up to  Locus_L where L is the number
  # of them
  Mm <- bind_cols(CL_plus, tibble(variant_id = paste0("Loc_", 1:ncol(G)))) %>%
    select(chrom, pos, variant_id)
  
  # Now, we have everything that we would need to run gscramble: Ge, Im, Mm, and RR.
  # So, we will do that in a function.
  
}

# write it out
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
write_to_plink(M = A, prefix = prefix)

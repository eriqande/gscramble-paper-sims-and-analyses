require(tidyverse)
require(gscramble)

#' This simulates the allele freqs from a and b and L
#' @param L the number of SNPs
#' @param a the first beta param
#' @param b the second beta param
#' @param t MAF frequency cutoff
#' @examples
#' p <- sim_freqs(1000, 1, 10, 0.05)
sim_freqs <- function(L, a, b, t) {
  p <- rbeta(L, a, b)
  p <- p[p>t & p<(1-t)]
  while(length(p) <= L) {
    pn <- rbeta(L, a, b)
    pn <- pn[pn>t & pn<(1-t)]
    p <- c(p, pn)
  }
  p[1:L]
}

#' This simulates the 2 samples of size N from the allele freqs
#' 
#' Some care is taken to make sure that the the samples are
#' not monomorphic for any site.  If need be, it keeps simulating
#' new allele freqs to add more sites that are monomorphic.
#' @param fp freq parameters. A list with L, a, b, and t.
#' @param N the number of individuals in each sample.  The first
#' N will be from the first sample and the second N from the second
#' sample.
#' @examples
#' fp <- list(L = 1000, a = 1, b = 10, t = 0.05)
#' inds <- sim_inds(fp, 25)
sim_inds <- function(fp, N) {
  ret <- matrix()
  while(ncol(ret) < fp$L) {
    p <- do.call(sim_freqs, fp)
    # turn those into genotype freqs
    g <- rbind(
      (1 - p) ^ 2,
      2 * p * (1 - p),
      p ^ 2
    )
    # sample 2N genotypes (0, 1, 2) from those
    mat <- apply(g, 2, function(x) sample(0:2, 2 * N, replace = TRUE, prob = x))
    
    # add those to the existing polymorphic ones
    if(nrow(ret) == 1) {
      ret <- mat
    } else {
      ret <- cbind(ret, mat)
    }
    
    # drop loci that are monomorphic
    ret <- ret[, colSums(ret) > 0 & colSums(ret) < 4 * N]
  }
  
  ret[, 1:fp$L]
}


#' Given a genotype matrix like that from sim_inds, this simulates the admixed individuals
#'
#' This will create n individuals with q = Qs and n individuals with q = 1 - Qs
#' for each value in the Qs vector,
#' and it will put them on the bottom of the G matrix. NOTE: if the Q-value is
#' 0.5 it will only do n of them, not 2n.
#' @param G the genotype matrix (0, 1, 2)
#' @param N the number from each baseline sample.  This is here mostly
#' to check that the matrix G has the right number of rows.
#' @param Qs the vector of Q-values desired.  They should be <= 0.5.
#' These are the fraction of gene copies from the first "population".
#' @param n the number of individuals
#' @examples
#' fp <- list(L = 1000, a = 1, b = 10, t = 0.05)
#' inds <- sim_inds(fp, 25)
#' Qs <- c(0, 0.125, 0.25, 0.375, 0.5)
#' smat <- sim_admixed(inds, 25, Qs, 3)
sim_admixed <- function(G, N, Qs, n) {
  # get the Q values we will simulate each n at:
  Q2 <- sort(unique(c(Qs, 1 - Qs)))
  names(Q2) <- Q2
  
  # get the allele freqs in the two samples
  stopifnot(nrow(G) == 2 * N)
  f1 <- colMeans(G[1:N,] / 2)
  f2 <- colMeans(G[(N+1):(2*N),] / 2)
  L <- length(f1)
  # make n copies of each q value and then cycle over those
  Q3 <- rep(Q2, each = n)
  glist <- lapply(Q3, function(q) {
    # get the allele frequency of the first gene copy (depends on which
    # population it came from)
    f <- ifelse(runif(L) < q, f1, f2)
    # then sample the allelic types.  It is a 1 if an runif is less than the freq
    g1 <- as.integer(runif(L) < f)
    # then do the same for the second gene copy
    f <- ifelse(runif(L) < q, f1, f2)
    g2 <- as.integer(runif(L) < f)
    # make the genotype
    g1 + g2
  })
  
  ad_mat <- do.call(rbind, glist)
  rownames(ad_mat) <- names(Q3)
  
  # finally, put them on the end with the rownames there
  rownames(G) <- c(rep("s1", N), rep("s2", N))
  
  rbind(G, ad_mat)
}


# Cool! Now we just need a function to write the output of
# sim_admixed() to PLINK format.  We will use gscramble2plink()
# for that. 


#' Write the output of sim_admixed to plink
#' @param M the output of sim_admixed.
#' @param prefix the PLINK file prefix to write it to
write_to_plink <- function(M, prefix = "plink") {
  Im <- tibble(
    indiv = rownames(M)
  ) %>%
    mutate(
      group = case_when(
        indiv == "s1" ~ "reference_1",
        indiv == "s2" ~ "reference_2",
        TRUE ~ str_c("q-", indiv)
      )
    ) %>%
    mutate(indiv = str_c("M_", 1:n(),"_", indiv)) %>%
    select(group, indiv)
  
  # Now, it is slightly more involved to get these into a
  # "two-columns-per-locus" matrix
  firsts <- M
  firsts[M == 0] <- 1
  firsts[M == 1] <- 1
  firsts[M == 2] <- 2
  seconds <- M
  seconds[M == 0] <- 1
  seconds[M == 1] <- 2
  seconds[M == 2] <- 2
  tc_mat <- matrix(rbind(firsts, seconds), nrow = nrow(M))
  
  #### this was for testing that it was done correctly ####
  # tibble(
  #   row = sample(1:nrow(M), 1000, replace = TRUE),
  #   col = sample(1:ncol(M), 1000, replace = TRUE)
  # ) %>%
  #   mutate(
  #     one = map2_dbl(row, col, function(x, y) M[x, y]),
  #     two = map2(row, col, function(x, y) tc_mat[x, c(y*2-1, y*2)]),
  #     str = map_chr(two, function(x) sprintf("%s-%s", x[1], x[2]))
  #   ) %>%
  #   count(one, str)
  
  # finally, some (made up / bogus) marker info
  M_meta <- tibble(
    chrom = "1",
    pos = (1:ncol(M)),
    variant_id = str_c("Locus_", 1:ncol(M))
  )
  
  # write the plink map and ped files
  gscramble2plink(Im, M_meta, tc_mat, prefix)
  
  # write the pop file:
  Im %>%
    mutate(
      pop = case_when(
        group == "reference_1" ~ "s1",
        group == "reference_2" ~ "s2",
        TRUE ~ "-"
      )
    ) %>%
    pull(pop) %>%
    cat(., sep = "\n", file = paste0(prefix, ".pop"))
  
  # and, finally, write the key
  write_tsv(Im, file = file.path(dirname(prefix), "key.tsv"), col_names = FALSE)
}

  
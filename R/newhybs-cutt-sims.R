# this is the script to be called by Snakemake for each separate simulation


# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


library(tidyverse)
library(gscramble)



if(!exists("snakemake")) {
  input_list <- list(
    rda = "results/cutthroat-sim-prep/objects.rda",
    nhcats = "input/newhybs_categories.txt"
  )
  output_list = list(
    nhdat = "results/newhyb_sims/markers_65-rep_1/nh_dat.txt",
    nhout = "results/newhyb_sims/markers_65-rep_1/aa-PofZ.txt"
  )
  params_list = list(
    nmark = "65"
  )
} else {
  input_list = snakemake@input
  output_list = snakemake@output
  params_list = snakemake@params
}


load(input_list$rda)


# get the number of markers
Nmark <- as.integer(params_list$nmark)

# set seed from a hash
seed_base <- output_list$nhout  # use a seed that depends on the rep
seed <- rlang::hash(seed_base) %>%
  str_replace_all("[^0-9]", "") %>%
  str_sub(end = 7) %>%
  as.integer()

set.seed(seed)



# get the directory and create it
dd <- dirname(output_list$nhout)
dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# make the segments
Segments2 <- segregate(
  request = InputTibble,
  RR = cutt_RecRates,
  MM = cutt_M_meta
)


# get the markers
Markers <- segments2markers(
  Segs = Segments2,
  Im = cutt_I_meta,
  Mm = cutt_M_meta,
  G = cutt_Geno
)


# rename indivs according to their hybrid category
Markers2 <- rename_indivs(Markers, hyb_key)



# make a newhybrids file., chooseing 65 or 30 loci as appropriate
if(Nmark == 65) {
  gs2nh_ret <- gscramble2newhybrids(
    Markers2, 
    cutt_M_meta, 
    z = c("CCT", "SH"),
    s = "CCT|SH",
    outfile = output_list$nhdat
  )
} else if(Nmark == 30) {
  gs2nh_ret <- gscramble2newhybrids(
    Markers2, 
    cutt_M_meta, 
    z = c("CCT", "SH"),
    s = "CCT|SH",
    outfile = output_list$nhdat,
    retain = markers30
  )
} else {
  stop("Nmark must be 30 or 65")
}

# get an absolute path for the newhybs categories
nhcats_path <- normalizePath(input_list$nhcats)

# get the basename of the input file
nhdat <- basename(output_list$nhdat)

# get two seeds for newhybrids
seed1 = floor(runif(1) * 1e7)
seed2 = floor(runif(1) * 1e7)

# make a call for system()
CALL <- glue::glue("
cd {dd};
../../../bin/newhybsng-$(uname) -d {nhdat} -c {nhcats_path} --no-gui -s {seed1} {seed2} --num-sweeps 25000 --burn-in 5000
")

system(CALL)

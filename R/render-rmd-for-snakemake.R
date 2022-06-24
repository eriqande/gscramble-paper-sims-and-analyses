

# this is a simple R script let's let me sanely
# render an Rmd with Snakemake, and put the html output where I
# want to, and also to have other output files that are specified
# in the snakefile rule block.  Also, chunks get evaluated in the
# current working directory.  

# Assumes:
#   - the Rmd input file is named "Rmd" in the snakemake input object
#   - the output html file is named html in the output object
#   - that renv is being used for R package management


# this loads the package library managed by renv
source("renv/activate.R");

# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# snakemake doesn't seem to make the necessary directories
# for output when running R?  Not sure what is up there...
# dump <- lapply(snakemake@output, function(f) {
#   dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
#   print(dirname(f))
# })

outfile <- basename(snakemake@output$html)
outdir <- dirname(snakemake@output$html)

rmarkdown::render(
  input = snakemake@input$Rmd,
  output_file = outfile,
  output_dir = outdir,
  knit_root_dir = getwd()
)



library(gscramble)


#### Make a simple F2 GSP  ####
F2 <- create_GSP("A", "B", F2 = TRUE)

dir.create("results/figures/gsps-for-paper", recursive = TRUE, showWarnings = FALSE)
gsp2dot(
  F2, 
  path = "results/figures/gsps-for-paper/F2", 
  indiv_node_label_font_size = 36, 
  sample_node_label_font_size = 36, 
  edge_label_font_size = 36
)

system(" rsvg-convert -f pdf -o results/figures/gsps-for-paper/F2.pdf results/figures/gsps-for-paper/F2.svg")

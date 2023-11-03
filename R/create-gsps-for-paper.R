

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




#### Make a GSP with and F1, F2, BC1 and BC_2  ####
AllThem <- create_GSP("A", "B", F1 = TRUE, F2 = TRUE, F1B = TRUE, F1B2 = TRUE)

dir.create("results/figures/gsps-for-paper", recursive = TRUE, showWarnings = FALSE)
gsp2dot(
  AllThem, 
  path = "results/figures/gsps-for-paper/AllThem", 
  indiv_node_label_font_size = 36, 
  sample_node_label_font_size = 36, 
  edge_label_font_size = 36
)

system(" rsvg-convert -f pdf -o results/figures/gsps-for-paper/AllThem.pdf results/figures/gsps-for-paper/AllThem.svg")

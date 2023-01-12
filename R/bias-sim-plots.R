# script to make plots from the bias simulation results


library(tidyverse)

file <- "results/compiled/Q-values-from-sims.tsv.gz"


# read in and extract parts from filenames, etc.
full_dat <- read_tsv(
  file, 
  col_names = c("path", "group", "sample", "qstr")
) %>%
  extract(
    path,
    into = c("freq_spec", "Q_spec", "L", "N", "n", "Rep", "scond"),
    regex = ".*freq_([A-Z])/Qs_([0-9])/L_([0-9]+)/N_([0-9]+)/n_([0-9]+)/Rep_([0-9]+)/(unsupervised|supervised)_Q\\.tsv",
    convert = TRUE
  ) %>%
  separate(
    qstr,
    into = c("q1", "q2"),
    sep = " ",
    convert = TRUE
  )

  
# first we will plot the supervised ones.  We can toss
# the reference individuals' Q-values in that case
supes <- full_dat %>%
  filter(scond == "supervised" & !str_detect(group, "reference")) %>%
  mutate(
    simQ = as.numeric(str_replace(group, "^q-", ""))
  ) %>%
  mutate(
    simQ = factor(simQ, levels = sort(unique(simQ)))
  )


sup_plot <- ggplot(supes, aes(x = simQ, y = q1, colour = factor(n))) +
  geom_abline(intercept = -1/8, slope = 1/8, linetype = "dotted") +
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.2) +
  facet_grid(N ~ L) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
  ) +
  scale_colour_manual(values = c(`3` = "#1b9e77", `12` = "#d95f02", `24` = "#7570b3")) +
  guides(
    colour = guide_legend(title = "n")
  ) +
  xlab("Simulated Q value") +
  ylab("Estimated Q value")


ggsave(sup_plot, filename = "results/figures/supervised_Qs.pdf", width = 10, height = 6)

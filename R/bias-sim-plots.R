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
    into = c("scenario", "freq_spec", "Q_spec", "L", "N", "n", "Rep", "scond"),
    regex = "results/(.*)/freq_([A-Z])/Qs_(1|1_sim_the_knowns)/L_([0-9]+)/N_([0-9]+)/n_([0-9]+)/Rep_([0-9]+)/(unsupervised|supervised)_Q\\.tsv",
    convert = TRUE
  ) %>%
  separate(
    qstr,
    into = c("q1", "q2"),
    sep = " ",
    convert = TRUE
  )

# we are going to nest this into different groups
nested_full <- full_dat %>%
  group_by(scenario, Q_spec, scond) %>%
  nest()
  

## Here are some functions that will get re-used a lot ##
MakePlot <- function(D, supervised) {
  if(supervised == TRUE) {
    # for the supervised ones we remove the reference individuals
    input <- D %>%
    filter(!str_detect(group, "reference")) %>%
      mutate(
        simQ = as.numeric(str_replace(group, "^q-", ""))
      ) %>%
      mutate(
        simQ = factor(simQ, levels = sort(unique(simQ)))
      )
    
    intercept = -1 / 8
  } else {
    # the unsupervised ones need some label switching taken care of.
    # we do some quick-and-dirty labeling here---since it is not supervised,
    # we don't know which cluster is ref-1 and which is ref-2.  So, we use a simple
    # rule to apply things in each iteration. Basically, we want ref-2 to be high in q1.
    input0 <- D %>%
      mutate(
        group2 = str_replace(group, "reference_", "ref-"),
        group2 = case_when(
          group2 == "ref-1" ~ "ref-2",
          group2 == "ref-2" ~ "ref-1",
          TRUE ~ group2
        )
      ) %>%
      mutate(
        simQ = factor(group2, levels = c("ref-1", sort(unique(group2[!str_detect(group2, "ref")])), "ref-2"))
      )
    
    # if the highest mean Q1 value is to reference-2 individuals,
    # then q-vals don't need flipping.  And it also turns out that ref-1 is akin
    # to q-0 here, so we need to flip the ref labels, too.
    flippers <- input0 %>%
      filter(str_detect(group2, "ref")) %>%
      group_by(L, N, n, Rep, group2) %>%
      summarise(meanQ1 = mean(q1)) %>%
      filter(meanQ1 == max(meanQ1)) %>%
      mutate(flip_it = (group2 == "ref-1")) %>%
      ungroup() %>%
      select(-group2, -meanQ1)
    
    input <- input0 %>%
      left_join(flippers, by = c("L", "N", "n", "Rep")) %>%
      mutate(
        nq1 = ifelse(flip_it, q2, q1),
        nq2 = ifelse(flip_it, q1, q2),
        q1 = nq1,
        q2 = nq2
      ) %>%
      mutate(
        simQ = fct_recode(simQ, `ref-1` = "ref-2", `ref-2` = "ref-1" )
      )
    
    intercept = -2 / 8
    
  }
  
  # then, make the plot
  ggplot(input, aes(x = simQ, y = q1, colour = factor(n))) +
    geom_abline(intercept = intercept, slope = 1/8, linetype = "dotted") +
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
  
  
}
  
  
# first we will plot the supervised ones.  We can toss
# the reference individuals' Q-values in that case
sup_plot <- nested_full %>%
  filter(scenario == "bias_sims", scond == "supervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = TRUE)


# Then we can do the unsupervised ones.
# Note that we do some quick-and-dirty labeling here---since it is not supervised,
# we don't know which cluster is ref-1 and which is ref-2.  So, we use a simple
# rule to apply things in each iteration. Basically, we want ref-2 to be high in q1.
unsup_plot <- nested_full %>%
  filter(scenario == "bias_sims", scond == "unsupervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = FALSE)


#### Now do the sim_the_knowns ####

# these are the cases in which all the samples were re-simulated, including
# the "reference" one in either the simulated or unsimulated cases.
stk_sup_plot <- nested_full %>%
  filter(scenario == "bias_sims", scond == "supervised", Q_spec == "1_sim_the_knowns") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = TRUE)




stk_unsup_plot <- nested_full %>%
  filter(scenario == "bias_sims", scond == "unsupervised", Q_spec == "1_sim_the_knowns") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = FALSE)




#### OK, now we want to get the results from running the same sorts of data sets through gscamble

## Gscramble supervised ##

gs_sup_plot <- nested_full %>%
  filter(scenario == "gscramble", scond == "supervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = TRUE) +
  geom_point(colour = "blue", alpha = 0.05)


## gscramble unsupervised

gs_unsup_plot <- nested_full %>%
  filter(scenario == "gscramble", scond == "unsupervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = FALSE) +
  geom_point(colour = "blue", alpha = 0.05)
  



## gscramble supervised with differentiated source pops

gs_diff_sup_plot <- nested_full %>%
  filter(scenario == "gscramble_diff", scond == "supervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = TRUE) +
  geom_point(colour = "blue", alpha = 0.05)


## gscramble unsupervised with differentiated source pops

gs_diff_unsup_plot <- nested_full %>%
  filter(scenario == "gscramble_diff", scond == "unsupervised", Q_spec == "1") %>%
  pull(data) %>%
  pluck(1) %>%
  MakePlot(supervised = FALSE) +
  geom_point(colour = "blue", alpha = 0.05)


### Write them to files ###
ggsave(sup_plot, filename = "results/figures/supervised_Qs.pdf", width = 10, height = 6)
ggsave(unsup_plot, filename = "results/figures/unsupervised_Qs.pdf", width = 10, height = 6)

ggsave(stk_sup_plot, filename = "results/figures/sim_the_knowns_supervised_Qs.pdf", width = 10, height = 6)
ggsave(stk_unsup_plot, filename = "results/figures/sim_the_knowns_unsupervised_Qs.pdf", width = 10, height = 6)


ggsave(gs_sup_plot, filename = "results/figures/gscram_supervised_Qs.pdf", width = 5, height = 3)
ggsave(gs_unsup_plot, filename = "results/figures/gscram_unsupervised_Qs.pdf", width = 5, height = 3)


ggsave(gs_diff_sup_plot, filename = "results/figures/gscram_diff_supervised_Qs.pdf", width = 5, height = 3)
ggsave(gs_diff_unsup_plot, filename = "results/figures/gscram_diff_unsupervised_Qs.pdf", width = 5, height = 3)



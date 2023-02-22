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
    regex = ".*freq_([A-Z])/Qs_(1|1_sim_the_knowns)/L_([0-9]+)/N_([0-9]+)/n_([0-9]+)/Rep_([0-9]+)/(unsupervised|supervised)_Q\\.tsv",
    convert = TRUE
  ) %>%
  separate(
    qstr,
    into = c("q1", "q2"),
    sep = " ",
    convert = TRUE
  )

# break it into the two different groups. 

# keep this for later!
sim_the_knowns <- full_dat %>%
  filter(Q_spec == "1_sim_the_knowns")

# call full_dat the
# one with Q_spec == 1, to fit with the first part of the code below.
full_dat <- full_dat %>%
  filter(Q_spec == "1")
  
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


# Then we can do the unsupervised ones.
# Note that we do some quick-and-dirty labeling here---since it is not supervised,
# we don't know which cluster is ref-1 and which is ref-2.  So, we use a simple
# rule to apply things in each iteration. Basically, we want ref-2 to be high in q1.
unsupes0 <- full_dat %>%
  filter(scond == "unsupervised") %>%
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
flippers <- unsupes0 %>%
  filter(str_detect(group2, "ref")) %>%
  group_by(L, N, n, Rep, group2) %>%
  summarise(meanQ1 = mean(q1)) %>%
  filter(meanQ1 == max(meanQ1)) %>%
  mutate(flip_it = (group2 == "ref-1")) %>%
  ungroup() %>%
  select(-group2, -meanQ1)
  


unsupes <- unsupes0 %>%
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


  
  
unsup_plot <- ggplot(unsupes, aes(x = simQ, y = q1, colour = factor(n))) +
  #geom_abline(intercept = -1/8, slope = 1/8, linetype = "dotted") +
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
  

ggsave(unsup_plot, filename = "results/figures/unsupervised_Qs.pdf", width = 10, height = 6)


#### Now do the sim_the_knowns ####

# these are the cases in which all the samples were re-simulated, including
# the "reference" one in either the simulated or unsimulated cases.

stk_supes <- sim_the_knowns %>%
  filter(scond == "supervised" & !str_detect(group, "reference")) %>%
  mutate(
    simQ = as.numeric(str_replace(group, "^q-", ""))
  ) %>%
  mutate(
    simQ = factor(simQ, levels = sort(unique(simQ)))
  )


stk_sup_plot <- ggplot(stk_supes, aes(x = simQ, y = q1, colour = factor(n))) +
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

ggsave(stk_sup_plot, filename = "results/figures/sim_the_knowns_supervised_Qs.pdf", width = 10, height = 6)









stk_unsupes0 <- sim_the_knowns %>%
  filter(scond == "unsupervised") %>%
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
stk_flippers <- stk_unsupes0 %>%
  filter(str_detect(group2, "ref")) %>%
  group_by(L, N, n, Rep, group2) %>%
  summarise(meanQ1 = mean(q1)) %>%
  filter(meanQ1 == max(meanQ1)) %>%
  mutate(flip_it = (group2 == "ref-1")) %>%
  ungroup() %>%
  select(-group2, -meanQ1)



stk_unsupes <- stk_unsupes0 %>%
  left_join(stk_flippers, by = c("L", "N", "n", "Rep")) %>%
  mutate(
    nq1 = ifelse(flip_it, q2, q1),
    nq2 = ifelse(flip_it, q1, q2),
    q1 = nq1,
    q2 = nq2
  ) %>%
  mutate(
    simQ = fct_recode(simQ, `ref-1` = "ref-2", `ref-2` = "ref-1" )
  )




stk_unsup_plot <- ggplot(stk_unsupes, aes(x = simQ, y = q1, colour = factor(n))) +
  #geom_abline(intercept = -1/8, slope = 1/8, linetype = "dotted") +
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


ggsave(stk_unsup_plot, filename = "results/figures/sim_the_knowns_unsupervised_Qs.pdf", width = 10, height = 6)



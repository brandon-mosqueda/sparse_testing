rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

source("~/Desktop/notes/utils.R")
source("~/Desktop/notes/genomic_plots.R")

Results <- merge_results("results", "summary", data_process_callback = function(Data, file) {
    Data$File <- file
    Data
  }) %>%
  filter(Fold != "Global")

Global <- Results %>%
  group_by(File) %>%
  reframe(global_row_summary(tibble(APC, NRMSE))) %>%
  mutate(Fold = "Global")

Results <- filter(Results, Fold == 1) %>%
  select(-APC, -NRMSE) %>%
  merge(Global, by = "File") %>%
  select(
    Dataset,
    Trait,
    Model,
    TestingProportion,
    WithInteraction,
    AnalysisType,
    APC,
    APC_SE,
    NRMSE,
    NRMSE_SE
  ) %>%
  mutate(Interaction = ifelse(WithInteraction, "Interaction", "No Interaction")) %>%
  relocate(Interaction, .before = APC)

generate_genomic_plots(
  Results,
  x = "Model",
  metrics = c("APC", "NRMSE"),
  facet_col = "Interaction",
  group_vars = "TestingProportion",
  results_dir = "plots"
)

write_csv(Results, file = "results/results.csv")

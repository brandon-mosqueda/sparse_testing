rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

source("~/Desktop/notes/utils.R")
source("~/Desktop/notes/genomic_plots.R")

Results <- merge_results("results", "summary",
  data_process_callback = function(Data, file) {
    Data$File <- file
    Data
  }) %>%
  filter(Fold != "Global")

# We recompute the global summary because in the original, the intervals SD were
# not computed.
Global <- Results %>%
  group_by(File) %>%
  reframe(global_row_summary(tibble(APC, NRMSE))) %>%
  mutate(Fold = "Global")

# We keep only the first Fold to extract the metadata and merge it with the new
# Global summary
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
  mutate(Interaction = ifelse(
    WithInteraction,
    "Interaction",
    "No Interaction"
  )) %>%
  relocate(Interaction, .before = APC)

generate_genomic_plots(
  Results,
  x = "Model",
  metrics = c("APC", "NRMSE"),
  facet_col = "Interaction",
  group_vars = c("Dataset", "TestingProportion"),
  results_dir = "plots"
)

Results %>%
  select(-AnalysisType) %>%
  round_df(4) %>%
  write_csv("results/results.csv")

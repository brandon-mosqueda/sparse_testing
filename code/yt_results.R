rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

source("~/Desktop/notes/utils.R")
source("~/Desktop/notes/genomic_plots.R")

mkdir("results/yt_eyt")

# Individual predictions -------------------------------------------------------
Predictions <- merge_results("results", "YT_EYT.+predictions.csv",
  data_process_callback = function(Data, file) {
    summary_file <- file.path(dirname(file), "summary.csv")
    file_sep <- str_split(file, "/")[[1]]

    gs_summaries(Data) %>%
      pluck("env") %>%
      mutate(
        Dataset = file_sep[3],
        Trait = file_sep[4],
        Model = file_sep[2],
        TestingProportion = file_sep[7],
        Interaction = file_sep[5]
      ) %>%
      write_csv(summary_file)

    Data %>%
      mutate(
        Dataset = file_sep[3],
        Trait = file_sep[4],
        Model = file_sep[2],
        TestingProportion = file_sep[7],
        Interaction = file_sep[5]
      )
  }) %>%
  select(
    Dataset,
    Trait,
    Model,
    TestingProportion,
    Interaction,
    Fold,
    Env,
    Line,
    Predicted,
    Observed
  ) %>%
  round_df(4)

write_csv(Predictions, "results/yt_eyt/all_predictions.csv")

# Summary by env ---------------------------------------------------------------
ByEnv <- merge_results("results", "YT_EYT.+summary.csv") %>%
  select(
    Dataset,
    Trait,
    Model,
    TestingProportion,
    Interaction,
    Env,
    APC,
    APC_SE,
    NRMSE,
    NRMSE_SE
  )

ByEnv %>%
  round_df(4) %>%
  write_csv("results/yt_eyt/by_env.csv")

ByEnv %>%
  filter(!(Env %in% c("Across-env", "Across-env-new"))) %>%
  generate_genomic_plots(
    x = "Model",
    fill_by = "Env",
    metrics = c("APC", "NRMSE"),
    facet_col = "Interaction",
    group_vars = c("Dataset", "TestingProportion"),
    results_dir = "plots/by_env"
  )

ByEnv %>%
  filter(Env %in% c("Across-env", "Across-env-new")) %>%
  generate_genomic_plots(
    x = "Model",
    fill_by = "Env",
    metrics = c("APC", "NRMSE"),
    facet_col = "Interaction",
    group_vars = c("Dataset", "TestingProportion"),
    results_dir = "plots/across_env"
  )

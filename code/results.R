rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

source("~/Desktop/notes/utils.R")
source("~/Desktop/notes/genomic_plots.R")

# Individual predictions -------------------------------------------------------
Predictions <- merge_results("results", "predictions.csv",
  data_process_callback = function(Data, file) {
    Summary <- read_csv(file.path(dirname(file), "summary.csv"))

    Data$Dataset <- Summary$Dataset[1]
    Data$Trait <- Summary$Trait[1]
    Data$Model <- Summary$Model[1]
    Data$TestingProportion <- Summary$TestingProportion[1]
    Data$WithInteraction <- as.logical(Summary$WithInteraction[1])

    Data
  }) %>%
  mutate(
    Interaction = ifelse(
      WithInteraction,
      "Interaction",
      "No Interaction"
    )
  ) %>%
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

library(openxlsx)

# Create a work book
book <- createWorkbook(title = "All predictions")

addWorksheet(book, "YT_23-24")
addWorksheet(book, "YT_22-23")

writeData(
  book,
  "YT_23-24",
  Predictions %>% filter(Dataset == "YT_23-24"),
  rowNames = FALSE
)
writeData(
  book,
  "YT_22-23",
  Predictions %>% filter(Dataset == "YT_22-23"),
  rowNames = FALSE
)

# Save workbook
saveWorkbook(book, file = "results/all_predictions.xlsx", overwrite = TRUE)

# Summary by env ---------------------------------------------------------------
ByEnv <- merge_results("results", "predictions.csv",
  data_process_callback = function(Predictions, file) {
    OldSummary <- read_csv(file.path(dirname(file), "summary.csv"))

    Summary <- Predictions %>%
      group_by(Env, Fold) %>%
      summarise(
        APC = pearson(Predicted, Observed),
        NRMSE = nrmse(Observed, Predicted)
      ) %>%
      group_by(Env) %>%
      summarise_if(
        is.numeric,
        list(
          MEAN = ~ mean(., na.rm = TRUE),
          SE = ~ sd(., na.rm = TRUE) / sqrt(n())
        )
      ) %>%
      select(-contains("Fold")) %>%
      rename_with(~ gsub("_MEAN", "", .)) %>%
      mutate(
        Dataset = OldSummary$Dataset[1],
        Trait = OldSummary$Trait[1],
        Model = OldSummary$Model[1],
        TestingProportion = OldSummary$TestingProportion[1],
        WithInteraction = as.logical(OldSummary$WithInteraction[1])
      )

    Summary
  })

ByEnv <- ByEnv %>%
  mutate(Interaction = ifelse(
    WithInteraction,
    "Interaction",
    "No Interaction"
  )) %>%
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

generate_genomic_plots(
  ByEnv,
  x = "Model",
  fill_by = "Env",
  metrics = c("APC", "NRMSE"),
  facet_col = "Interaction",
  group_vars = c("Dataset", "TestingProportion"),
  results_dir = "plots/by_env"
)

ByEnv %>%
  round_df(4) %>%
  write_csv("results/by_env.csv")

# Across fold ------------------------------------------------------------------
AcrossFold <- merge_results("results", "predictions.csv",
  data_process_callback = function(Predictions, file) {
    OldSummary <- read_csv(file.path(dirname(file), "summary.csv"))

    Summary <- Predictions %>%
      group_by(Fold) %>%
      summarise(
        APC = pearson(Predicted, Observed),
        NRMSE = nrmse(Observed, Predicted)
      ) %>%
      summarise_if(
        is.numeric,
        list(
          MEAN = ~ mean(., na.rm = TRUE),
          SE = ~ sd(., na.rm = TRUE) / sqrt(n())
        )
      ) %>%
      select(-contains("Fold")) %>%
      rename_with(~ gsub("_MEAN", "", .)) %>%
      mutate(
        Dataset = OldSummary$Dataset[1],
        Trait = OldSummary$Trait[1],
        Model = OldSummary$Model[1],
        TestingProportion = OldSummary$TestingProportion[1],
        WithInteraction = as.logical(OldSummary$WithInteraction[1])
      )

    Summary
  }
)

AcrossFold <- AcrossFold %>%
  mutate(Interaction = ifelse(
    WithInteraction,
    "Interaction",
    "No Interaction"
  )) %>%
  select(
    Dataset,
    Trait,
    Model,
    TestingProportion,
    Interaction,
    APC,
    APC_SE,
    NRMSE,
    NRMSE_SE
  )

generate_genomic_plots(
  AcrossFold,
  x = "Model",
  metrics = c("APC", "NRMSE"),
  facet_col = "Interaction",
  group_vars = c("Dataset", "TestingProportion"),
  results_dir = "plots/across_fold"
)

AcrossFold %>%
  round_df(4) %>%
  write_csv("results/across_fold.csv")

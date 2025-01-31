rm(list = ls())

library(SKM)
library(dplyr)
library(doParallel)
library(foreach)

setwd("~/data_science/sparse_testing")

source("code/utils.R")
source("code/cross_validators_env.R")

# Globals definition --------------------------------------------------
with_interaction <- FALSE
model_name <- "M1"

dataset <- "YT_EYT_22_23"
testing_proportions <- c(0.15, 0.25, 0.5, 0.75, 0.85)
iterations_number <- 10000
burn_in <- 5000
base_results_dir <- "results"

# dataset <- "Test_YT_EYT"
# testing_proportions <- c(0.5)
# iterations_number <- 10
# burn_in <- 5
# base_results_dir <- "trash"

# Data preparation -------------------------------------------------------------
load(sprintf("data/%s.RData", dataset), verbose = TRUE)

Pheno <- Pheno %>%
  arrange(Env, Line) %>%
  droplevels() %>%
  mutate(Index = seq(n()))

Env <- model.matrix(~ 0 + Env, data = Pheno)
Line <- model.matrix(~ 0 + Line, data = Pheno)

lines <- levels(Pheno$Line)
Geno <- cholesky(Geno[lines, lines])
GenoLine <- Line %*% Geno

ETA <- list(
  Env = list(X = Env, model = "FIXED"),
  Line = list(X = GenoLine, model = "BRR")
)

if (with_interaction) {
  ETA$LinexEnv <- list(
    X = model.matrix(~ 0 + GenoLine:Env),
    model = "BRR"
  )
}

y <- Pheno$GY
# The trails with EYT are always included in the training set
testing_indices <- stringr::str_which(Pheno$Trial, "EYT", negate = TRUE)
# PhenoTesting is only used to generate the partitions
PhenoTesting <- Pheno %>% slice(testing_indices) %>% droplevels()

cluster <- makeCluster(10, outfile = "")
registerDoParallel(cluster)

# Model evaluation -------------------------------------------------------------
for (testing_proportion in testing_proportions) {
  echo("*** Testing proportion: %s ***", testing_proportion)

  results_dir <- file.path(
    base_results_dir,
    model_name,
    data_info$name,
    "GY",
    if (with_interaction) "with_interaction" else "no_interaction",
    "unitrait",
    testing_proportion
  )

  echo("\t-------- Generating folds --------")

  training_proportion <- 1 - testing_proportion
  training_lines_num <- ceiling(lunique(PhenoTesting$Line) * training_proportion)

  folds <- folds_by_model(
    model = model_name,
    lines = PhenoTesting$Line,
    envs = PhenoTesting$Env,
    folds_num = 10,
    training_lines_num = training_lines_num,
    training_proportion = training_proportion,
    verbose = TRUE
  )

  mkdir(results_dir)
  save(folds, file = file.path(results_dir, "folds.RData"))

  echo("\t-------- Cross validation --------")

  Predictions <- foreach(
    i = seq_along(folds),
    .combine = rbind,
    .packages = "dplyr"
  )  %dopar% {
    SKM::echo("\t\t*** Fold: %s ***", i)
    fold <- folds[[i]]

    # This indices are in terms of PhenoTesting
    testing_indices <- env_testing_indices(fold)
    # We extract the real indices from PhenoTesting
    final_testing_indices <- PhenoTesting %>%
      slice(testing_indices) %>%
      pull(Index)

    y_na <- replace(y, final_testing_indices, NA)

    model <- BGLR::BGLR(
      y_na,
      ETA = ETA,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- model$yHat[testing_indices]

    data.frame(
      Fold = i,
      Line = Pheno$Line[testing_indices],
      Env = Pheno$Env[testing_indices],
      Predicted = predictions,
      Observed = y[testing_indices]
    )
  }

  summaries <- gs_summaries(Predictions)

  # Join the summaries
  Summary <- summaries$env %>%
    mutate(Type = "Summary 1") %>%
    bind_rows(
      summaries$env2 %>%
        mutate(Type = "Summary 2")
    )

  write_csv(Summary, file = file.path(results_dir, "summary.csv"))
  write_csv(Predictions, file = file.path(results_dir, "predictions.csv"))
}

stopCluster(cluster)

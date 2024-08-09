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

dataset <- "YT_22-23"
testing_proportions <- c(0.15, 0.25, 0.5, 0.75, 0.85)
iterations_number <- 10000
burn_in <- 5000
base_results_dir <- "results"

# dataset <- "Test"
# testing_proportions <- c(0.5)
# base_results_dir <- "trash"
# iterations_number <- 10
# burn_in <- 5

# Data preparation --------------------------------------------------
load(sprintf("data/%s.RData", dataset), verbose = TRUE)

Pheno <- Pheno %>% arrange(Env, Line) %>% droplevels()
Env <- model.matrix(~ 0 + Env, data = Pheno)
Line <- model.matrix(~ 0 + Line, data = Pheno)

geno_lines <- levels(Pheno$Line)
Geno <- cholesky(Geno[geno_lines, geno_lines])
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

y <- Pheno[[data_info$traits]]

cluster <- makeCluster(10, outfile = "")
registerDoParallel(cluster)

# Model evaluation --------------------------------------------------
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
  training_lines_num <- ceiling(lunique(Pheno$Line) * training_proportion)

  folds <- folds_by_model(
    model = model_name,
    lines = Pheno$Line,
    envs = Pheno$Env,
    folds_num = 10,
    training_lines_num = training_lines_num,
    training_proportion = training_proportion,
    verbose = TRUE
  )

  mkdir(results_dir)
  save(folds, file = file.path(results_dir, "folds.RData"))

  echo("\t-------- Cross validation --------")

  Predictions <- data.frame()

  Predictions <- foreach(
    i = seq_along(folds),
    .combine = rbind,
    .packages = "dplyr"
  )  %dopar% {
    SKM::echo("\t\t*** Fold: %s ***", i)
    fold <- folds[[i]]

    testing_indices <- env_testing_indices(fold)
    y_na <- replace(y, testing_indices, NA)

    temp_dir <- file.path(tempdir(), i, runif(1))
    SKM::mkdir(temp_dir)

    model <- BGLR::BGLR(
      y_na,
      ETA = ETA,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = temp_dir
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

  Summary <- Predictions %>%
    group_by(Fold) %>%
    summarise(
      APC = pearson(Predicted, Observed),
      NRMSE = nrmse(Observed, Predicted)
    ) %>%
    mutate(Fold = as.character(Fold))

  Global <- Summary %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    mutate(Fold = "Global")

  Summary <- bind_rows(Summary, Global) %>%
    mutate(
      Dataset = data_info$name,
      Trait = data_info$traits,
      Model = model_name,
      TestingProportion = testing_proportion,
      WithInteraction = with_interaction,
      AnalysisType = "Unitrait"
    )

  write_csv(Summary, file = file.path(results_dir, "summary.csv"))
  write_csv(Predictions, file = file.path(results_dir, "predictions.csv"))
}

stopCluster(cluster)

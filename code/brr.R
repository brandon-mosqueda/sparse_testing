rm(list = ls())

library(SKM)
library(dplyr)
library(doParallel)
library(foreach)

setwd("~/data_science/sparse_testing")

source("code/utils.R")

# Globals definition --------------------------------------------------
iterations_number <- 1000
burn_in <- 500
dataset <- "YT_22-23"

# iterations_number <- 10
# burn_in <- 5
# dataset <- "Test"

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

y <- Pheno[[data_info$traits]]

cluster <- makeCluster(5, outfile = "")
registerDoParallel(cluster)

folds <- cv_random_line(Pheno$Line, 10, 0.3)

Predictions <- foreach(
  fold = folds,
  .combine = rbind,
  .packages = "dplyr"
) %dopar% {
  SKM::echo("\t\t*** Fold: %s ***", fold$num)
  y_na <- replace(y, fold$testing, NA)

  temp_dir <- file.path(tempdir(), fold$num, runif(1))
  SKM::mkdir(temp_dir)

  model <- BGLR::BGLR(
    y_na,
    ETA = ETA,
    nIter = iterations_number,
    burnIn = burn_in,
    verbose = FALSE,
    saveAt = temp_dir
  )
  predictions <- model$yHat[fold$testing]

  data.frame(
    Fold = fold$num,
    Line = Pheno$Line[fold$testing],
    Env = Pheno$Env[fold$testing],
    Predicted = predictions,
    Observed = y[fold$testing]
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

Summary <- bind_rows(Summary, Global)

write_csv(Summary, file = "results/test.csv")
write_csv(Predictions, file = "results/test_predictions.csv")

stopCluster(cluster)

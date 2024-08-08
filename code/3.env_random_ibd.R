rm(list = ls())

library(SKM)
library(dplyr)

source("utils.R")
source("cross_validators_env.R")

load("data/SNP.RData", verbose = TRUE)

# Globals definition --------------------------------------------------
with_interaction <- FALSE
model_name <- "BRR"
cv_method <- "M3"
testing_proportions <- c(0.1, 0.25, 0.5, 0.75, 0.85)

digits <- 4

# Data preparation --------------------------------------------------
Pheno <- Pheno %>% arrange(Env, Line)
Env <- model.matrix(~ 0 + Env, data = Pheno)
Line <- model.matrix(~ 0 + Line, data = Pheno)

geno_lines <- sort(rownames(Geno))
Geno <- Geno[geno_lines, geno_lines] %>% cholesky()
GenoLine <- Line %*% Geno

LinexEnv <- model.matrix(~ 0 + GenoLine:Env)

ETA <- list(
  list(x = Env, model = "fixed"),
  list(x = GenoLine, model = model_name),
  list(x = LinexEnv, model = model_name)
)
Y <- Pheno[, data_info$responses, drop = FALSE]

# Model evaluation --------------------------------------------------
for (testing_proportion in testing_proportions) {
  echo("*** Testing proportion: %s ***", testing_proportion)

  results_dir <- file.path(
    "results",
    "env",
    cv_method,
    data_info$name,
    model_name,
    if (with_interaction) "with_interaction" else "no_interaction",
    "multitrait",
    testing_proportion
  )

  echo("\t-------- Generating folds --------")

  training_proportion <- 1 - testing_proportion

  folds <- cv_env_random_ibd(
    lines = Pheno$Line,
    envs = Pheno$Env,
    training_proportion = training_proportion,
    folds_num = 10,
    verbose = TRUE
  )

  mkdir(results_dir)
  save(folds, file = file.path(results_dir, "folds.RData"))

  echo("\t-------- Cross validation --------")

  Results <- data.frame()

  for (i in seq_along(folds)) {
    echo("\t\t*** Fold: %s ***", i)
    fold <- folds[[i]]

    testing_indices <- env_testing_indices(fold)
    Y_NA <- Y
    Y_NA[testing_indices, , drop = FALSE] <- NA

    model <- bayesian_model(
      ETA,
      Y_NA,
      iterations_number = 10000,
      burn_in = 5000
    )
    predictions <- predict(model)

    FoldSummary <- env_multivariate_summary(
        Y[testing_indices, ],
        predictions,
        digits
      ) %>%
      mutate(Fold = i)

    Results <- rbind(Results, FoldSummary)
  }

  write_csv(Results, file = file.path(results_dir, "predictions.csv"))

  Results <- prepare_results_by_fold(
    Results = Results,
    dataset = data_info$name,
    model_name = model_name,
    cv_method = cv_method,
    testing_proportion = testing_proportion,
    with_interaction = with_interaction,
    analysis_type = "Multitrait"
  )

  write_csv(Results, file = file.path(results_dir, "by_fold.csv"))

  Global <- multitrait_global_env_summary(Results, digits)

  write_csv(Global, file = file.path(results_dir, "global.csv"))
}

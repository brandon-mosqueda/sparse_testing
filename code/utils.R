library(dplyr)
library(SKM)

# Multitrait --------------------------------------------------

testing_to_NA <- function(Pheno, fold) {
  for (trait in names(fold)) {
    testing_indices <- fold[[trait]]$testing

    Pheno[[trait]][testing_indices] <- NA
  }

  return(Pheno)
}

env_multivariate_summary <- function(observed, predicted, digits) {
  Summary <- data.frame()

  for (trait in names(observed)) {
    summary <- numeric_summary(
        observed[[trait]],
        predicted[[trait]]$predicted
      ) %>%
      as.data.frame.list() %>%
      mutate(
        Trait = trait,
        nrmse = nrmse(
          observed[[trait]],
          predicted[[trait]]$predicted,
          type = "mean"
        )
      ) %>%
      relocate(Trait, 1)

    Summary <- rbind(Summary, summary)
  }

  Summary <- Summary %>% round_df(digits)

  return(Summary)
}

print.CVMultivariate <- function(folds) {
  folds_num <- length(folds)
  traits <- names(folds[[1]])
  traits_num <- length(traits)

  for (i in seq(folds_num)) {
    cat("*** Fold", i, "***\n")
    lines <- c(
      attr(folds[[i]][[1]]$training, "lines"),
      attr(folds[[i]][[1]]$testing, "lines")
    )
    lines_num <- length(lines)

    Matrix <- matrix(1, lines_num, traits_num)
    rownames(Matrix) <- lines
    colnames(Matrix) <- traits

    fold <- folds[[i]]

    for (trait in traits) {
      Matrix[attr(fold[[trait]]$testing, "lines"), trait] <- 0
    }

    print(Matrix)
    cat("\n")
  }
}

multitrait_global_env_summary <- function(Results, digits) {
  traits <- unique(Results$Trait)

  Summary <- data.frame()

  for (trait in traits) {
    Data <- Results %>%
      filter(Trait == trait) %>%
      droplevels()

    Summary <- Summary %>%
      rbind(global_env_summary(Data, digits))
  }

  return(Summary)
}

# Utils --------------------------------------------------

divide <- function(n, divisor) {
  times <- rep(floor(n / divisor), divisor)
  module <- n %% divisor

  if (module > 0) {
    increase_one_indices <- seq(1, module)
    times[increase_one_indices] <- times[1] + 1
  }

  return(times)
}

env_testing_indices <- function(fold) {
  result <- c()
  for (env in names(fold)) {
    result <- c(result, fold[[env]]$testing)
  }

  return(result)
}

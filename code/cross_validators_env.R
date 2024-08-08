random_line <- function(lines, envs, training_lines_num) {
  unique_lines_num <- length(unique(lines))
  envs_num <- length(envs)

  training_lines <- sample(unique_lines_num, training_lines_num)

  TrainingMatrix <- matrix(
    rep(training_lines, envs_num),
    ncol = training_lines_num,
    nrow = envs_num,
    byrow = TRUE
  )

  return(TrainingMatrix)
}

partial_random_line <- function(lines,
                                envs,
                                training_lines_num) {
  unique_lines <- unique(lines)
  unique_lines_num <- length(unique_lines)
  envs_num <- length(envs)
  testing_lines_num <- unique_lines_num - training_lines_num

  training_lines <- sample(unique_lines_num, training_lines_num)
  testing_lines <- setdiff(seq(unique_lines_num), training_lines)

  fold <- list()
  # It is assumed the first env has the correct number of testing lines
  fold[[envs[1]]] <- list(training = training_lines, testing = testing_lines)

  to_change_testing_elements <- divide(testing_lines_num, envs_num - 1)
  # The first NA is because the first environment is completely assigned
  # but this value is not used
  to_change_testing_elements <- c(NA, to_change_testing_elements)

  max_change_testing <- max(to_change_testing_elements, na.rm = TRUE)
  if (training_lines_num > max_change_testing) {
    common_training_lines <- seq(
      1,
      training_lines_num - max_change_testing
    )

    training_lines <- training_lines[common_training_lines]
  } else {
    training_lines <- NULL
  }

  index <- 1

  for (i in seq(2, envs_num)) {
    env_testing_lines <- testing_lines[
      seq(index, index + to_change_testing_elements[i] - 1)
    ]
    fold[[envs[i]]] <- list(training = c(training_lines, env_testing_lines))
    fold[[envs[i]]]$testing <- setdiff(
      seq(unique_lines_num),
      fold[[envs[i]]]$training
    )

    index <- index + to_change_testing_elements[i]
  }

  return(fold)
}

random_ibd <- function(lines, envs, training_proportion = 0.5) {
  unique_lines_num <- length(unique(lines))
  unique_envs_num <- length(unique(envs))
  envs_num <- length(envs)
  k <- ceiling(unique_lines_num * training_proportion)
  # Equivalent to
  # k <- ceiling(
  #   (unique_lines_num  * unique_envs_num * training_proportion) /
  #   unique_envs_num
  # )

  BlocksDesign <- matrix(NA, nc = k, nr = envs_num)

  line_frequency <- ceiling(k * unique_envs_num / unique_lines_num)
  lines_frequencies <- rep(line_frequency, unique_lines_num)

  for (env_num in seq(envs_num)) {
    if (sum(lines_frequencies > 0) >= k) {
      env_block <- sample(unique_lines_num, k, prob = lines_frequencies)
      BlocksDesign[env_num, ] <- env_block

      lines_frequencies[env_block] <- lines_frequencies[env_block] - 1
    } else {
      # In last fold
      if (envs_num - env_num == 0) {
        BlocksDesign[env_num, ] <- sample(unique_lines_num, k)
      } else {
        stop("Not enough GIDs")
      }
    }
  }

  return(BlocksDesign)
}

generate_multivariate_folds <- function(lines,
                                        envs,
                                        type = "ibd",
                                        folds_num = 10,
                                        training_proportion = 0.5,
                                        training_lines_num = length(lines) / 2,
                                        verbose = TRUE,
                                        seed = NULL) {
  # Set seed for reproducible results preserving the previous random state
  old_random_state <- NULL
  if (!is.null(seed)) {
    old_random_state <- get_rand_state()

    set.seed(seed)
  }
  on.exit(set_rand_state(old_random_state))

  type <- tolower(type)
  lines <- as.character(lines)
  unique_lines <- unique(lines)
  unique_lines_num <- length(unique_lines)

  envs <- as.character(envs)
  unique_envs <- unique(envs)
  unique_envs_num <- length(unique_envs)

  k <- ceiling(unique_lines_num * training_proportion)
  # Equivalent to
  # k <- ceiling(
  #   (unique_lines_num  * unique_envs_num * training_proportion) /
  #   unique_envs_num
  # )

  folds <- list()

  for (fold_num in seq(folds_num)) {
    if (verbose) {
      cat("*** Fold", fold_num, "***\n")
    }
    folds[[fold_num]] <- list()

    if (type == "ibd") {
      blocks_design <- crossdes::find.BIB(
        trt = unique_lines_num,
        b = unique_envs_num,
        k = k
      )
    } else if (type == "random_ibd") {
      blocks_design <- random_ibd(lines, unique_envs, training_proportion)
    } else if (type == "random_line") {
      blocks_design <- random_line(lines, unique_envs, training_lines_num)
    }

    for (env_num in seq(unique_envs_num)) {
      env <- unique_envs[env_num]
      env_indices <- which(envs == env)

      training_lines <- unique_lines[blocks_design[env_num, ]]
      training_lines_indices <- which(lines %in% training_lines)
      training_indices <- intersect(training_lines_indices, env_indices)
      attr(training_indices, "lines") <- training_lines

      testing_lines <- setdiff(unique_lines, training_lines)
      testing_lines_indices <- which(lines %in% testing_lines)
      testing_indices <- intersect(testing_lines_indices, env_indices)
      attr(testing_indices, "lines") <- testing_lines

      folds[[fold_num]][[env]] <- list(
        training = training_indices,
        testing = testing_indices
      )
    }
  }

  class(folds) <- "CVMultivariate"

  return(folds)
}

cv_env_ibd <- function(lines,
                       envs,
                       folds_num = 10,
                       training_proportion = 2,
                       verbose = TRUE,
                       seed = NULL) {
  return(generate_multivariate_folds(
    lines = lines,
    envs = envs,
    type = "ibd",
    folds_num = folds_num,
    training_proportion = training_proportion,
    verbose = verbose,
    seed = seed
  ))
}

cv_env_random_ibd <- function(lines,
                              envs,
                              folds_num = 10,
                              training_proportion = 2,
                              verbose = TRUE,
                              seed = NULL) {
  return(generate_multivariate_folds(
    lines = lines,
    envs = envs,
    type = "random_ibd",
    folds_num = folds_num,
    training_proportion = training_proportion,
    verbose = verbose,
    seed = seed
  ))
}

cv_env_random_line <- function(lines,
                               envs,
                               training_lines_num,
                               folds_num = 10,
                               verbose = TRUE,
                               seed = NULL) {
  return(generate_multivariate_folds(
    lines = lines,
    envs = envs,
    type = "random_line",
    folds_num = folds_num,
    training_lines_num = training_lines_num,
    verbose = verbose,
    seed = seed
  ))
}

cv_env_partial_random_line <- function(lines,
                                       envs,
                                       training_lines_num,
                                       folds_num = 10,
                                       verbose = TRUE,
                                       seed = NULL) {
  # Set seed for reproducible results preserving the previous random state
  old_random_state <- NULL
  if (!is.null(seed)) {
    old_random_state <- get_rand_state()

    set.seed(seed)
  }
  on.exit(set_rand_state(old_random_state))

  lines <- as.character(lines)
  unique_lines <- unique(lines)

  envs <- as.character(envs)
  unique_envs <- unique(envs)

  folds <- list()

  for (fold_num in seq(folds_num)) {
    if (verbose) {
      cat("*** Fold", fold_num, "***\n")
    }
    folds[[fold_num]] <- list()
    fold <- partial_random_line(lines, unique_envs, training_lines_num)

    for (env in names(fold)) {
      env_indices <- which(envs == env)

      training_lines <- unique_lines[fold[[env]]$training]
      training_lines_indices <- which(lines %in% training_lines)
      training_indices <- intersect(training_lines_indices, env_indices)
      attr(training_indices, "lines") <- training_lines

      testing_lines <- setdiff(unique_lines, training_lines)
      testing_lines_indices <- which(lines %in% testing_lines)
      testing_indices <- intersect(testing_lines_indices, env_indices)
      attr(testing_indices, "lines") <- testing_lines

      folds[[fold_num]][[env]] <- list(
        training = training_indices,
        testing = testing_indices
      )
    }
  }

  class(folds) <- "CVMultivariate"

  return(folds)
}

cv_augmented_ibd <- function(lines,
                             envs,
                             training_proportion = 0.5,
                             folds_num = 10,
                             tolerance = 10,
                             verbose = TRUE,
                             seed = NULL) {
  # Set seed for reproducible results preserving the previous random state
  old_random_state <- NULL
  if (!is.null(seed)) {
    old_random_state <- get_rand_state()

    set.seed(seed)
  }
  on.exit(set_rand_state(old_random_state))

  lines <- as.character(lines)
  unique_lines <- unique(lines)
  unique_lines_num <- length(unique_lines)

  envs <- as.character(envs)
  unique_envs <- unique(envs)
  unique_envs_num <- length(unique_envs)

  testing_proportion <- 1 - training_proportion
  training_length <- unique_lines_num * unique_envs_num * training_proportion
  ko <- round(training_length / unique_envs_num)

  # El numero de lines por este numero no sea mayor que 100
  fictitious_proportion_checks <- min(100 / unique_lines_num, 1)
  tolerance_fictitious_proportion <- seq(100, 10, length.out = tolerance) /
    unique_lines_num
  tolerance_fictitious_proportion <- ifelse(
    tolerance_fictitious_proportion > 1,
    1,
    tolerance_fictitious_proportion
  )

  # No cambiar
  fictitious_proportion_checks_training <- 0.5

  folds <- list()

  for (fold_num in seq(folds_num)) {
    if (verbose) {
      cat("*** Fold", fold_num, "***\n")
    }

    blocks_design <- NULL
    i <- 1
    while (i != tolerance) {
      i <- i + 1
      fictitious_lines <- sample(
        unique_lines,
        size = round(fictitious_proportion_checks * unique_lines_num)
      )

      fictitious_lines_num <- length(fictitious_lines)
      no_fictitious_lines <- setdiff(lines, fictitious_lines)

      fictitious_checks_training_num <- round(
        fictitious_lines_num *
          unique_envs_num *
          fictitious_proportion_checks_training
      )
      k <- round(fictitious_checks_training_num / unique_envs_num)

      tryCatch({
        blocks_design <<- crossdes::find.BIB(
          trt = fictitious_lines_num,
          b = unique_envs_num,
          k = k
        )
        break
      }, error = function(error) {
        echo("Error with crossdes::find.BIB")
        print(error)
        if (i == tolerance) {
          stop("Maximum tolerance reached")
        }
        fictitious_proportion_checks <<- tolerance_fictitious_proportion[i]
      })
      unique_lines <- fictitious_lines
    }

    no_fictitious_lines_x2 <- c(
      sample(no_fictitious_lines),
      sample(no_fictitious_lines)
    )
    to_agument_num <- ko - k

    folds[[fold_num]] <- list()

    for (i in seq(unique_envs_num)) {
      env <- unique_envs[i]
      env_indices <- which(envs == env)

      training_lines_augmented <- unique_lines[blocks_design[i, ]]

      if (!SKM::is_empty(no_fictitious_lines_x2)) {
        j <- max(to_agument_num * i, 1)
        k <- max(to_agument_num * (i - 1) + 1, 1)
        pos_aug <- seq(k, j)

        training_lines_augmented <- c(
          training_lines_augmented,
          no_fictitious_lines_x2[pos_aug]
        )
      }
      training_lines_indices <- which(lines %in% training_lines_augmented)
      training_indices <- intersect(training_lines_indices, env_indices)
      attr(training_indices, "lines") <- training_lines_augmented

      testing_lines <- setdiff(unique_lines, training_lines_augmented)
      testing_lines_indices <- which(lines %in% testing_lines)
      testing_indices <- intersect(testing_lines_indices, env_indices)
      attr(testing_indices, "lines") <- testing_lines

      folds[[fold_num]][[env]] <- list(
        training = training_indices,
        testing = testing_indices
      )
    }

    fictitious_proportion_checks <- min(100 / unique_lines_num, 1)
  }

  class(folds) <- "CVMultivariate"

  return(folds)
}

folds_by_model <- function(model,
                           lines,
                           envs,
                           folds_num,
                           training_lines_num,
                           training_proportion,
                           verbose) {
  model <- tolower(model)

  if (model == "m1") {
    folds <- cv_env_random_line(
      lines = Pheno$Line,
      envs = Pheno$Env,
      folds_num = 10,
      training_lines_num = training_lines_num,
      verbose = verbose
    )
  } else if (model == "m2") {
    folds <- cv_env_partial_random_line(
      lines = Pheno$Line,
      envs = Pheno$Env,
      folds_num = 10,
      training_lines_num = training_lines_num,
      verbose = verbose
    )
  } else if (model == "m3") {
    folds <- cv_env_random_ibd(
      lines = Pheno$Line,
      envs = Pheno$Env,
      training_proportion = training_proportion,
      folds_num = 10,
      verbose = verbose
    )
  } else if (model == "m4") {
    folds <- cv_augmented_ibd(
      lines = Pheno$Line,
      envs = Pheno$Env,
      training_proportion = training_proportion,
      folds_num = 10,
      verbose = verbose
    )
  } else {
    stop("Model not found")
  }

  return(folds)
}

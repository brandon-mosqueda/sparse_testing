rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

# YT_22-23 ---------------------------------------------------------------------
Pheno <- SKM::read_csv(
    "original_files/YT_22-23/input_pheno_YT_22-23.txt",
    sep = "\t"
  ) %>%
  pivot_longer(
    cols = c(B2IR, B5IR, BDRT, BLHT),
    names_to = "Env",
    values_to = "GY"
  ) %>%
  rename(Line = "GID") %>%
  mutate(Env = factor(Env), Line = factor(Line)) %>%
  group_by(Env) %>%
  # Impute some of the few missing values
  mutate(GY = replace(GY, is.na(GY), mean(GY, na.rm = TRUE))) %>%
  ungroup() %>%
  arrange(Env, Line)

Markers <- SKM::read_csv(
    "original_files/YT_22-23/YT_22-23_filtred_SNP.txt",
    sep = "\t"
  ) %>%
  rename(Line = "V1") %>%
  as.data.frame()

final_lines <- intersect(levels(Pheno$Line), Markers$Line)

Pheno <- Pheno %>% filter(Line %in% final_lines) %>% droplevels()
Markers <- Markers %>% filter(Line %in% final_lines)

rownames(Markers) <- as.character(Markers$Line)
Markers$Line <- NULL
Markers <- scale(as.matrix(Markers))

Geno <- tcrossprod(Markers) / ncol(Markers)

data_info <- list(
  name = "YT_22-23",
  traits = "GY"
)

save(Pheno, Geno, data_info, file = "data/YT_22-23.RData")

# YT_23-24 ---------------------------------------------------------------------

Pheno <- readxl::read_excel("data/YT_23-24.xlsx") %>%
  # Convert B2IR, B5IR, BDRT, BLHT columns to rows
  pivot_longer(
    cols = c(B2IR, B5IR, BDRT, BLHT),
    names_to = "Env",
    values_to = "GY"
  ) %>%
  rename(Line = "GID") %>%
  mutate(Env = factor(Env), Line = factor(Line)) %>%
  arrange(Env, Line)

Markers <- SKM::read_csv(
    "data/Numerical_input_YT_23_24.txt",
    sep = "\t"
  ) %>%
  rename(Line = "V1") %>%
  as.data.frame()

final_lines <- intersect(levels(Pheno$Line), Markers$Line)

Pheno <- Pheno %>% filter(Line %in% final_lines) %>% droplevels()
Markers <- Markers %>% filter(Line %in% final_lines)

rownames(Markers) <- as.character(Markers$Line)
Markers$Line <- NULL
Markers <- scale(as.matrix(Markers))

Geno <- tcrossprod(Markers) / ncol(Markers)

data_info <- list(
  name = "YT_23-24",
  traits = "GY"
)

save(Pheno, Geno, data_info, file = "data/YT_23-24.RData")

# Test -------------------------------------------------------------------------
load("data/YT_22-23.RData", verbose = TRUE)

Pheno <- Pheno %>%
  group_by(Env) %>%
  arrange(Line) %>%
  slice(seq(50)) %>%
  ungroup() %>%
  droplevels()

lines <- levels(Pheno$Line)
Geno <- Geno[lines, lines]

data_info$name <- "Test"
save(Pheno, Geno, data_info, file = "data/Test.RData")

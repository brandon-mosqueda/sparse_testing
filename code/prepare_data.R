rm(list = ls())

setwd("~/data_science/sparse_testing")

library(SKM)
library(tidyverse)

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

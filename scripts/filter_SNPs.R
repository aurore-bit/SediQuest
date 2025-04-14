library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
score_b <- as.numeric(args[1])  # The first argument is the score_b
path_to_burden <- args[2]
path_to_control <- args[3]
output <- args[4]
score_n <- as.numeric(args[5])

burden_score <- fread(path_to_burden) %>%
   dplyr::filter(b_3 != ".") %>%
   dplyr::mutate(b_3 = as.numeric(as.character(b_3))) %>%
   dplyr::filter(b_3 >= score_b) %>%
   dplyr::mutate(n_3 = as.numeric(as.character(n_3))) %>%
   dplyr::filter(n_3 >= score_n)

new_bed <- fread(path_to_control) %>%
   dplyr::rename(chrom = V1, pos0 = V2, pos = V3) %>%
  inner_join(burden_score, by = c("chrom", "pos0", "pos")) %>%
   dplyr::select(chrom, pos0, pos, V4, V5, V6, V7, V8, V9)

fwrite(new_bed, output, sep = "\t", col.names = FALSE)

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
score_b <- as.numeric(args[1])  # The first argument is the score_b
temp_sites <- args[2]
path_to_burden <- args[3]
path_to_control <- args[4]

burden_score <- fread(path_to_burden) %>%
  filter(b_3 != ".") %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 >= score_b)

new_bed <- fread(path_to_control) %>%
  rename(pos0 = V2,
         pos = V3) %>%
  filter(pos0 %in% burden_score$pos0 & pos %in% burden_score$pos)

fwrite(new_bed, temp_sites, sep = "\t", col.names = FALSE)

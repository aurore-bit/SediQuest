
#TO HAVE MD score ACCROSS BURDEN SCORES#

library(knitr)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggplot2)

#install.packages("soiltestcorr")

library(soiltestcorr)



args <- commandArgs(trailingOnly = TRUE)
coverage_path <- args[1]
burden_path <- args[2]
dir_to_save <- args[3]

burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, b_3) %>%
 mutate(chrom = as.character(chrom))

burden_number <- burden %>%
  select(pos, b_3) %>%
  filter(b_3 != ".") %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) %>%
  group_by(b_3) %>%
  summarize(count = n(), .groups = "drop") %>%  
  distinct(b_3, .keep_all = TRUE) 

#coverage
coverage <- fread(coverage_path, header = FALSE) %>%
  select(V1, V2, V4) %>%
  rename(chrom = V1,
         pos = V2,
         coverage = V4) %>%
  filter(coverage > 0) %>%
   mutate(chrom = as.character(chrom)) %>%
  left_join(burden, by = c("chrom","pos"))  %>%
  filter(b_3 != ".") %>%
  filter(!is.na(b_3)) %>%
  filter(!is.na(coverage)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(coverage, b_3) %>%
  group_by(b_3) %>%
  summarise(total_coverage = sum(coverage), .groups = 'drop') %>%
  left_join(burden_number, by = c("b_3")) %>%
  group_by(b_3) %>%
  mutate(total = total_coverage/count) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) 


safe_quadratic_plateau <- function(coverage, b_3, total) {
  tryCatch(
    quadratic_plateau(coverage, b_3, total, tidy = TRUE),
    error = function(e) {
      tibble(
        STVt = NA_character_
      )
    }
  )
}

table_qp <- safe_quadratic_plateau(coverage, b_3, total)



MD_value <- table_qp %>%
  mutate(
    STVt_num = as.numeric(STVt),  # convert once
    MD_value = case_when(
      is.na(STVt_num) ~ NA_real_,  # keep NA
      STVt_num > 30 ~ 0,           # if >30, return 0
      TRUE ~ ceiling(STVt_num)     # otherwise ceiling
    )
  ) %>%
  select(MD_value)

write.table(MD_value, dir_to_save, 
            col.names = FALSE,
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE,
            na = "NA")

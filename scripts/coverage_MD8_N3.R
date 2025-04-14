#TO HAVE COVERAGE FOR PROBESET TO COMPARE ACCROSS ALL PROBESETS MITO AND NUCLEAR RESULTS#

#TAKE THE LESS PERMISSIVE PROBE => N=3 AND B >= 8

library(knitr)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
coverage_path <- args[1]
burden_path <- args[2]
dir_to_save <- args[3]
n_score <-  args[4]
b_score <- args[5]

indexlibid <- sub(".*/mappedbams/([^/]+)/.*$", "\\1", dir_to_save)

number_of_snps <- fread(burden_path) %>%
  mutate(n_3 = as.factor(n_3))  %>%
  select(pos, n_3, b_3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(n_3 == 3, b_3 >= 8) %>%
  summarize(count = n(), .groups = "drop")

burden <- fread(burden_path) %>%
  mutate(n_3 = as.factor(n_3))  %>%
  select(pos, n_3, b_3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(n_3 == 3, b_3 >= 8) 
#%>%
  #summarize(count = n(), .groups = "drop") %>%  
  #distinct(b_3, n_3, .keep_all = TRUE) 

#coverage
coverage <- fread(coverage_path, header = FALSE) %>%
  select(V1, V2, V3) %>%
  rename(pos = V2,
         coverage = V3) %>%
  filter(coverage > 0) %>%
  full_join(burden) %>%
  filter(!is.na(coverage)) %>%
  filter(!is.na(b_3)) %>%
  summarise(total_coverage = sum(coverage), .groups = 'drop') %>%
  #mutate(indexlibid = indexlibid) %>%
  mutate(total_coverage = as.numeric(total_coverage)) %>%
  mutate(norm_cov = (total_coverage/number_of_snps$count)*1000) %>%
  select(norm_cov)


path_to_save <- paste0(dir_to_save, indexlibid,"_uniqL35MQ25_MD",b_score,"_N", n_score, "_deam.cov_burden.txt")
write.table(coverage, file = path_to_save, row.names = FALSE,  col.names = FALSE)
library(knitr)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
deam_path <- args[1]
summary_path <- args[2]
dir_to_save <- args[3]
score_b <- as.numeric(args[4])
n_score <- as.numeric(args[5])


#get deamination stats

pattern_deam <- "\\.summary_damage.txt$"

file_list_deam <- list.files(deam_path, pattern = pattern_deam, full.names = TRUE, recursive = TRUE)

#create extra column to have burden score information in the final table
all_deam_burden <- rbindlist(lapply(file_list_deam, function(file) {
  data <- fread(file)
  Burden_score <- sub(".*/(\\d+_filter)/.*", "\\1", file)
  data$Burden_score <- Burden_score
  return(data)
}), fill = TRUE)

all_deam <- all_deam_burden %>%
  filter(!is.na(`#ID`)) %>%
  filter(!is.na(`5'CT`)) %>%
  select(-V1) %>%
  mutate(Kraken = ifelse(grepl("Primates", `#ID`), "Primates", "all"),
         N_score = ifelse(grepl("n_0", `#ID`), "ALL", 
                          ifelse(grepl("n_3", `#ID`), "Only n_3", NA))) %>%
  mutate(Library = sub("^([^\\_]+).*", "\\1", `#ID`)) %>%
  select(-`#ID`)


#get count reads summary

pattern_summary <- "\\.pipeline_summary.txt$"

file_list_summaries <- list.files(summary_path, pattern = pattern_summary, full.names = TRUE, recursive = TRUE)


all_summaries_burden <- rbindlist(lapply(file_list_summaries, function(file) {
  data <- fread(file)
  Burden_score <- sub(".*/(\\d+_filter)/.*", "\\1", file)
  data$Burden_score <- Burden_score
  return(data)
}), fill = TRUE)

save_deam <- paste0(dir_to_save,"all_deam.txt")
save_summary <- paste0(dir_to_save,"all_summaries.txt")

write.table(all_summaries,save_summary , sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_deam, save_deam, sep = "\t", row.names = FALSE, quote = FALSE)



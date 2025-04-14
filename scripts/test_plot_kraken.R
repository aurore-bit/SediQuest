library(knitr)
library(data.table)
library(tidyverse)

deam_stats <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/"

pattern <- "\\.summary_damage.txt$"

file_list_deam <- list.files(deam_stats, pattern = pattern, full.names = TRUE, recursive = TRUE)

all_deam <- rbindlist(lapply(file_list_deam, function(file) {
  data <- fread(file)
  Burden_score <- sub(".*/(\\d+_filter)/.*", "\\1", file)
  
  # Add the extracted Burden_score as a new column
  data$Burden_score <- Burden_score
  
  return(data)
}), fill = TRUE)


all_deam <- all_deam %>%
  filter(!is.na(`#ID`)) %>%
  filter(!is.na(`5'CT`)) %>%
  select(-V1) %>%
  mutate(Kraken = ifelse(grepl("Primates", `#ID`), "Primates", "all"),
         N_score = ifelse(grepl("n_0", `#ID`), "ALL", 
                          ifelse(grepl("n_3", `#ID`), "Only n_3", NA))) %>%
  mutate(Library = sub("^([^\\_]+).*", "\\1", `#ID`)) %>%
  select(-`#ID`)

  
summaries <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/"

pattern <- "\\.pipeline_summary.txt$"
file_list_summaries <- list.files(summaries, pattern = pattern, full.names = TRUE, recursive = TRUE)

all_summaries <- rbindlist(lapply(file_list_summaries, fread), fill=TRUE)


write.table(all_summaries_deam, "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/all_summaries_deam.txt", sep = "\t", row.names = FALSE, quote = FALSE)


    

all_summaries <- rbindlist(lapply(file_list_summaries, fread), fill = TRUE) 

# Load the sample data
sample <- fread("/mnt/expressions/Aurore/sediment_pipeline_v0/config/samples.csv")

# Construct the path dynamically using paste0 to concatenate
summaries <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/", sample$indexlibid, "/TW1/rmdupL35MQ25/target/0_filter/deam")

# Define the pattern for kraken_spc files
pattern <- "\\.kraken_spc$"

# List files matching the pattern
file_list_summaries <- list.files(summaries, pattern = pattern, full.names = TRUE, recursive = TRUE)

# Read in all summary files and combine them into a single data.table
all_summaries <- rbindlist(lapply(file_list_summaries, fread), fill = TRUE)

sites_b_0_n_3.filtered.txt
data <- fread("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/Lib.L.5384/TW1/rmdupL35MQ25/target/0_filter/sites_b_0_n_3.filtered.txt")


read_to_burden <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/Lib.L.5239/TW1/rmdupL35MQ25/target/0_filter/deam/Lib.L.5239.read_summary.txt.gz"

kraken_byread <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/Lib.L.5239/TW1/rmdupL35MQ25/target/0_filter/deam/Lib.L.5239.byread"

burden_path <- "/mnt/expressions/benjamin_vernot/soil_capture_2017/site_categories_for_capture/twist_1240k_2023_09/twist.1240k.burden.txt"

kraken_target <- fread(kraken_byread, header =FALSE)

burden_to_read <- fread(read_to_burden) %>%
  select(read_id,chrom,pos)

burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, n_3, b_3) %>%
  filter(n_3==3)


family_colors <- c(
  "Other" = "#FFDAB9",           # Pastel Peach
  "Cetacea" = "#AEC6CF",         # Pastel Blue
  "Rodentia" = "#FFB6A2",        # Pastel Orange
  "Glires" = "#B2E2B0",          # Pastel Green
  "Ruminantia" = "#F4A6A1",      # Pastel Red
  "Primates" = "#B6A7D1",        # Pastel Purple
  "Feliformia" = "#D7A8A0",      # Pastel Brown
  "Afrotheria" = "#F1C3D5",      # Pastel Pink
  "Equidae" = "#A1E1E6",         # Pastel Cyan
  "Chiroptera" = "#E0A6B9",      # Pastel Crimson
  "Caniformia" = "#E2D97F",      # Pastel Olive Green
  "Carnivora" = "#F1C27D"        # Pastel Yellow-orange
)





#####FAMILY LEVEL#######
#non cumulatif
kraken_output <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  select(read_id, level_1, n_3, b_3) %>%
  group_by(b_3, n_3) %>%
  mutate(global_count = n()) %>%
  ungroup() %>%
  group_by(b_3, n_3, level_1) %>%
  mutate(read_count = n()) %>%
  select(-read_id) %>%  
  distinct()  %>%
  group_by(b_3, n_3, level_1) %>%
  mutate(percentage = (read_count/global_count)*100) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3==3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30)

plot_kraken_mapped_bam <- ggplot(kraken_output, aes(x = b_3, y = percentage, fill = level_1)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = species_colors) +  # Manually set colors for each species
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of Each Species at Each Step",
    x = "burden_score",
    y = "Percentage_nc",
    fill = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_kraken_mapped_bam


#cumulatif
#to have total number of reads for each b_3 and n_3 (only n_3=3 here)
kraken_cumulatif <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_1)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3 == 3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_1, n_3, b_3) %>%
  group_by(b_3, n_3) %>%
  mutate(global_count = n()) %>%
  ungroup() %>%
  distinct(b_3, n_3, global_count) %>%
  arrange(b_3) %>%
  mutate(cumulative_global_count = rev(cumsum(rev(global_count)))) %>%
  select(b_3, n_3, cumulative_global_count, global_count)



kraken_read_count <-  kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_1)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3==3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_1, n_3, b_3) %>%
  group_by(b_3, n_3, level_1) %>%
  mutate(read_count = n()) %>%
  ungroup() %>%
  select(-read_id)%>%
  distinct()


new_kraken_output <- expand.grid(
  b_3 = unique(kraken_read_count$b_3),  
  level_1 = unique(kraken_read_count$level_1)  
)

final_output <- new_kraken_output %>%
  left_join(kraken_read_count %>%
              select(b_3, level_1, read_count), by = c("b_3", "level_1")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_1) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count))))%>%
  ungroup() %>%
  full_join(kraken_cumulatif, by = c("b_3"))%>%
  group_by(b_3, n_3, level_1) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup()
  
plot_kraken_mapped_bam <- ggplot(final_output, aes(x = b_3, y = percentage, fill = level_1)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of kraken family assignation after target filtering for different burden score",
    x = "burden_score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90) 
plot_kraken_mapped_bam



#####everyhting  LEVEL#######
# cumulatif
kraken_cumulatif <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3 == 3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3) %>%
  mutate(global_count = n()) %>%
  ungroup() %>%
  distinct(b_3, n_3, global_count) %>%
  arrange(b_3) %>%
  mutate(cumulative_global_count = rev(cumsum(rev(global_count)))) %>%
  select(b_3, n_3, cumulative_global_count, global_count)



kraken_read_count <-  kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3==3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3, level_2) %>%
  mutate(read_count = n()) %>%
  ungroup() %>%
  select(-read_id)%>%
  distinct()


new_kraken_output <- expand.grid(
  b_3 = unique(kraken_read_count$b_3),  
  level_2 = unique(kraken_read_count$level_2)  
)

final_output <- new_kraken_output %>%
  left_join(kraken_read_count %>%
              select(b_3, level_2, read_count), by = c("b_3", "level_2")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_2) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count))))%>%
  ungroup() %>%
  full_join(kraken_cumulatif, by = c("b_3"))%>%
  group_by(b_3, n_3, level_2) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup() 
  #filter(percentage!="0")

plot_kraken_mapped_bam <- ggplot(final_output, aes(x = b_3, y = percentage, fill = level_2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of kraken family assignation after target filtering for different burden score",
    x = "burden_score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    legend.position = "none"  # Remove the legend
  ) + 
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90)       # Rotate text to be vertical

plot_kraken_mapped_bam

#####SPECIES LEVEL#######
# cumulatif
kraken_cumulatif_SPECIES <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  filter(grepl("_", level_2)) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3 == 3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3) %>%
  mutate(global_count = n()) %>%
  ungroup() %>%
  distinct(b_3, n_3, global_count) %>%
  arrange(b_3) %>%
  mutate(cumulative_global_count = rev(cumsum(rev(global_count)))) %>%
  select(b_3, n_3, cumulative_global_count, global_count)



kraken_read_count_SPECIES <-  kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  filter(grepl("_", level_2)) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3==3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3, level_2) %>%
  mutate(read_count = n()) %>%
  ungroup() %>%
  select(-read_id)%>%
  distinct()


new_kraken_output <- expand.grid(
  b_3 = unique(kraken_read_count_SPECIES$b_3),  
  level_2 = unique(kraken_read_count_SPECIES$level_2)  
)

kraken_SPECIES <- new_kraken_output %>%
  left_join(kraken_read_count_SPECIES %>%
              select(b_3, level_2, read_count), by = c("b_3", "level_2")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_2) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count))))%>%
  ungroup() %>%
  full_join(kraken_cumulatif_SPECIES, by = c("b_3"))%>%
  group_by(b_3, n_3, level_2) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup() %>%
 filter(percentage > "1") 
 

plot_kraken_SPECIES <- ggplot(kraken_SPECIES, aes(x = b_3, y = percentage, fill = level_2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of kraken family assignation after target filtering for different burden score",
    x = "burden_score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90)       # Rotate text to be vertical

plot_kraken_SPECIES
  


library(knitr)
library(data.table)
library(tidyverse)

read_to_burden <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/A0101/1240k/rmdupL35MQ25/target/0_filter/A0101.read_summary.txt.gz"

kraken_byread <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/A0101/1240k/rmdupL35MQ25/target/0_filter/A0101.byread"

burden_path <- "/mnt/expressions/Aurore/sediment_pipeline_test/1240k/1240k.burden.txt"
kraken_target <- fread(kraken_byread, header =FALSE)

burden_to_read <- fread(read_to_burden) %>%
  select(read_id,chrom,pos)

burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, n_3, b_3)


kraken_cumulatif_fam_spe <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3 == 3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3) %>%
  mutate(global_count = n()) %>%
  ungroup() %>%
  distinct(b_3, n_3, global_count) %>%
  arrange(b_3) %>%
  mutate(cumulative_global_count = rev(cumsum(rev(global_count)))) %>%
  select(b_3, n_3, cumulative_global_count, global_count)



kraken_read_count_fam_spe <-  kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read, by = "read_id") %>%  
  full_join(burden, by = c("chrom", "pos")) %>%
  filter(b_3 != ".") %>%
  filter(!is.na(read_id)) %>%
  filter(!is.na(level_2)) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  filter(n_3==3) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(read_id, level_2, n_3, b_3) %>%
  group_by(b_3, n_3, level_2) %>%
  mutate(read_count = n()) %>%
  ungroup() %>%
  select(-read_id)%>%
  distinct()


new_kraken_output <- expand.grid(
  b_3 = unique(kraken_read_count_fam_spe$b_3),  
  level_2 = unique(kraken_read_count_fam_spe$level_2)  
)

#regrouping every family and species <1% in "Other" and removing number of reads for this cattegory
kraken_fam_spe <- new_kraken_output %>%
  left_join(kraken_read_count_fam_spe %>%
              select(b_3, level_2, read_count), by = c("b_3", "level_2")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_2) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count)))) %>%
  ungroup() %>%
  full_join(kraken_cumulatif_fam_spe, by = c("b_3")) %>%
  group_by(b_3, n_3, level_2) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup() %>%
  mutate(level_2 = ifelse(percentage < 1, "Other", level_2)) %>% 
  group_by(b_3,level_2) %>%
  mutate(cumulative_count = ifelse(level_2 == "Other", sum(cumulative_count), cumulative_count))%>% 
  mutate(percentage = ifelse(level_2 == "Other", sum(percentage), percentage))%>%
  ungroup() %>%
  distinct(b_3,level_2, cumulative_count, percentage)

#distinct(level_2, b_3, cumulative_count, .keep_all = TRUE) 

plot_kraken_fam_spe <- ggplot(kraken_fam_spe, aes(x = b_3, y = percentage, fill = level_2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of kraken species if available otherwise family assignation after target filtering by burden score (and n=3)",
    x = "burden_score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90,)   # Rotate text to be vertical

plot_kraken_fam_spe





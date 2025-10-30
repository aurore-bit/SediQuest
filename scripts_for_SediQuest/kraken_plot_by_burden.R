#new kraken plot with number of reads with cumulatif info for burden score for both target and deam
#at the order and family/species level only for n score 3

library(knitr)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
burden_path <- args[1]
info_kraken_read_target_path <- args[2]
classification_kraken_byread_target_path <- args[3]
info_kraken_read_deam_path <- args[4]
classification_kraken_byread_deam_path <- args[5]
dir_to_save <- args[6]
n_score <- as.numeric(args[7])

indexlibid <- sub("^([^/]+)/([^/]+)/.*$", "\\2", dir_to_save)

burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, n_3, b_3)


####ORDER####
order_colors <- c(
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

##TARGET##
classification_kraken_byread_target <- fread(classification_kraken_byread_target_path, header =FALSE)

info_kraken_read_target <- fread(info_kraken_read_target_path) %>%
  select(read_id,chrom,pos)

#cumulatif
#to have total number of reads for each b_3 and n_3 (only n_3=3 here)
kraken_cumulatif_target <- classification_kraken_byread_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(info_kraken_read_target, by = "read_id") %>%  
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



kraken_read_count_target <-  classification_kraken_byread_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(info_kraken_read_target, by = "read_id") %>%  
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
  b_3 = unique(kraken_read_count_target$b_3),  
  level_1 = unique(kraken_read_count_target$level_1)  
)

final_output_kraken_target <- new_kraken_output %>%
  left_join(kraken_read_count_target %>%
              select(b_3, level_1, read_count), by = c("b_3", "level_1")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_1) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count))))%>%
  ungroup() %>%
  full_join(kraken_cumulatif_target, by = c("b_3"))%>%
  group_by(b_3, n_3, level_1) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup()

plot_kraken_target <- ggplot(final_output_kraken_target, aes(x = b_3, y = percentage, fill = level_1)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  scale_fill_manual(values = order_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Kraken order assignation in cumulative percentage before target filtering",
    x = "MD score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90) 

#path_to_save <- paste0(dir_to_save,indexlibid, "_kraken_target_byburden.pdf")
#ggsave(path_to_save, plot = plot_kraken_target, device = "pdf")

##DEAM##
classification_kraken_byread_deam <- fread(classification_kraken_byread_deam_path, header =FALSE)

info_kraken_read_deam <- fread(info_kraken_read_deam_path) %>%
  select(read_id,chrom,pos)

#cumulatif
#to have total number of reads for each b_3 and n_3 (only n_3=3 here)
kraken_cumulatif_deam <- classification_kraken_byread_deam%>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(info_kraken_read_deam, by = "read_id") %>%  
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



kraken_read_count_deam <-  classification_kraken_byread_deam %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  full_join(info_kraken_read_deam, by = "read_id") %>%  
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
  b_3 = unique(kraken_read_count_deam$b_3),  
  level_1 = unique(kraken_read_count_deam$level_1)  
)

final_output_deam <- new_kraken_output %>%
  left_join(kraken_read_count_deam %>%
              select(b_3, level_1, read_count), by = c("b_3", "level_1")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_1) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count))))%>%
  ungroup() %>%
  full_join(kraken_cumulatif_deam, by = c("b_3"))%>%
  group_by(b_3, n_3, level_1) %>%
  mutate(percentage = (cumulative_count / sum(cumulative_global_count)) * 100) %>%
  filter(b_3 < 30) %>%
  ungroup()

plot_kraken_deam <- ggplot(final_output_deam, aes(x = b_3, y = percentage, fill = level_1)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  scale_fill_manual(values = order_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Kraken order assignation in cumulative percentage after deam filtering",
    x = "MD score",
    y = "Percentage",
    fill = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90) 


plot_all_kraken_burden <- grid.arrange(plot_kraken_target, plot_kraken_deam, ncol = 2)  

path_to_save <- paste0(dir_to_save,indexlibid,"_kraken_order_byburden.pdf")
ggsave(path_to_save, plot = plot_all_kraken_burden, device = "pdf", width = 14, height = 7, 
       units = "in", dpi = 300)


####FAMILY/SPECIES####

#TARGET#
kraken_target <- fread(classification_kraken_byread_target_path, header =FALSE)

burden_to_read_target <- fread(info_kraken_read_target_path) %>%
  select(read_id,chrom,pos)

kraken_cumulatif_fam_spe_target <- kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read_target, by = "read_id") %>%  
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



kraken_read_count_fam_spe_target <-  kraken_target %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read_target, by = "read_id") %>%  
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


new_kraken_output_target <- expand.grid(
  b_3 = unique(kraken_read_count_fam_spe_target$b_3),  
  level_2 = unique(kraken_read_count_fam_spe_target$level_2)  
)

#regrouping every family and species <1% in "Other" and removing number of reads for this cattegory
kraken_fam_spe_target <- new_kraken_output_target %>%
  left_join(kraken_read_count_fam_spe_target %>%
              select(b_3, level_2, read_count), by = c("b_3", "level_2")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_2) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count)))) %>%
  ungroup() %>%
  full_join(kraken_cumulatif_fam_spe_target, by = c("b_3")) %>%
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

plot_kraken_fam_spe_target <- ggplot(kraken_fam_spe_target, aes(x = b_3, y = percentage, fill = level_2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Kraken species/family assignation before target filtering",
    x = "MD score",
    y = "Percentage",
    fill = "Species/Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90,)   # Rotate text to be vertical



#DEAM#

kraken_deam <- fread(classification_kraken_byread_deam_path, header =FALSE)

burden_to_read_deam <- fread(info_kraken_read_deam_path) %>%
  select(read_id,chrom,pos)

kraken_cumulatif_fam_spe_deam <- kraken_deam %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read_deam, by = "read_id") %>%  
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



kraken_read_count_fam_spe_deam <-  kraken_deam %>%
  rename(read_id = V1,
         level_1 = V2,
         level_2 = V3) %>%
  #filter(grepl("_", level_2)) %>%
  full_join(burden_to_read_deam, by = "read_id") %>%  
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


new_kraken_output_deam <- expand.grid(
  b_3 = unique(kraken_read_count_fam_spe_deam$b_3),  
  level_2 = unique(kraken_read_count_fam_spe_deam$level_2)  
)

#regrouping every family and species <1% in "Other" and removing number of reads for this cattegory
kraken_fam_spe_deam <- new_kraken_output_deam %>%
  left_join(kraken_read_count_fam_spe_deam %>%
              select(b_3, level_2, read_count), by = c("b_3", "level_2")) %>%
  mutate(read_count = ifelse(is.na(read_count), 0, read_count)) %>%
  arrange(b_3) %>%
  group_by(level_2) %>%
  mutate(cumulative_count = rev(cumsum(rev(read_count)))) %>%
  ungroup() %>%
  full_join(kraken_cumulatif_fam_spe_deam, by = c("b_3")) %>%
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

plot_kraken_fam_spe_deam <- ggplot(kraken_fam_spe_deam, aes(x = b_3, y = percentage, fill = level_2)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  #scale_fill_manual(values = family_colors) +  # Apply the custom colors
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Kraken species/family assignation after deam filtering",
    x = "MD score",
    y = "Percentage",
    fill = "Species/Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  geom_text(aes(label = cumulative_count),  # Add the number of reads as labels
            position = position_stack(vjust = 0.5),  # Place the text in the center of each segment
            color = "white",  # Text color
            size = 3,         # Adjust text size
            angle = 90,)   # Rotate text to be vertical


plot_all_kraken_burden_fam_spe <- grid.arrange(plot_kraken_fam_spe_target, plot_kraken_fam_spe_deam, ncol = 2)  

path_to_save <- paste0(dir_to_save,indexlibid,"_kraken_fam_spe_byburden.pdf")
ggsave(path_to_save, plot = plot_all_kraken_burden_fam_spe, device = "pdf", width = 14, height = 7, 
       units = "in", dpi = 300)
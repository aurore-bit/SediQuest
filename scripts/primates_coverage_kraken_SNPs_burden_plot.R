#script to create all summary plots

library(data.table)
library(tidyverse)
library(ggplot2)

#knitr:opts_knit$set(root.dir = normalizePath(".."))

data <- fread("/mnt/expressions/Aurore/sediment_pipeline_v0/config/samples.csv")

burden_path <- paste0("/mnt/expressions/Aurore/sediment_pipeline_test/",
                      data$probeset, "/",
                      data$probeset, ".burden.txt")

coverage_path <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/",
                        data$indexlibid, "/",
                        data$probeset, "/target/rmdupL35MQ25/",
                        data$indexlibid, ".cov")

# Read the burden data using the dynamically created file path
burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, n_3, b_3)

##############
#Kraken plot
##############
#for bam
mapped_bam_path <- "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/"

# List all .kraken_phylo files in the directory and subdirectories
file_list_mapped_bam <- list.files(mapped_bam_path, pattern = "\\.kraken_spc$", full.names = TRUE, recursive = TRUE)

kraken_mapped_bam <- rbindlist(lapply(file_list_mapped_bam, fread)) %>%
  mutate(V2 = as.numeric(V2)) %>%
  group_by(V1) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = V1, values_from = V2, values_fn = list) %>%
  ungroup() %>%
  select(-row_id) %>%
  mutate(across(where(is.list), ~ unlist(.))) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(step = rep(c("split", "mapped", "target", "unique&qual&length", "deam"), length.out = n())) %>%
  pivot_longer(cols = -c(step, total), names_to = "Species", values_to = "Value") %>%
  mutate(Percentage = (Value / total) * 100) %>%
  select(Percentage, Species, step)

species_colors <- c(
  "Primates" = "lightblue",
  "Glires" = "green",
  "Laurasiatheria" = "red",
  "Afrotheria" = "purple",
  "Mammalia" = "orange",
  "Bacteria" = "brown",
  "Sauropsida" = "yellow",
  "root" = "pink",
  "unclassified" = "gray"
)


plot_kraken_mapped_bam <- ggplot(kraken_mapped_bam, aes(x = step, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack bars for each step
  scale_fill_manual(values = species_colors) +  # Manually set colors for each species
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "Percentage of Each Species at Each Step",
    x = "Step",
    y = "Percentage",
    fill = "Species"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

path_to_save <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/", data$indexlibid, "/", data$probeset, "/",data$indexlibid, "_kraken_mapped.pdf")
ggsave(path_to_save, plot = plot_kraken_mapped_bam, device = "pdf")


##############
#Coverage plot
##############
burden_number <- burden %>%
  select(pos, n_3, b_3) %>%
  filter(b_3 != ".") %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) %>%
  group_by(n_3, b_3) %>%
  summarize(count = n(), .groups = "drop") %>%  
  distinct(b_3, n_3, .keep_all = TRUE)

coverage <-  fread(coverage_path, header = FALSE) %>%
  select(V1, V2, V4) %>%
  rename(pos = V2,
         coverage = V4) %>%
  filter(coverage > 0) %>%
  left_join(burden, by = "pos")  %>%
  filter(b_3!=".") %>%
  filter(!is.na(b_3)) %>%
  filter(!is.na(coverage)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(coverage, n_3,b_3) %>%
  group_by(b_3, n_3) %>%
  summarise(total_coverage = sum(coverage), .groups = 'drop') %>%
  left_join(burden_number, by = c("b_3", "n_3")) %>%
  group_by(b_3, n_3) %>%
  mutate(total_coverage = total_coverage/count) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) 
  
cov_plot <- ggplot(coverage, aes(x = b_3, y = total_coverage, color = as.factor(n_3), group = n_3)) + 
  geom_line(linetype = "dashed") +  # Dashed lines
  geom_point() +                    # Points at each data point
  labs(x = "B Score (b_3)", 
       y = "Coverage", 
       color = "N Score") + 
  ylim(0,1)+
  theme_minimal()

path_to_save <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/", data$indexlibid, "/", data$probeset, "/",data$indexlibid, "_burden_primates_plot.pdf")
ggsave(path_to_save, plot = coverage, device = "pdf")


##############
#Primates % plot
##############

#kraken info
path_to_kraken <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/", data$indexlibid, "/", data$probeset, "/target/rmdupL35MQ25/", data$indexlibid, ".byread")
kraken <-  fread(path_to_kraken, header=F) %>%
  select(V1,V2) %>%
  rename(read=V1)

#take positions and reads names assuming . value are 0 for b score
#to get the number of read for each burden scores
path_to_primate <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/", data$indexlibid, "/", data$probeset, "/target/rmdupL35MQ25/",data$indexlibid, ".bed")
primates <- fread(path_to_primate) %>%
  select(V1, V2, V3, V4) %>%
  rename(chrom = V1,
         pos0 = V2,
         pos = V3,
         read = V4) %>%
  full_join(burden, by = c("chrom", "pos0", "pos")) %>%
  filter(b_3 != ".") %>%
  full_join(kraken, by = "read") %>%
  filter(!is.na(read)) %>%
  group_by (b_3, n_3) %>%
  mutate(read_count = n()) %>%
  group_by(b_3,n_3) %>%
  mutate(primate_count = sum(V2 == "Primates", na.rm = TRUE)) %>%
  mutate(percentage_primate = (primate_count/read_count)*100) %>%
  distinct(b_3, n_3, .keep_all = TRUE) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30)

#burden plot
burden_primates_plot <- ggplot(primates, aes(x = b_3, y = percentage_primate, color = as.factor(n_3), group = n_3)) + 
  geom_line() +
  labs(x = "B Score (b_3)", 
       y = "Percentage of Primates", 
       color = "N Score") + 
  theme_minimal()

#save plot
path_to_save <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/", data$indexlibid, "/", data$probeset, "/", data$indexlibid, "_burden_primates_plot.pdf")
ggsave(path_to_save, plot = burden_primates_plot, device = "pdf")


##############
#SNPs
##############

#number of SNPs for each b and n score
burden_SNPs <- burden %>%
  select(pos, n_3, b_3) %>%
  filter(b_3 != ".") %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) %>%
  group_by(n_3, b_3) %>%
  summarize(count = n(), .groups = "drop") %>%  # Use summarize instead of mutate for count
  distinct(b_3, n_3, .keep_all = TRUE)


burden_SNPs_plot <- ggplot(burden_number, aes(x = b_3, y = count, fill = as.factor(n_3))) +
  geom_bar(stat = "identity") +
  labs(x = "B Score (b)", 
       y = "Total Number of Positions", 
       fill = "N Score") +  # Correctly labels the fill aesthetic
  theme_minimal()

#save plot
path_to_save <- paste0("/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0_summary/" , data$indexlibid, "/",data$probeset, "/", data$indexlibid, "_burden_SNPs_plot.pdf")
ggsave(path_to_save, plot = plot, device = "pdf")

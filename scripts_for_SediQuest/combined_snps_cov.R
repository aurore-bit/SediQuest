
#TO HAVE COVERAGE ACCROSS BURDEN SCORES#

library(knitr)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
coverage_path <- args[1]
burden_path <- args[2]
dir_to_save <- args[3]

indexlibid <- sub("^([^/]+)/([^/]+)/.*$", "\\2", dir_to_save)


burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, b_3)


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
  rename(pos = V2,
         coverage = V4) %>%
  filter(coverage > 0) %>%
  left_join(burden, by = "pos")  %>%
  filter(b_3 != ".") %>%
  filter(!is.na(b_3)) %>%
  filter(!is.na(coverage)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  select(coverage, b_3) %>%
  group_by(b_3) %>%
  summarise(total_coverage = sum(coverage), .groups = 'drop') %>%
  left_join(burden_number, by = c("b_3")) %>%
  group_by(b_3) %>%
  mutate(total_coverage = total_coverage/count) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) 



data <- coverage %>%
  #left_join(burden_number, by = c("b_3", "n_3", "count"))
  left_join(burden_number, by = c("b_3"))


# plot coverage and number of SNPs
comb_cov_snps <- ggplot(data) +
  geom_line(aes(x = b_3, y = total_coverage)) +
  labs(x = "Mammalian diversity Score") +
  theme(
 #   axis.title.y = element_text(color = "blue", size = 13),  # Left axis title (Coverage)
  #  axis.title.y.right = element_text(color = "red", size = 13),  # Right axis title (burden_SNPs)
    axis.text = element_text(size = 12),  # Axis text size
    axis.title.x = element_text(size = 13)  # X-axis title size
  ) +
  labs(
    x = "MD score",
    y = "Coverage",
   )+
 ggtitle("Coverage by Mammalian diversity score")  +
  coord_cartesian(xlim = c(0,30),ylim = c(0,1))


path_to_save <- paste0(dir_to_save, indexlibid,"_cov_MD.pdf")
ggsave(path_to_save, plot = comb_cov_snps, device = "pdf")
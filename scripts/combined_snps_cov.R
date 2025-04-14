
#TO HAVE SNPS NUMBER AND COVERAGE IN THE SAME GRAPH#

library(knitr)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
coverage_path <- args[1]
burden_path <- args[2]
dir_to_save <- args[3]
n_score <- as.numeric(args[4])

indexlibid <- sub("^([^/]+)/([^/]+)/.*$", "\\2", dir_to_save)

burden <- fread(burden_path) %>%
  select(chrom, pos0, pos, n_3, b_3)


burden_number <- burden %>%
  mutate(n_3 = as.factor(n_3))  %>%
  select(pos, n_3, b_3) %>%
  filter(b_3 != ".") %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) %>%
  group_by(n_3, b_3) %>%
  summarize(count = n(), .groups = "drop") %>%  
  distinct(b_3, n_3, .keep_all = TRUE) 

#coverage
coverage <- fread(coverage_path, header = FALSE) %>%
  select(V1, V2, V4) %>%
  rename(pos = V2,
         coverage = V4) %>%
  filter(coverage > 0) %>%
  left_join(burden, by = "pos")  %>%
  filter(b_3!=".") %>%
  filter(!is.na(b_3)) %>%
  filter(!is.na(coverage)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  mutate(n_3 = as.factor(n_3))  %>%
  select(coverage, n_3,b_3) %>%
  group_by(b_3, n_3) %>%
  summarise(total_coverage = sum(coverage), .groups = 'drop') %>%
  left_join(burden_number, by = c("b_3", "n_3")) %>%
  group_by(b_3, n_3) %>%
  mutate(total_coverage = total_coverage/count) %>%
  mutate(n_3 = as.factor(n_3)) %>%
  mutate(b_3 = as.numeric(as.character(b_3))) %>%
  filter(b_3 < 30) 


#number of position max (will change from one capture set to another)
max_burden <- max(burden_number$count, na.rm = TRUE)

data <- coverage %>%
  left_join(burden_number, by = c("b_3", "n_3", "count"))

coeff <- max_burden 

# plot coverage and number of SNPs
comb_cov_snps <- ggplot(data) +
  geom_line(aes(x = b_3, y = total_coverage, color = as.factor(n_3), linetype = as.factor(n_3))) +
  #geom_line(aes(x = b_3, y = count / coeff, color = as.factor(n_3),linetype = as.factor(n_3)) )+ 
  
  scale_y_continuous(
    name = "Coverage",  
   # limits = c(0, 1),   
    sec.axis = sec_axis(~ . * coeff, name = "Number of positions")  
  ) +
  scale_color_discrete(name = "N Score") +
  labs(x = "B Score") +
  theme(
 #   axis.title.y = element_text(color = "blue", size = 13),  # Left axis title (Coverage)
  #  axis.title.y.right = element_text(color = "red", size = 13),  # Right axis title (burden_SNPs)
    axis.text = element_text(size = 12),  # Axis text size
    axis.title.x = element_text(size = 13)  # X-axis title size
  ) +
  scale_linetype_manual(values = c("1" = "dashed", "2" = "dashed", "3" = "solid"))  + 
  ggtitle("Coverage by burden score and n score")  +
  coord_cartesian(xlim = c(0,30),ylim = c(0,1))


path_to_save <- paste0(dir_to_save, indexlibid,"_SNPs_number_cov_combined.pdf")
ggsave(path_to_save, plot = comb_cov_snps, device = "pdf")
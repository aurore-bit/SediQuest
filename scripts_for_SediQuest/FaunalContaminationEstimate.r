#Calculate faunal estimate based on kraken simulations

library(knitr)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(soiltestcorr)


args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

#Number of Primates reads in my sample

kraken_primates <- fread(input, header = FALSE) %>%
  separate(
    V2,
    into = paste0("level_", 1:27),
    sep = ";",
    fill = "right"
  ) %>%
  filter(!is.na(level_24)) %>%
  mutate(
    order = case_when(
      grepl("Primates", level_24) ~ "Primates",
      TRUE ~ "Other"
    )) %>%
  group_by(order) %>%
  mutate(global_count = n()) %>%
  distinct(order, global_count) 
 # full_join(coverage_no_damage, by = "V1") %>%
 # group_by(V2, order) %>%
 # mutate(order_count = n()) %>%
  #distinct(V2, order_count, order, global_count) %>%
 # mutate(V2 = as.numeric(V2))
  

#using the flying lemur no damage kraken output values as 100 contamination and Neandertal as 0
contamination <- kraken_primates %>%
  ungroup() %>%
  tidyr::complete(order = c("Primates", "Other"),
    fill = list(global_count = 0)) %>%
  pivot_wider(
    names_from = order,
    values_from = global_count,
     values_fill = list(global_count = 0)
  ) %>%
  mutate(
    ratio = Primates / (Primates + Other)
  ) %>%
  mutate(
  contam = case_when(
    ratio == 1 ~ 0,
    ratio == 0 ~ 1,
    TRUE ~ (0.9981 - ratio) / (0.9981 - 0.09650)
  )
) %>%
  select(contam)


write.table(contamination, output, 
            col.names = FALSE,
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

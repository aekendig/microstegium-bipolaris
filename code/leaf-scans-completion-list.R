##### info ####

# file: leaf-scans-completion-list-2019
# author: Amy Kendig
# date last edited: 6/30/19
# goal: check that leaf scans downloaded from One Drive (initial storage space) are complete


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# retrieve list of files
jun19 <- tibble(
  files = list.files("../Leaf Scans/Leaf Scans June 2019")
)

# data collection sheet
djn <- read_csv("./data/focal-disease-jun-2019-density-exp.csv")


#### edit data ####

# separate experiments
jn19de <- jun19 %>%
  filter(substr(files, 1, 1) == "D")

jn19le <- jun19 %>%
  filter(substr(files, 1, 1) == "L")

# new columns
jn19de <- jn19de %>%
  mutate(
    site = substr(files, 1, 2),
    plot = substr(files, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.numeric(),
    treatment = substr(files, 5, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "W" = "water"),
    sp = ifelse(grepl("Ev", files, fixed = T), "Ev", "Mv"),
    ID = case_when(sp == "Ev" ~ substr(files, 10, 11) %>% gsub(c("_"), "", .),
                   sp == "Mv" ~ ifelse(grepl("Edge", files, fixed = T), 
                                       "Edge", 
                                       substr(files, 10, 11) %>% gsub(c("_"), "", .))),
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
  ) 

# check that the above worked as expected
unique(jn19de$site)
unique(jn19de$plot) # one is missing plot number
unique(jn19de$treatment)
unique(jn19de$sp)
unique(jn19de$ID)
filter(jn19de, !(ID %in% c("1", "2", "3", "A", "Edge")))
unique(jn19de$focal)


#### analyze data ####

# extract plots
plots <- djn %>%
  select(site, plot, treatment) %>%
  unique()

# compare plots to scans
plots %>% anti_join(jn19de)

# check that all edge scans are included
plots %>% anti_join(filter(jn19de, ID == "Edge"))

# see why scans weren't collected
djn %>%
  filter(scan == 0) %>%
  select(field_notes) %>%
  unique()

# check that all individual scans are included
djn %>% 
  filter(scan == 1) %>%
  anti_join(filter(jn19de, ID != "Edge"))

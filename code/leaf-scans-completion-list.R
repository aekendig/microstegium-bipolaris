##### info ####

# file: leaf-scans-completion-list-2019
# author: Amy Kendig
# date last edited: 8/7/19
# goal: check that leaf scans downloaded from One Drive (initial storage space) are complete


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# retrieve list of scans
scans_jn <- tibble(
  files = list.files("../Leaf Scans/Leaf Scans June 2019")
)
scans_jl <- tibble(
  files = list.files("..../Downloads/Leaf Scans June 2019")
)

# data collection sheet
field_jn <- read_csv("./data/focal-disease-jun-2019-density-exp.csv")


#### edit data ####

## functions for processing scan data by experiment

# density experiment
dens_fun <- function(dat) {
  dat_out <- dat %>%
    filter(substr(files, 1, 1) == "D") %>%
    mutate(
      site = substr(files, 1, 2),
      plot = substr(files, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.numeric(),
      treatment = substr(files, 5, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "W" = "water"),
      sp = ifelse(grepl("Ev", files, fixed = T), "Ev", "Mv"),
      ID = case_when(sp == "Ev" ~ substr(files, 10, 11) %>% gsub(c("_"), "", .),
                     sp == "Mv" ~ ifelse(grepl("Edge", files, fixed = T) | grepl("edge", files, fixed = T) , 
                                         "Edge", 
                                         substr(files, 10, 11) %>% gsub(c("_"), "", .))),
      focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
    ) 
  
  return(dat_out)
}

# litter experiment
litt_fun <- function(dat) {
  dat_out <- dat %>%
    filter(substr(files, 1, 1) == "L") %>%
    mutate(
      site = substr(files, 8, 9),
      plot = substr(files, 11, 12),
      treatment = substr(plot, 1, 1) %>% 
        recode("R" = "removal", "A" = "addition", "C" = "control"),
      replicate = substr(plot, 2, 2)
      ) 
  
  return(dat_out)
}

# create experiment-specific datasets
scans_jn_d <- dens_fun(scans_jn)
scans_jn_l <- litt_fun(scans_jn)


#### analyze data ####

# manually check litter data because it's small
data.frame(scans_jn_l)

## functions to check that all scans are present, there are no unexpected scnas, and there are no duplicates

# field notes
notes_fun <- function(field_dat) {
  scan_notes <- field_dat %>%
    filter(scan == 0) %>%
    select(site, plot, treatment, sp, ID, field_notes) %>%
    data.frame()
  
  print(scan_notes)
}

# duplicate scans
dup_fun <- function(scan_dat) {
  dup_scans <- scan_dat %>%
    group_by(site, plot, treatment, sp, ID) %>%
    summarise(reps = n()) %>%
    filter(reps > 1) %>%
    data.frame()

  ifelse(nrow(dup_scans) > 0, print(dup_scans), "no duplicate scans")
}

# edge Mv
edge_fun <- function(field_dat, scan_dat) {
  
  plots_missing <- field_dat %>%
    select(site, plot, treatment) %>%
    unique() %>%
    anti_join(filter(scan_dat, focal == 0)) %>%
    data.frame()
  
  ifelse(nrow(plots_missing) > 0, print(plots_missing), print("all plots are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(focal == 0) %>%
    anti_join(select(field_dat, site, plot, treatment) %>% unique()) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, print(unexpected_scans), "no unexpected scans")
}

# focal Mv
mv_fun <- function(field_dat, scan_dat) {
  
  leaves_missing <- field_dat %>%
    filter(sp == "Mv" & scan == 1) %>%
    select(site, plot, treatment, ID) %>%
    anti_join(filter(scan_dat, sp == "Mv" & focal == 1)) %>%
    data.frame()
  
  ifelse(nrow(leaves_missing) > 0, print(leaves_missing), print("all leaves are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(sp == "Mv" & focal == 1) %>%
    anti_join(filter(field_dat, sp == "Mv" & scan == 1) %>% select(site, plot, treatment, ID)) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, print(unexpected_scans), "no unexpected scans")
}

# focal Ev
ev_fun <- function(field_dat, scan_dat) {
  
  leaves_missing <- field_dat %>%
    filter(sp == "Ev" & scan == 1) %>%
    select(site, plot, treatment, ID) %>%
    anti_join(filter(scan_dat, sp == "Ev")) %>%
    data.frame()
  
  ifelse(nrow(leaves_missing) > 0, print(leaves_missing), print("all leaves are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(sp == "Ev") %>%
    anti_join(filter(field_dat, sp == "Ev" & scan == 1) %>% select(site, plot, treatment, ID)) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, print(unexpected_scans), "no unexpected scans")
}


# June data
notes_fun(field_jn)
dup_fun(scans_jn_d)
edge_fun(field_jn, scans_jn_d)
mv_fun(field_jn, scans_jn_d)
ev_fun(field_jn, scans_jn_d)

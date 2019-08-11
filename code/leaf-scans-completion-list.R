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
  files = list.files("~/OneDrive - University of Florida/Flory Lab/Leaf Scans/July 2019")
)

# data collection sheet
field_jn <- read_csv("./data/focal-disease-jun-2019-density-exp.csv")
field_jl <- read_csv("./data/focal-disease-jul-2019-density-exp.csv")


#### edit data ####

## functions for processing scan data by experiment

# density experiment
dens_fun <- function(dat) {
  dat_out <- dat %>%
    filter(grepl("Litter", files, fixed = T) == F) %>%
    mutate(
      site = substr(files, 1, 2),
      plot = substr(files, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.numeric(),
      treatment = substr(files, 5, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "f" = "fungicide", "W" = "water", "w" = "water"),
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
    filter(grepl("Litter", files, fixed = T) == T) %>%
    mutate(
      site = ifelse(substr(files, 1, 1) == "L", substr(files, 8, 9), substr(files, 1, 2)),
      plot = substr(files, 11, 12),
      treatment = substr(plot, 1, 1) %>% 
        recode("R" = "removal", "A" = "addition", "C" = "control"),
      block = substr(plot, 2, 2),
      replicate = str_sub(files, -13) %>% gsub("[^[:digit:]]", "", .) %>% as.numeric()
      ) 
  
  return(dat_out)
}

# create experiment-specific datasets
scans_jn_d <- dens_fun(scans_jn)
scans_jn_l <- litt_fun(scans_jn)

scans_jl_d <- dens_fun(scans_jl)
scans_jl_l <- litt_fun(scans_jl)


#### analyze data ####

# manually check June litter data because it's small
data.frame(scans_jn_l)

# check for replication and duplication in July litter dataset
head(scans_jl_l)
scans_jl_l %>%
  group_by(site, treatment, block) %>%
  summarise(reps = length(unique(replicate)), files = n())

## functions to check that all scans are present, there are no unexpected scans, and there are no duplicates

# field notes
notes_fun <- function(field_dat) {
  scan_notes <- field_dat %>%
    filter(scan == 0) %>%
    select(site, plot, treatment, sp, ID, field_notes) %>%
    data.frame()
  
  print(list("leaves not collected for scans:", scan_notes))
}

# duplicate scans
dup_fun <- function(scan_dat) {
  dup_scans <- scan_dat %>%
    group_by(site, plot, treatment, sp, ID) %>%
    summarise(reps = n()) %>%
    filter(reps > 1) %>%
    data.frame()

  ifelse(nrow(dup_scans) > 0, print(list("duplicate scans:", dup_scans)), "no duplicate scans")
}

# edge Mv
edge_fun <- function(field_dat, scan_dat) {
  
  plots_missing <- field_dat %>%
    select(site, plot, treatment) %>%
    unique() %>%
    anti_join(filter(scan_dat, focal == 0), by = c("site", "plot", "treatment")) %>%
    data.frame()
  
  ifelse(nrow(plots_missing) > 0, print(list("plots missing:", plots_missing)), print("all plots are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(focal == 0) %>%
    anti_join(select(field_dat, site, plot, treatment) %>% unique(), by = c("site", "plot", "treatment")) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, return(list("unexpected scans:", unexpected_scans)), "no unexpected scans")
}

# focal Mv
mv_fun <- function(field_dat, scan_dat) {
  
  leaves_missing <- field_dat %>%
    filter(sp == "Mv" & scan == 1) %>%
    select(site, plot, treatment, ID) %>%
    anti_join(filter(scan_dat, sp == "Mv" & focal == 1), by = c("site", "plot", "treatment", "ID")) %>%
    data.frame()
  
  ifelse(nrow(leaves_missing) > 0, print(list("leaves missing:", leaves_missing)), print("all leaves are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(sp == "Mv" & focal == 1) %>%
    anti_join(filter(field_dat, sp == "Mv" & scan == 1) %>% select(site, plot, treatment, ID), by = c("site", "plot", "treatment", "ID")) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, return(list("unexpected scans:", unexpected_scans)), "no unexpected scans")
}

# focal Ev
ev_fun <- function(field_dat, scan_dat) {
  
  leaves_missing <- field_dat %>%
    filter(sp == "Ev" & scan == 1) %>%
    select(site, plot, treatment, ID) %>%
    anti_join(filter(scan_dat, sp == "Ev"), by = c("site", "plot", "treatment", "ID")) %>%
    data.frame()
  
  ifelse(nrow(leaves_missing) > 0, print(list("leaves missing:", leaves_missing)), print("all leaves are accounted for"))
  
  unexpected_scans <- scan_dat %>%
    filter(sp == "Ev") %>%
    anti_join(filter(field_dat, sp == "Ev" & scan == 1) %>% select(site, plot, treatment, ID), by = c("site", "plot", "treatment", "ID")) %>%
    data.frame()
  
  ifelse(nrow(unexpected_scans) > 0, return(list("unexpected scans:", unexpected_scans)), "no unexpected scans")
}


# June data
notes_fun(field_jn)
dup_fun(scans_jn_d)
edge_fun(field_jn, scans_jn_d)
mv_fun(field_jn, scans_jn_d)
ev_fun(field_jn, scans_jn_d)

# July data
notes_fun(field_jl)
dup_fun(scans_jl_d)
edge_fun(field_jl, scans_jl_d)
mv_fun(field_jl, scans_jl_d)

ev_fun(field_jl, scans_jl_d)
# D1 6F EvA - bag empty
# D1 7F Ev1 - bag empty
# D2 6W is unedited, but all leaves are present, double check with Laney about cropping order (she uploaded it)
# D4 1W Ev3 - bag empty
# D4 1W EvA - leaf included, marked as not collected, check early Aug datasheet
# D4 4W Ev1 - bag empty
# D4 4W EvA - leaf included, marked as not collected, check early Aug datasheet
# Change D3 8F Ev3 on field datasheet to scan collected?

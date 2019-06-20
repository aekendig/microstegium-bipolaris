##### info ####

# file: leaf-scans-data-processing
# author: Amy Kendig
# date last edited: 6/12/19
# goal: combine raw 2019 leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: LeafScanAnalysis_mv_cw_011719.ijm


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import all raw data files
may <- read_csv("./data/leaf-scans-may-2019-density-exp.csv")

# import data files to check collection
ev_may <- read_csv("./data/ev-disease-may-2019-density-exp.csv")


#### edit data ####

# indicate whether images were edited
# remove edited from name
# remove month from name
may <- may %>%
  mutate(edited = ifelse(grepl("edited", slice, fixed = T), 1, 0),
         slice = gsub("_edited", "", slice) %>% gsub("_May", "", .) %>% gsub("_may", "", .),
         month = "May")
nrow(may) # 834
sum(may$edited) # 153

# create list of edited scans
may_ed <- may %>%
  filter(edited == 1) %>%
  select(slice)

# check that all are contained in unedited list
may_ed %>% filter(!(slice %in% filter(may, edited == 0)$slice))
# one sample is not - looked through files and can't figure out why, but it's fine because we have the edited one

# remove unedited scans
may <- may %>%
  filter(edited == 1 | !(slice %in% may_ed$slice))
nrow(may) # 684
834 - 153 + 3

# combine
dat <- may

# new columns
dat <- dat %>%
  mutate(
    plant = gsub(".*:","",slice),
    part = gsub(":.*$", "", slice) %>% recode(lesions = "lesion", greens = "green"),
    site = substr(plant, 1, 2),
    plot = substr(plant, 4, 5) %>% gsub("[^[:digit:]]", "", .) %>% as.numeric(),
    treatment = substr(plant, 5, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "W" = "water"),
    sp = ifelse(grepl("Ev", slice, fixed = T) | grepl("EV", slice, fixed = T), "Ev", "Mv"),
    ID = case_when(sp == "Ev" ~ substr(plant, 10, 12) %>% gsub(c("_"), "", .),
                   sp == "Mv" ~ ifelse(grepl("Edge", slice, fixed = T), 
                                       "Edge", 
                                       substr(plant, 10, 12) %>% gsub(c("_"), "", .))),
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
    age = case_when(sp == "Ev" ~ ifelse(ID == "A" | (focal == 0 & plot > 7), "adult", "seedling"),
                    sp == "Mv" ~ "seedling")
  ) 

# check that the above worked as expected
unique(dat$plant)
unique(dat$part)
unique(dat$site)
unique(dat$plot)
unique(dat$treatment)
unique(dat$sp)
unique(dat$ID)
unique(dat$focal)
unique(dat$age)
head(data.frame(dat), n = 10)
tail(data.frame(dat), n = 10)

# spread by part
datw <- dat %>%
  select(-c(slice)) %>%
  gather(variable, value, -(edited:age)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value)

# indicate values to remove
datw <- datw %>%
  mutate(remove = ifelse(plant %in% c("D3_3F_Ev_A", "D3_6F_Ev_A", "D3_10W_Ev_R8") & month == "May",1, 0))


#### check values ####

# leaf area
datw %>%
  filter(remove == 0) %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(sp == "Ev" & leaf_area.pix > 7.5e5) %>% data.frame()

# percent lesion area 
datw %>%
  filter(remove == 0) %>%
  ggplot(aes(x = (lesion_area.pix - green_area.pix)/leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~sp, scales = "free")

# manually check extremes
datw %>%
  filter(remove == 0) %>%
  filter(sp == "Ev" & (lesion_area.pix - green_area.pix)/leaf_area.pix > 0.6) %>% data.frame()

# check that all plots were collect for Mv
datw %>%
  filter(ID == "Edge") %>%
  group_by(month, site, treatment) %>%
  summarise(n = length(unique(plant)))

# check that all Ev were collected
ev_may %>%
  filter(scan == 1) %>%
  left_join(datw) %>%
  filter(is.na(plant)) %>%
  data.frame()
# missing scans - will ask Laney about these

datw %>%
  filter(sp == "Ev") %>%
  left_join(ev_may) %>%
  filter(is.na(leaves_tot)) %>% 
  data.frame()
# shouldn't have scans - could be the ones above, but mislabelled


#### save data ####

leaf <- datw
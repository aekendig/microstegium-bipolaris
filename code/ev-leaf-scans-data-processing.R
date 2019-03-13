##### info ####

# file: ev-leaf-scans-data-processing
# author: Amy Kendig
# date last edited: 3/12/19
# goal: combine raw 2018 Elymus leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: LeafScanAnalysis_ev_ak_021819.ijm


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd("./data")

# import all raw data files
jul <- read_csv("ev-leaf-scans-jul-2018-density-exp.csv")
jul_ed <- read_csv("edited-ev-leaf-scans-jul-2018-density-exp.csv")
aug <- read_csv("ev-leaf-scans-aug-2018-density-exp.csv")
aug_ed <- read_csv("edited-ev-leaf-scans-aug-2018-density-exp.csv")
sep <- read_csv("ev-leaf-scans-sep-2018-density-exp.csv")
sep_ed <- read_csv("edited-ev-leaf-scans-sep-2018-density-exp.csv")


#### edit data ####

# indicate whether images were edited
jul$edited <- 0
jul_ed$edited <- 1
aug$edited <- 0
aug_ed$edited <- 1
sep$edited <- 0
sep_ed$edited <- 1

# remove edited from name
jul_ed$slice2 <- gsub("_edited", "", jul_ed$slice)
aug_ed$slice2 <- gsub("_edited", "", aug_ed$slice)
sep_ed$slice2 <- gsub("_edited", "", sep_ed$slice)

# remove unedited images and combine with edited
nrow(jul)
jul <- jul %>%
  filter(!(slice %in% jul_ed$slice2)) %>%
  full_join(select(jul_ed, -c(slice2)))
nrow(jul)

nrow(aug)
aug <- aug %>%
  filter(!(slice %in% aug_ed$slice2)) %>%
  full_join(select(aug_ed, -c(slice2)))
nrow(aug)

nrow(sep)
sep <- sep %>%
  filter(!(slice %in% sep_ed$slice2)) %>%
  full_join(select(sep_ed, -c(slice2)))
nrow(sep)

# assign month
jul$month <- "July"
aug$month <- "August"
sep$month <- "September"

# combine
dat <- rbind(jul, aug, sep)

# new columns
dat <- dat %>%
  mutate(
    plant = gsub(".*:","",slice) %>% gsub("_edited","",.),
    part = gsub(":.*$", "", slice),
  ) %>%
  mutate(
    part = recode(part, lesions = "lesion"),
    site = substr(plant, 1, 2),
    plot = substr(plant, 4, 5) %>% gsub("[^[:digit:]]", "", .),
    treatment = substr(plant, 4, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "W" = "water"),
    age = ifelse(grepl("EvA", plant), "adult", "seedling"),
    ID = substr(plant, 9, 13) %>% gsub(c("E|v|A|_"), "", .)
  ) %>%
  mutate(
    remove = ifelse(ID == "Bp", 1, 0), # remove extra leaf collected to exemplify Bipolaris infection in July
    ID = ifelse(ID == "" | ID == "Bp", "A", ID)
  ) %>%
  mutate(
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
  )

# check that the above worked as expected
unique(dat$plant)
unique(dat$part)
unique(dat$site)
unique(dat$plot)
unique(dat$treatment)
unique(dat$age)
unique(dat$ID)
unique(filter(dat, age == "adult")$ID)
unique(filter(dat, age == "seedling")$ID)
unique(dat$focal)

# spread by part
datw <- dat %>%
  select(-c(slice)) %>%
  gather(variable, value, -(edited:focal)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value)


#### check values ####

# leaf area
datw %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~age)

# manually check extremes
datw %>%
  filter(leaf_area.pix > 1.25e6)

datw %>%
  filter(leaf_area.pix == min(leaf_area.pix))

# lesion area 
datw %>%
  ggplot(aes(x = lesion_area.pix)) +
  geom_histogram() +
  facet_wrap(~age)

# manually check extremes
datw %>%
  filter(lesion_area.pix > 6e5)

# leaf objects
datw %>%
  ggplot(aes(x = leaf_objects)) +
  geom_histogram() +
  facet_wrap(~age)

# lesion objects
datw %>%
  ggplot(aes(x = lesion_objects)) +
  geom_histogram() +
  facet_wrap(~age)

# manually check extremes
datw %>%
  filter(lesion_objects == 0)

datw %>%
  filter(lesion_objects > 150)


#### export data ####

# set working directory
setwd("./processed-data")

# save file
write_csv(datw, "ev-leaf-scans-2018-desns-exp.csv")



#### modify data ####

datw <- datw %>%
  mutate(
    month = factor(month, levels = c("July", "August", "September")),
    damage = lesion_area.pix / leaf_area.pix
  )


#### figures ####

datw %>%
  ggplot(aes(x = treatment, y = damage)) +
  geom_boxplot()

datw %>%
  filter(treatment == "water") %>%
  ggplot(aes(x = month, y = damage)) +
  stat_summary(fun.data = "mean_cl_boot") +
  facet_grid(site~plot)

datw %>%
  filter(treatment == "fungicide") %>%
  ggplot(aes(x = month, y = damage)) +
  stat_summary(fun.data = "mean_cl_boot") +
  facet_grid(site~plot)
##### info ####

# file: ev-leaf-scans-data-processing
# author: Amy Kendig
# date last edited: 8/2/19
# goal: combine raw 2018 Elymus leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: LeafScanAnalysis_ev_ak_021819.ijm


#### set up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)

# import all raw data files
jul <- read_csv("./data/ev-leaf-scans-jul-2018-density-exp.csv")
jul_ed <- read_csv("./data/edited-ev-leaf-scans-jul-2018-density-exp.csv")
aug <- read_csv("./data/ev-leaf-scans-late-aug-2018-density-exp.csv")
aug_ed <- read_csv("./data/edited-ev-leaf-scans-late-aug-2018-density-exp.csv")
sep <- read_csv("./data/ev-leaf-scans-sep-2018-density-exp.csv")
sep_ed <- read_csv("./data/edited-ev-leaf-scans-sep-2018-density-exp.csv")


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
dat2 <- dat %>%
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
unique(dat2$plant)
unique(dat2$part)
unique(dat2$site)
unique(dat2$plot)
unique(dat2$treatment)
unique(dat2$age)
unique(dat2$ID)
unique(filter(dat2, age == "adult")$ID)
unique(filter(dat2, age == "seedling")$ID)
unique(dat2$focal)

# spread by part
datw <- dat2 %>%
  select(-c(slice)) %>%
  gather(variable, value, -(edited:focal)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  filter(!is.na(leaf_area.pix))


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


#### save data ####

eleaf <- datw




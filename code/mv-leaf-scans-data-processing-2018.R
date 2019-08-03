##### info ####

# file: mv-leaf-scans-data-processing
# author: Amy Kendig
# date last edited: 8/2/19
# goal: combine raw 2018 Microstegium leaf scan data and check for errors
# background: leaf scans were analyzed using FIJI, script: LeafScanAnalysis_mv_cw_011719.ijm


#### set-up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)

# import all raw data files
jul <- read_csv("./data/mv-leaf-scans-jul-2018-density-exp.csv")
jul_ed <- read_csv("./data/edited-mv-leaf-scans-jul-2018-density-exp.csv")
aug <- read_csv("./data/mv-leaf-scans-late-aug-2018-density-exp.csv")
aug_ed <- read_csv("./data/edited-mv-leaf-scans-late-aug-2018-density-exp.csv")
sep <- read_csv("./data/mv-leaf-scans-sep-2018-density-exp.csv")
sep_ed <- read_csv("./data/edited-mv-leaf-scans-sep-2018-density-exp.csv")


#### edit data ####

# indicate whether images were edited
jul$edited <- 0
jul_ed$edited <- 1
aug$edited <- 0
aug_ed$edited <- 1
sep$edited <- 0
sep_ed$edited <- 1

# remove edited from name
jul_ed$slice2 <- gsub("_edit", "", jul_ed$slice)
aug_ed$slice2 <- gsub("_edit", "", aug_ed$slice)
sep_ed$slice2 <- gsub("_edited", "", sep_ed$slice)

# remove unedited images and combine with edited
nrow(jul)
jul2 <- jul %>%
  filter(!(slice %in% jul_ed$slice2)) %>%
  full_join(select(jul_ed, -c(slice2)))
nrow(jul2)

nrow(aug)
aug2 <- aug %>%
  filter(!(slice %in% aug_ed$slice2)) %>%
  full_join(select(aug_ed, -c(slice2)))
nrow(aug2)

nrow(sep)
sep2 <- sep %>%
  filter(!(slice %in% sep_ed$slice2)) %>%
  full_join(select(sep_ed, -c(slice2)))
nrow(sep2)

# assign month
jul2$month <- "July"
aug2$month <- "August"
sep2$month <- "September"

# combine
dat <- rbind(jul2, aug2, sep2)

# new columns
dat2 <- dat %>%
  mutate(
    plant = gsub(".*:","",slice) %>% gsub("_edited","",.) %>% gsub("_edit","",.),
    part = gsub(":.*$", "", slice) %>% recode(lesions = "lesion", greens = "green"),
    site = substr(plant, 1, 2),
    treatment = substr(plant, 4, 6) %>% gsub("[^[:alpha:]]", "", .) %>% recode("F" = "fungicide", "W" = "water"),
    experiment = case_when(treatment == "T" ~ "transect",
                           TRUE ~ "density"),
    treatment = na_if(treatment, "T")
  ) 

# check that the above worked as expected
unique(dat2$plant)
unique(dat2$part)
unique(dat2$site)
unique(dat2$treatment)

# list of manually checked images with high senescence or overestimated damage
# this is from the R files folder in the Leaf Scans folder. I don't necessarily agree with these being removed, but they could be edited for more accurate lesion estimation
#datr <- tibble(
#  month = c(rep("July", 8), "August", rep("September", 5)),
#  plant = c("D1_4W_Mv_3","D1_10F_Mv_1","D4_9F_Mv_2", "D3_4F_Mv_1","D3_5F_Mv_3","D3_7F_Mv_1","D3_10W_Mv_3","D4_1F_Mv_1", "D2_1F_Mv2", "D2_4W_Mv_1","D3_4W_Mv_1","D3_6W_Mv_2","D4_7W_Mv_2","D4_8W_Mv_1"),
#  remove = 1
#)
datr <- tibble(month = "July", plant = "D3_4F_Mv_1", remove = 1)

# spread by part, remove green if leaf NA (one of the edited images had no green)
datw <- dat2 %>%
  select(-c(slice)) %>%
  gather(variable, value, -(edited:experiment)) %>%
  unite(temp, part, variable) %>%
  spread(temp, value) %>%
  filter(!is.na(leaf_area.pix)) %>%
  full_join(datr) %>%
  mutate(remove = replace_na(remove, 0))

# split by experiment and edit
dat_e <- datw %>%
  filter(experiment == "density") %>%
  mutate(
    plot = substr(plant, 4, 5) %>% gsub("[^[:digit:]]", "", .),
    ID = substr(plant, 9, 11) %>% gsub("[^[:digit:]]", "", .) %>% na_if(""),
    focal = ifelse(ID %in% c("1", "2", "3"), 1, 0)
  )

dat_t <- datw %>%
  filter(experiment == "transect") %>%
  mutate(
    transect = substr(plant, 5, 5),
    distance.m = substr(plant, 7, 8) %>% gsub("[^[:digit:]]", "", .)
  )

# check that the above worked as expected
unique(dat_e$plot)
unique(dat_e$ID)
filter(dat_e, is.na(ID)) %>% select(plant, focal) %>% unique() %>% data.frame()
unique(dat_e$focal)

unique(dat_t$transect)
unique(dat_t$distance.m)

# combine samples from the same plant

dat_t$plant # not needed for transect data

nrow(dat_e)
dat_e2 <- dat_e %>%
  group_by(edited, month, site, treatment, experiment, plot, ID, focal, remove) %>%
  summarise(
    green_area.pix = sum(green_area.pix, na.rm = T),
    green_objects = sum(green_objects, na.rm = T),
    leaf_area.pix = sum(leaf_area.pix, na.rm = T),
    leaf_objects = sum(leaf_objects, na.rm = T),
    lesion_area.pix = sum(lesion_area.pix, na.rm = T),
    lesion_objects = sum(lesion_objects, na.rm = T)
  ) %>%
  ungroup()
nrow(dat_e2)

# check for duplicates
dat_e2 %>% 
  group_by(month, site, plot, treatment, ID) %>%
  mutate(reps = n()) %>%
  filter(reps > 1)


#### check values ####

# leaf area
dat_e2 %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram() +
  facet_wrap(~focal, scales = "free")

dat_t %>%
  ggplot(aes(x = leaf_area.pix)) +
  geom_histogram()

# manually check extremes
dat_e2 %>%
  filter(focal == "0" & leaf_area.pix > 5.8e6) %>% data.frame()

dat_e2 %>%
  filter(focal == "1" & leaf_area.pix > 7e5) %>% data.frame()

dat_e2 %>%
  filter(focal == "1" & leaf_objects > 1) %>% data.frame()

# lesion area 
dat_e2 %>%
  ggplot(aes(x = lesion_area.pix)) +
  geom_histogram() +
  facet_wrap(~focal, scales = "free")

dat_t %>%
  ggplot(aes(x = lesion_area.pix)) +
  geom_histogram()

# manually check extremes
dat_e2 %>%
  filter(focal == "0" & lesion_area.pix > 1.5e6) %>% data.frame()

dat_e2 %>%
  filter(focal == "1" & lesion_area.pix > 3e5) %>% data.frame()

# lesion objects
dat_e2 %>%
  ggplot(aes(x = lesion_objects)) +
  geom_histogram() +
  facet_wrap(~focal, scales = "free")

# manually check extremes
dat_e2 %>%
  filter(focal == "0" & lesion_objects > 600) %>% data.frame()

dat_e2 %>%
  filter(focal == "1" & lesion_objects > 90) %>% data.frame()


#### save data ####

mleaf <- dat_e2
tleaf <- dat_t
##### info ####

# file: ev-survival-data-processing-2018
# author: Amy Kendig
# date last edited: 4/30/19
# goal: calculate Ev survival for summer 2018 and winter 2018/2019


#### set up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)

# run seed data files
source("./code/ev-seeds-data-processing-2018.R")

# clear everything except seed data
rm(list = setdiff(ls(), "eseeds"))

# import survival data
fjn <- read_csv("./data/focal-size-disease-jun-2018-density-exp.csv") # focal June, height and notes
fjl <- read_csv("./data/focal-size-disease-jul-2018-density-exp.csv") # focal July, height and notes
fea <- read_csv("./data/focal-status-early-aug-2018-density-exp.csv") # focal early August, status
ela <- read_csv("./data/all-disease-seeds-late-aug-2018-density-exp.csv") # all late August, no green and leaves tot
es <- read_csv("./data/ev-disease-seeds-sep-2018-density-exp.csv") # all Ev Sep, leaves tot
ew <- read_csv("./data/ev-winter-survival-apr-2019-density-exp.csv") # all Ev over winter
bw <- read_csv("./data/bg-ev-counts-apr-2019-density-exp.csv") # background counts over winter - not going to use these because I don't think the fall data are that reliable (individuals were difficult to find)
ejn <- read_csv("./data/ev-size-disease-jun-2018-litter-exp.csv") # litter Ev June, height
ejl <- read_csv("./data/ev-size-disease-jul-2018-litter-exp.csv") # litter Ev July, height
esl <- read_csv("./data/ev-disease-sep-2018-litter-exp.csv") # litter Ev Sep, leaves tot


#### edit data ####

# create survival columns

# June
unique(fjn$field_notes)
unique(subset(fjn, field_notes == "dead")$height.cm)
unique(subset(fjn, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjn <- fjn %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# July
unique(fjl$field_notes)
unique(subset(fjl, field_notes == "dead")$height.cm)
unique(subset(fjl, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjl <- fjl %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# early August
unique(fea$status)
fea <- fea %>%
  mutate(survival = ifelse(status %in% c("dead", "tree fell on plot"), 0, 1),
         age = case_when(sp == "Ev" & ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1)

# late August
unique(ela$no_green)
unique(subset(ela, no_green == 1)$leaves_tot)
unique(subset(ela, is.na(no_green))$leaves_tot)
head(subset(ela, is.na(no_green))) %>% data.frame() #one plant is not Ev (remove from all: D2.7.Water.Ev.A), the other lost its tag and is Mv 
unique(subset(ela, seeds_collected == 1)$no_green) # 1 or 0, will correct with actual seeds
ela <- ela %>%
  mutate(survival = ifelse(no_green == 0, 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# September
unique(subset(es, is.na(leaves_tot))$field_notes) # three are unsure if plant survived
es <- es %>%
  mutate(survival = ifelse(!is.na(leaves_tot), 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# over winter
ew <- ew %>%
  mutate(focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0),
         age = case_when(ID == "A" ~ "adult",
                         focal == 0 & plot > 7 ~ "adult",
                         TRUE ~ "seedling"),
         ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID))

# litter
ejn <- ejn %>%
  mutate(survival = ifelse(!is.na(height.cm), 1, 0),
         age = "adult",
         focal = 1)

ejl <- ejl %>%
  mutate(survival = ifelse(!is.na(height.cm), 1, 0),
         age = "adult",
         focal = 1)

unique(esl$leaves_tot)
esl <- esl %>%
  mutate(survival = ifelse(!is.na(leaves_tot), 1, 0),
         age = "adult",
         focal = 1)

# seeds
sum(is.na(eseeds$seeds)) # all included survived

# select data and add month
fjn <- fjn %>%
  filter(sp == "Ev") %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "June")

fjl <- fjl %>%
  filter(sp == "Ev") %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "July")

fea <- fea %>%
  filter(sp == "Ev") %>%
  select(site, plot, treatment, sp, age, ID, focal, survival) %>%
  mutate(month = "early August")

ela <- ela %>%
  filter(sp == "Ev") %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "late August")

es <- es %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "September")

ew <- ew %>%
  select(site, plot, treatment, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "April")

ejn <- ejn %>%
  select(site, plot, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "June")

ejl <- ejl %>%
  select(site, plot, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "July")

esl <- esl %>%
  select(site, plot, sp, age, ID, focal, survival, field_notes) %>%
  mutate(month = "September") 

# combine seeds across dates
eseeds <- eseeds %>%
  filter(ID_unclear == 0) %>%
  mutate(survival_seeds = 1) %>%
  select(site, plot, treatment, sp, age, ID, focal, survival_seeds) %>%
  unique()

# combine all data
d <- full_join(fjn, fjl) %>%
  full_join(fea) %>%
  full_join(ela) %>%
  full_join(es) %>%
  full_join(ew) %>%
  full_join(ejn) %>%
  full_join(ejl) %>%
  full_join(esl) %>%
  full_join(eseeds)

# remove compromised data
d <- d %>%
  filter(
    !(site == "D4" & plot == 8 & treatment == "fungicide") & 
      !(site == "D4" & plot == 10 & treatment == "fungicide") & 
      !(site == "D2" & plot == 7 & treatment == "water" & ID == "A") &
      !(site == "L2" & plot == 3)
  )

# change survival values if seeds were produced
d <- d %>%
  mutate(survival_seeds = replace_na(survival_seeds, 0),
         survival = ifelse(survival_seeds == 1, 1, survival))

# change survival values if plant was noted as surviving at a later month
sum(is.na(d$survival)) #14
d2 <- d %>%
  mutate(month = recode(month, "late August" = "lAug", "early August" = "eAug")) %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  mutate(September = ifelse(April %in% 1, 1, September),
         lAug = ifelse(September %in% 1, 1, lAug),
         eAug = ifelse(lAug %in% 1, 1, eAug),
         July = ifelse(eAug %in%1, 1, July),
         June = ifelse(July %in% 1, 1, June)) %>%
  gather("month", "survival", -c(site, plot, treatment, sp, age, ID, focal, survival_seeds)) %>%
  mutate(month = recode(month, "lAug" = "late August", "eAug" = "early August")) %>%
  full_join(select(d, site, plot, treatment, sp, age, ID, focal, month, field_notes)) %>%
  mutate(month = factor(month, levels = c("June", "July", "early August", "late August", "September", "April")))

# make sure it worked
d2 %>% 
  mutate(plant = paste(site, plot, treatment, ID, age, focal, sep = ".")) %>%
  ggplot(aes(x = as.numeric(month), y = survival, colour = plant)) +
  geom_line() +
  theme(legend.position = "none")
sum(is.na(d2$survival)) # 67

# check NA values
filter(d2, is.na(survival)) %>% select(field_notes) %>% unique()
filter(d2, is.na(survival) & is.na(field_notes)) %>% data.frame() # all individuals we didn't assess

# mean survival through summer
d2 %>%
  filter(month == "September" & !is.na(survival)) %>%
  group_by(age) %>%
  summarise(sum(survival) / length(survival))

# same, but remove litter plots
d2 %>%
  filter(month == "September" & !is.na(survival) & !is.na(treatment)) %>%
  group_by(age) %>%
  summarise(sum(survival) / length(survival))

# survival through summer and winter
d2 %>%
  filter(month == "April" & !is.na(survival)) %>%
  group_by(age) %>%
  summarise(sum(survival) / length(survival))

# survival through winter given summer survival
d2 %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April)) %>%
  group_by(age) %>%
  summarise(sum(April) / length(April))

# save data
esurv <- d2
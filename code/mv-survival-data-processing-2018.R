##### info ####

# file: mv-survival-data-processing-2018
# author: Amy Kendig
# date last edited: 6/17/19
# goal: calculate Mv survival for summer 2018


#### set up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)

# import survival data
fjn <- read_csv("./data/focal-size-disease-jun-2018-density-exp.csv") # focal June, height and notes
fjl <- read_csv("./data/focal-size-disease-jul-2018-density-exp.csv") # focal July, height and notes
fea <- read_csv("./data/focal-status-early-aug-2018-density-exp.csv") # focal early August, status
ala <- read_csv("./data/all-disease-seeds-late-aug-2018-density-exp.csv") # all late August, no green and leaves tot
ms <- read_csv("./data/mv-disease-sep-2018-density-exp.csv") # all Mv Sep, leaves tot


#### edit data ####

# create survival columns

# June
unique(fjn$field_notes)
unique(subset(fjn, field_notes == "dead")$height.cm)
unique(subset(fjn, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjn <- fjn %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         focal = 1)

# July
unique(fjl$field_notes)
unique(subset(fjl, field_notes == "dead")$height.cm)
unique(subset(fjl, is.na(height.cm))$field_notes) # dead or tree fell on plot (remove these plots from full dataset)
fjl <- fjl %>%
  mutate(survival = ifelse(is.na(height.cm), 0, 1),
         focal = 1)

# early August
unique(fea$status)
fea <- fea %>%
  mutate(survival = ifelse(status %in% c("dead", "tree fell on plot"), 0, 1),
         focal = 1)

# late August
unique(ala$no_green)
unique(subset(ala, no_green == 1)$leaves_tot)
unique(subset(ala, is.na(no_green))$leaves_tot)
head(subset(ala, is.na(no_green))) %>% data.frame() #one plant is not Ev (remove from all: D2.7.Water.Ev.A), the other lost its tag and is Mv (D3 2W Mv 1) 
ala <- ala %>%
  mutate(survival = ifelse(no_green == 0, 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0))

# September
unique(subset(ms, is.na(leaves_tot))$field_notes)
ms <- ms %>%
  mutate(survival = ifelse(!is.na(leaves_tot), 1, 0),
         focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0))

# select data and add month
fjn <- fjn %>%
  filter(sp == "Mv" & focal == 1) %>%
  select(site, plot, treatment, sp, ID, focal, survival, field_notes) %>%
  mutate(month = "June")

fjl <- fjl %>%
  filter(sp == "Mv" & focal == 1) %>%
  select(site, plot, treatment, sp, ID, focal, survival, field_notes) %>%
  mutate(month = "July")

fea <- fea %>%
  filter(sp == "Mv" & focal == 1) %>%
  select(site, plot, treatment, sp, ID, focal, survival) %>%
  mutate(month = "early August")

ala <- ala %>%
  filter(sp == "Mv" & focal == 1) %>%
  select(site, plot, treatment, sp, ID, focal, survival, field_notes) %>%
  mutate(month = "late August")

ms <- ms %>%
  filter(focal == 1) %>%
  select(site, plot, treatment, sp, ID, focal, survival, field_notes) %>%
  mutate(month = "September")

# combine all data
d <- full_join(fjn, fjl) %>%
  full_join(fea) %>%
  full_join(ala) %>%
  full_join(ms)

# remove compromised data
d <- d %>%
  filter(
    !(site == "D4" & plot == 8 & treatment == "fungicide") & 
      !(site == "D4" & plot == 10 & treatment == "fungicide") & 
      !(site == "D3" & plot == 2 & treatment == "water" & ID == "1" & month %in% c("late August", "September"))
  )

# change survival values if plant was noted as surviving at a later month
sum(is.na(d$survival)) # 0

d2 <- d %>%
  mutate(month = recode(month, "late August" = "lAug", "early August" = "eAug")) %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  mutate(lAug = ifelse(September %in% 1, 1, lAug),
         eAug = ifelse(lAug %in% 1, 1, eAug),
         July = ifelse(eAug %in%1, 1, July),
         June = ifelse(July %in% 1, 1, June)) %>%
  gather("month", "survival", -c(site, plot, treatment, sp, ID, focal)) %>%
  mutate(month = recode(month, "lAug" = "late August", "eAug" = "early August")) %>%
  full_join(select(d, site, plot, treatment, sp, ID, focal, month, field_notes)) %>%
  mutate(month = factor(month, levels = c("June", "July", "early August", "late August", "September")))

# make sure it worked
d2 %>% 
  mutate(plant = paste(site, plot, treatment, ID, focal, sep = ".")) %>%
  ggplot(aes(x = as.numeric(month), y = survival, colour = plant)) +
  geom_line() +
  theme(legend.position = "none")
sum(is.na(d2$survival)) # 2

# check NA values
filter(d2, is.na(survival)) %>% data.frame() # the one that lost its tag

# mean survival through summer
d2 %>%
  filter(month == "September" & !is.na(survival)) %>%
  summarise(sum(survival) / length(survival))

d2 %>%
  filter(month == "September" & !is.na(survival)) %>%
  group_by(site, treatment) %>%
  summarise(n = length(survival),
            surv = sum(survival) / length(survival))

# save data
msurv <- d2

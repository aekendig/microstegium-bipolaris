##### info ####

# file: fungal-isolation-by-treatment-2018
# author: Amy Kendig
# date last edited: 6/3/19
# goal: see if the treatment and sampling time affected isolation rate


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# run leaf scan files
source("./code/ev-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), "eleaf"))
source("./code/mv-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf", "eleaf")))

# import data
dj <- read_csv("./data/fungal-isolation-jul-2018-density-exp.csv")
da <- read_csv("./data/fungal-isolation-late-aug-2018-density-exp.csv")
ds <- read_csv("./data/fungal-isolation-sep-2018-density-exp.csv")
treat <- read_csv("./data/plot-treatments-2018-density-exp.csv")


#### edit data ####

# remove extra Ev leaves
eleaf <- eleaf %>%
  filter(remove == 0)

# examine conidiophore info
dj %>% select(observation, gigantea) %>% unique() %>% data.frame()
dj %>% filter(is.na(observation) & !is.na(gigantea)) %>% data.frame()
dj %>% filter(observation == "nothing") %>% data.frame()

da %>% select(observation, gigantea) %>% unique() %>% data.frame()

ds %>% select(observation, gigantea) %>% unique() %>% data.frame()

# remove NA symptoms/observations and non-NA gigantea
# no conidiophores, but WCH, should be yes
# senescing leaf
dj <- dj %>%
  mutate(gigantea = case_when(is.na(observation) ~ NA_character_,
                               observation == "no conidiophores, but WCH" ~ "Yes",
                               observation == "no conidiophores but WCH" ~ "Yes",
                               observation == "nothing" ~ NA_character_,
                               TRUE ~ gigantea))

# missing gigantea values
ds <- ds %>%
  mutate(gigantea = case_when(observation == "conidiophores" ~ "Yes",
                               observation == "no conidiophores" ~ "No",
                               TRUE ~ gigantea))

# combine data
d <- full_join(dj, da) %>% full_join(ds) %>% full_join(treat)

# summarise leaf scan data - indicate if leaves were collected from a plot
elsum <- eleaf %>%
  select(month, site, plot, treatment) %>%
  unique() %>%
  mutate(collect = 1,
         sp = "Ev",
         plot = as.numeric(plot),
         month = recode(month, August = "late August"))

mlsum <- mleaf %>%
  ungroup() %>%
  select(month, site, plot, treatment) %>%
  unique() %>%
  mutate(collect = 1,
         sp = "Mv",
         plot = as.numeric(plot),
         month = recode(month, August = "late August"))

# merge all data
# make an isolation column
# order months
d <- d %>%
  left_join(full_join(elsum, mlsum)) %>%
  mutate(isolation = ifelse(gigantea == "Yes", 1, 0),
         collect = replace_na(collect, 0),
         month = factor(month, levels = c("July", "late August", "September")),
         density_level = factor(density_level, levels = c("none", "low", "medium", "high")))


#### visualize ####

d %>%
  filter(collect == 1) %>%
  ggplot(aes(x = month, y = isolation, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", na.rm = T) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", na.rm = T, width = 0.2) +
  facet_wrap(~sp) +
  theme_bw()
# Ev highest in late August
# no data for Mv in late August, but higher in July than September

d %>% 
  filter(month == "late August" & sp == "Ev") %>%
  ggplot(aes(x = density_level, y = isolation, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", na.rm = T) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", na.rm = T, width = 0.2) +
  facet_wrap(~background, nrow = 4) +
  theme_bw()
# Elymus isolation affected by Mv density and Ev density , really high variation

d %>% 
  filter(month == "July" & sp == "Mv") %>%
  ggplot(aes(x = density_level, y = isolation, colour = treatment)) +
  stat_summary(fun.y = "mean", geom = "point", na.rm = T) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", na.rm = T, width = 0.2) +
  facet_wrap(~background, nrow = 4) +
  theme_bw()
# Mv isolation pretty independent of desnity, especially for water treatment


#### table of counts ####

d %>%
  group_by(sp, month) %>%
  summarise(iso_tot = sum(collect),
            iso_pos = sum(isolation, na.rm = T)) %>%
  mutate(iso_neg = iso_tot - iso_pos)

d %>%
  group_by(sp, month, treatment) %>%
  summarise(iso_tot = sum(collect),
            iso_pos = sum(isolation, na.rm = T)) %>%
  mutate(iso_neg = iso_tot - iso_pos)
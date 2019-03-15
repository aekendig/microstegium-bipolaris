##### info ####

# file: ev-seeds-data-processing
# author: Amy Kendig, Chris Wojan
# date last edited: 3/14/19
# goal: edit Ev seed data and check for errors
# background: spikelet counts for all samples, seed counts for a subset


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd("./data")

# import all raw data files
di <- read_csv("ev-spikelets-2018-density-litter-exp.csv")
di
de <- read_csv("ev-seed-subset-2018-density-exp.csv")
de


#### edit data ####

# check columns
unique(di$site)
unique(di$plot)
unique(di$plant)
filter(di, is.na(plant))
unique(di$collect_date)

unique(de$site)
unique(de$plot)
unique(de$treatment)
unique(de$age)
unique(de$ID)
unique(filter(de, is.na(ID))$age)

# add and modify columns in di
di <- di %>%
  mutate(
    plant = ifelse(is.na(plant), "A", plant),
    treatment = str_extract(plot, "[aA-zZ]+"),
    plot = str_extract(plot, "[0-9]+") %>% as.numeric(),
    month = case_when(
      collect_date == 20180711 ~ "July",
      collect_date > 20180711 & collect_date < 20180827 ~ "August",
      collect_date > 20180803 & collect_date < 20180924 ~ "late August",
      collect_date > 20180830 & collect_date < 20181022 ~ "September",
      collect_date == 20181022 ~ "October"
    )
  ) %>%
  mutate(
    age = ifelse(substr(plant, 1, 3) == "EVA" | plant == "A", "adult", "seedling"),
    ID = str_remove(plant, "EV"),
    treatment = recode(treatment, "W" = "water", "F" = "fungicide"),
  ) %>%
  mutate(
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
  ) %>%
  mutate(
    ID = ifelse(focal == 0 & age == "adult", str_remove(ID, "A"), ID)
  ) %>%
  rename(weigh_date = process_date)

# check that month worked
select(di, collect_date, month) %>% unique()

# add and modify columns in de
de <- de %>%
  mutate(
    treatment = recode(treatment, "W" = "water", "F" = "fungicide"),
    age = recode(age, "A" = "adult", "S" = "seedling"),
    ID = ifelse(is.na(ID), "A", ID),
    month = case_when(
      collect_date == 20180711 ~ "July",
      collect_date > 20180711 & collect_date < 20180827 ~ "August",
      collect_date > 20180803 & collect_date < 20180924 ~ "late August",
      collect_date > 20180830 & collect_date < 20181022 ~ "September",
      collect_date == 20181022 ~ "October"
    )
  ) %>%
  mutate(
    focal = ifelse(ID %in% c("1", "2", "3", "A"), 1, 0)
  ) %>%
  rename(seed_extraction_date = process_date) %>%
  mutate(collect_date = ifelse(site == "D3" & plot == 7 & treatment == "water" & month == "late August" & ID == "R2" & age == "seedling", 20180828, collect_date))

# merge
d <- full_join(di, de)
d

# notes
unique(d$notes)

# check some samples
d %>% filter(notes == "unknown plant")
d %>% filter(notes %in% c("1 OF 2", "2 OF 2")) %>% data.frame()

# add a removal column, edit notes, and check for duplicates
d <- d %>%
  mutate(
    remove = case_when(
      notes == "unknown plant" ~ 1,
      TRUE ~ 0), 
    notes = ifelse(notes %in% c("1 OF 2", "2 OF 2"), "2 samples", notes )
  ) %>%
  group_by(site, plot, plant, treatment, ID, age, month, focal) %>%
  mutate(reps = n())

# check for duplicate rows
filter(d, reps > 1)
# check with Chris about the labels on these envelopes

# check litter samples
filter(d, substr(site, 1, 1) == "L") %>% data.frame
filter(d, site == "L1" & month == "late August")
filter(d, site == "L2" & month == "late August")
filter(d, site == "L3" & month == "late August")
filter(d, site == "L4" & month == "late August")
filter(d, substr(site, 1, 1) == "L" & month == "September") %>% data.frame

# combine separated samples
group_cols <- d %>% select(site:weigh_date, notes:seed_extraction_date, remove, reps) %>% names
d2 <- d %>%
  filter(notes == "2 samples") %>%
  group_by(.dots = group_cols) %>%
  summarise(spikelets = sum(spikelets),
            spikelet_weight.mg = sum(spikelet_weight.mg),
            seeds = sum(seeds),
            seed_weight.mg = sum(seed_weight.mg)
   )  %>%
  full_join(filter(d, notes != "2 samples" | is.na(notes)))


#### check data ####

# spikelet number
d2 %>%
  ggplot(aes(x = spikelets)) +
  geom_histogram()

d2 %>%
  filter(spikelets > 200)
# adults during late August and September

# spikelet weight
d2 %>%
  ggplot(aes(x = spikelet_weight.mg)) +
  geom_histogram()

d2 %>%
  filter(spikelet_weight.mg > 3500)
# adults in late August

d2 %>%
  ggplot(aes(x = spikelets, y = spikelet_weight.mg, colour = age, shape = month)) +
  geom_point()

# large weight:number
d2 %>%
  filter(
    (spikelets < 20 & spikelet_weight.mg > 500) |
    (spikelets < 50 & spikelet_weight.mg > 1000) |
    (spikelet_weight.mg > 4000)
    )

# seed count
d2 %>%
  ggplot(aes(x = seeds)) +
  geom_histogram()

d2 %>%
  filter(seeds > 200)
# adult in September

# seed weight
d2 %>%
  ggplot(aes(x = seed_weight.mg)) +
  geom_histogram()

d2 %>%
  filter(seed_weight.mg > 600)
# same as above

# seed and spikelet
d2 %>%
  ggplot(aes(x = spikelets, y = seeds)) +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

# regression fit
m1 <- lm(seeds ~ spikelets, data = d2)
summary(m1) # R2 is 0.87

# seed and spikelet weight
d2 %>%
  ggplot(aes(x = spikelet_weight.mg, y = seed_weight.mg))  +
  geom_point(aes(colour = age, shape = treatment)) +
  geom_smooth(method = lm)

# regression fit
m2 <- lm(seed_weight.mg ~ spikelet_weight.mg, data = d2)
summary(m2) # R2 is 0.91

# seed count and spikelet weight
d2  %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = lm, size = 0.5, colour = "black")

# regression fit
m3 <- lm(seeds ~ spikelet_weight.mg, data = d2)
summary(m3) # R2 is 0.94

# seed count and spikelet weight again with site
d2 %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = treatment)) +
  geom_smooth(method = lm)

# residuals
d2 %>%
  filter(!is.na(seeds)) %>%
  ungroup() %>%
  mutate(resid = as.numeric(resid(m3))) %>%
  ggplot(aes(x = spikelet_weight.mg, y = resid)) +
  geom_point(aes(colour = site, shape = age)) + 
  facet_wrap(~treatment)
# may underestimate water more than fungicide

# seed count and spikelet weight by treatment
d2 %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = spikelet_weight.mg, y = seeds))  +
  geom_point(aes(colour = site, shape = age)) +
  facet_wrap(~treatment) +
  geom_smooth(method = lm)

# regression fit
m4 <- lm(seeds ~ spikelet_weight.mg * treatment, data = d2)
summary(m4)
m5 <- lm(seeds ~ spikelet_weight.mg + treatment, data = d2)
summary(m5)
AIC(m3, m4, m5)
# m3 is the best fit model and adding treatment hardly changes the slope of the relationship. it does decrease the estimate with water, but this effect seems unnecessary when looking at the figures. also, the simple model underestimates water.

# add ratio columns, modify date
d3 <- d2 %>%
  ungroup() %>%
  mutate(
    count_ratio = seeds / spikelets,
    weight_ratio = seed_weight.mg / spikelet_weight.mg, 
    cw_ratio = seeds / spikelet_weight.mg,
    collect_date = as.Date(as.character(collect_date), "%Y%m%d")
  ) %>%
  mutate(
    seeds = ifelse(is.na(seeds), coef(m3)[1] + spikelet_weight.mg * coef(m3)[2], seeds)
  )

# ratio over time
d3 %>%
  ggplot(aes(x = collect_date, y = cw_ratio)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) + 
  geom_smooth(method = "lm")
# not a strong trend

# ratio over spikelet count
d3 %>%
  ggplot(aes(x = spikelets, y = cw_ratio)) + 
  geom_point(aes(colour = age, shape = treatment), size = 2) +
  geom_smooth(method = "lm")
# hard tell because sample size gets smaller with more spikelets

# save data
eseeds = d3

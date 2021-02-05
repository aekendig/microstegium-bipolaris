##### info ####

# file: mv_seeds_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 2/4/21
# goal: edit Mv seed data and check for errors
# background: 10 samples collected per plot for the bags


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)

# import all raw data files
bag <- read_csv("./data/mv_bag_seed_2018_density_exp.csv")
bio <- read_csv("./intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")


#### edit data ####

# look at bag data
bag
unique(bag$bag_notes)
bag %>% filter(substr(bag_notes, 1, 7) == "Recount")

# modify columns
bag2 <- bag %>%
  mutate(recount = case_when(substr(bag_notes, 1, 9) == "Recounted" ~ 1,
                        TRUE ~ 0) %>% as.factor())

# look at biomass data 
bio

# modify columns
# soil sample is 0.05 x 0.25 m and quadrat was 0.49 x 0.25 m
bio2 <- bio %>%
  mutate(seeds_soil_total = seeds_soil * (0.49 / 0.05))


#### check data ####

# stems
bag2 %>%
  ggplot(aes(x = stems)) +
  geom_histogram()

# seeds
bag2 %>%
  ggplot(aes(x = seeds)) +
  geom_histogram()

bag2 %>% 
  filter(seeds > 175) %>%
  data.frame()
# highest one was double checked and others have similar weight, different counters

# seed weight
bag2 %>%
  ggplot(aes(x = seed_weight.g)) +
  geom_histogram()

bag2 %>% 
  filter(seed_weight.g > 0.175) %>%
  data.frame()
# high seed counts and different counters

# seeds and seed weight
bag2 %>%
  ggplot(aes(x = seeds, y = seed_weight.g)) +
  geom_point(data = filter(bag2, recount == "0")) +
  geom_point(data = filter(bag2, recount == "1"), colour = "red", size = 2)
# several of the extreme outliers were checked and none were entered incorrectly

filter(bag2, seed_weight.g == max(seed_weight.g, na.rm = T)) %>% data.frame()

# biomass seeds
bio2 %>%
  ggplot(aes(x = seeds_bio)) +
  geom_histogram()

filter(bio2, seeds_bio > 1500) # D1 10F

# soil seeds
bio2 %>%
  ggplot(aes(x = seeds_soil_total)) +
  geom_histogram()

filter(bio2, seeds_soil_total > 2500) # D3 10W and 1W

# seed relationship
bio2 %>%
  ggplot(aes(x = seeds_soil_total, y = seeds_bio, color = bio.g)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red")


#### combine data ####

# remove recount data
# if no stems were found in the bag, assume they had two
# combine bag data by plot
# use the average because not all plots had 10 bags
# add soil seeds
dat <- bag2 %>%
  filter(recount == "0") %>%
  mutate(stems = ifelse(stems == 0, 2, stems)) %>%
  group_by(site, plot, treatment) %>%
  summarise(seeds_per_bag = mean(seeds, na.rm = T),
            seeds_per_stem = mean(seeds/stems, na.rm = T)) %>%
  ungroup() %>%
  full_join(bio2 %>%
              select(site, plot, treatment, seeds_bio, seeds_soil_total) %>%
              rename(seeds_soil = seeds_soil_total))

# look at relationship between the seed sources
dat %>%
  select(seeds_per_bag, seeds_per_stem, seeds_bio, seeds_soil) %>%
  ggpairs()


#### output data ####
write_csv(dat, "./intermediate-data/mv_processed_seeds_2018_density_exp.csv")

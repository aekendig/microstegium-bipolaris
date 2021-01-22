##### info ####

# file: mv_seeds_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 1/21/21
# goal: edit Mv seed data and check for errors
# background: 10 samples collected per plot for the bags


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

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
                        TRUE ~ 0) %>% as.factor(),
    source = "bag")

# look at biomass data 
bio

# modify columns
# soil sample is 0.05 x 0.25 m and quadrat was 0.49 x 0.25 m
bio2 <- bio %>%
  mutate(seeds_soil_total = seeds_soil * (0.49 / 0.05),
         seeds = seeds_bio + seeds_soil_total,
         source = "bio")


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

# combine bag data by plant
# use the average because not all plots had 10 bags
# calculate expected value for 10 bags
# add soil seeds
dat <- bag2 %>%
  group_by(site, plot, treatment, source) %>%
  summarise(seeds = mean(seeds, na.rm = T) * 10) %>%
  ungroup() %>%
  full_join(bio2 %>%
              select(site, plot, treatment, seeds, source))

# make data wide
datw <- dat %>%
  spread(key = source, value = seeds) %>%
  rename(bag_seeds = bag, bio_seeds = bio) %>%
  mutate(total_seeds = bag_seeds + bio_seeds)

# look at relationship between the two seed sources
ggplot(datw, aes(x = bag_seeds, y = bio_seeds)) +
  geom_point()


#### output data ####
write_csv(datw, "./intermediate-data/mv_processed_seeds_2018_density_exp.csv")

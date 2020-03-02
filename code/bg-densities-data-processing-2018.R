##### info ####

# file: bg-densities-data-processing-2018
# author: Amy Kendig
# date last edited: 6/19/19
# goal: estimate actual background densities for summer 2018


#### set up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)

# run Ev survival file
source("./code/ev-survival-data-processing-2018.R")

# clear everything except final data
rm(list = setdiff(ls(), "esurv"))

# import data
jn <- read_csv("./data/bg_counts_jun_2018_density_exp.csv")
jl <- read_csv("./data/bg_counts_jul_2018_density_exp.csv")
ea <- read_csv("./data/bg_ev_counts_early_aug_2018_density_exp.csv")
la <- read_csv("./data/bg_ev_counts_late_aug_2018_density_exp.csv")


#### edit data ####

# June
jn
jn <- jn %>%
  mutate(month = "June",
         green = planted - missing - dead) %>%
  select(month, site, plot, treatment, planted, green)

# July
jl
jl <- jl %>%
  mutate(month = "July",
         green = planted - missing - dead) %>%
  select(month, site, plot, treatment, planted, green)

# early August
ea
ea <- ea %>%
  mutate(month = "early August",
         green = planted - missing - dead) %>%
  select(month, site, plot, treatment, planted, green)

# late August
la
la <- la %>%
  mutate(month = "late August") %>%
  select(month, site, plot, treatment, planted, green)

# month numbers
mn <- tibble(
  month = c("June", "July", "early August", "late August", "September"),
  week_num = c(4, 8, 12, 16, 20)
)

# tracked Ev
esurv
seed <- esurv %>%
  filter(focal == 0 & month != "April") %>%
  group_by(month, site, plot, treatment) %>%
  summarise(green = sum(survival, na.rm = T)) %>%
  left_join(mn) %>%
  ungroup()

# merge
d <- jn %>%
  full_join(jl) %>%
  full_join(ea) %>%
  full_join(la) %>%
  left_join(mn)


#### visualize ####

d %>%
  ggplot(aes(x = week_num, y = green)) +
  geom_line(aes(linetype = treatment)) +
  geom_point(data = seed, aes(shape = treatment)) +
  facet_grid(plot~site)

d %>%
  filter(plot > 4) %>%
  ggplot(aes(x = week_num, y = green)) +
  geom_line(aes(linetype = treatment)) +
  geom_point(data = seed, aes(shape = treatment)) +
  facet_grid(plot~site)

# use late August unless seed is higher


#### final data set ####

mvd <- d %>% 
  filter(plot %in% c(1, 2, 3) & month == "July") %>%
  mutate(counted_density = planted) %>%
  select(month, site, plot, treatment, counted_density)

evd <- seed %>%
  filter(month == "late August") %>%
  mutate(seed = green) %>%
  select(month, site, plot, treatment, seed) %>%
  full_join(filter(d, month == "late August")) %>%
  rowwise() %>%
  mutate(counted_density = max(seed, green, na.rm = T)) %>%
  select(month, site, plot, treatment, counted_density)

bgd <- full_join(mvd, evd)  

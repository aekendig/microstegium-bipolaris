##### info ####

# file: mv-seeds-data-processing
# author: Amy Kendig
# date last edited: 3/18/19
# goal: edit Mv seed data and check for errors
# background: 10 samples collected per plot

#### set up ####

# load packages
#library(tidyverse)

# import all raw data files
d <- read_csv("./data/mv-bag-seed-2018-density-exp.csv")
d


#### edit data ####

# check columns
unique(d$site)
unique(d$treatment)
unique(d$plot)
unique(d$bag)
unique(d$bag_notes)

# modify columns
d <- d %>%
  mutate(
    treatment = recode(treatment, "W" = "water", "F" = "fungicide"), 
    recount = case_when(substr(bag_notes, 1, 9) == "Recounted" ~ 1,
                        TRUE ~ 0) %>% as.factor(),
    sp = "Mv"
  ) 


#### check data ####

# stems
d %>%
  ggplot(aes(x = stems)) +
  geom_histogram()

# seeds
d %>%
  ggplot(aes(x = seeds)) +
  geom_histogram()

d %>% 
  filter(seeds > 175) %>%
  data.frame()
# highest one was double checked and others have similar weight, different counters

# seed weight
d %>%
  ggplot(aes(x = seed_weight.g)) +
  geom_histogram()

d %>% 
  filter(seed_weight.g > 0.175) %>%
  data.frame()
# high seed counts and different counters

# seeds and seed weight
d %>%
  ggplot(aes(x = seeds, y = seed_weight.g)) +
  geom_point(data = filter(d, recount == "0")) +
  geom_point(data = filter(d, recount == "1"), colour = "red", size = 2)
# several of the extreme outliers were checked and none were entered incorrectly

filter(d, seed_weight.g == max(seed_weight.g, na.rm = T)) %>% data.frame()

## Save data
mseeds = d
  
#### outputs ####

# mv_plot_biomass_seeds_2019_density_exp.csv
# plot_biomass_2019_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
bgBioD2Dat <- read_csv("data/bg_biomass_2019_density_exp.csv")
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 


#### background biomass ####

# combine plots
bgBioD2Dat2 <- bgBioD2Dat %>%
  group_by(site, plot, treatment) %>%
  summarise(biomass_bg = sum(biomass.g)) %>%
  ungroup()

# look at missing data
bgBioD2Dat2 %>%
  group_by(site) %>%
  count()
# none are missing -- there are no 1's

# add 0 data to 1 plots
bgBioD2Dat3 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 2),
               treatment = rep(c("water", "fungicide"), 4)) %>%
  mutate(plot = 1,
         biomass_bg = 0) %>%
  full_join(bgBioD2Dat2)

sum(is.na(bgBioD2Dat3$biomass_bg))


#### Ev seeds ####

# planted density
evPlantsD2Dat <- tibble(plot = 5:10,
                        planted = c(7, 11, 19, 3, 5, 9))

# format seeds
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>% # complete list of Ev
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  filter((plot %in% 5:7 & age == "seedling") | (plot %in% 8:10 & age == "adult")) %>% # only use the type planted as background
  group_by(site, plot, treatment, age) %>%
  summarise(seeds_per_plant = mean(seeds)) %>%
  ungroup() %>%
  left_join(evPlantsD2Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-c(seeds_per_plant, planted))

# check
filter(evSeedD2Dat2, is.na(seeds))

# figure
ggplot(evSeedD2Dat2, aes(treatment, seeds)) +
  geom_boxplot() +
  facet_wrap(~ age, scales = "free")


#### Mv seeds ####

# planted density
mvPlantsD2Dat <- tibble(plot = 1:4,
                        planted = c(3, 11, 19, 67))

# format seeds
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(age = "seedling")  %>%
  filter(plot %in% 1:4) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(seeds_per_plant = mean(seeds, na.rm = T)) %>% # plot average
  ungroup() %>%
  left_join(mvPlantsD2Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-c(seeds_per_plant, planted))

# check
filter(mvSeedD2Dat2, is.na(seeds))

# figure
ggplot(mvSeedD2Dat2, aes(treatment, seeds)) +
  geom_boxplot()


##### plot biomass ####

# use average of others in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat3 %>%
  left_join(mvBioD2Dat %>% # add focal biomass
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment) %>%
              summarise(biomass_foc_mv = sum(biomass_weight.g)) %>%
              ungroup() %>%
              full_join(evBioD2Dat %>%
                          filter(ID %in% c("1", "2", "3")) %>%
                          group_by(site, plot, treatment) %>%
                          mutate(weight_adj = mean(weight, na.rm = T)) %>%
                          ungroup() %>%
                          mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                                    TRUE ~ weight)) %>%
                          group_by(site, plot, treatment) %>%
                          summarise(biomass_foc_evS = sum(weight)) %>%
                          ungroup()) %>%
              full_join(evBioD2Dat %>%
                          filter(ID == "A") %>%
                          select(site, plot, treatment, weight) %>%
                          rename(biomass_foc_evA = weight))) %>%
  mutate(biomass.g_m2 = case_when(plot %in% 2:4 ~ biomass_bg + biomass_foc_mv,
                                  plot %in% 5:7 ~ biomass_bg + biomass_foc_evS,
                                  plot %in% 8:10 ~ biomass_bg + biomass_foc_evA),
         biomass_mv = case_when(plot %in% 2:4 ~ biomass.g_m2,
                                TRUE ~ biomass_foc_mv),
         biomass_evS = case_when(plot %in% 5:7 ~ biomass.g_m2,
                                TRUE ~ biomass_foc_evS),
         biomass_evA = case_when(plot %in% 8:10 ~ biomass.g_m2,
                                 TRUE ~ biomass_foc_evA))

# check
filter(plotBioD2Dat, is.na(biomass.g_m2))

# figure
ggplot(plotBioD2Dat, aes(treatment, biomass.g_m2)) +
  geom_boxplot()


##### combine data ####

dat <- mvSeedD2Dat2 %>%
  full_join(plotBioD2Dat %>%
              filter(plot %in% 1:4) %>%
              select(site, plot, treatment, biomass_mv, biomass_bg, biomass_foc_mv))


#### output ####

write_csv(dat, "intermediate-data/mv_plot_biomass_seeds_2019_density_exp.csv")
write_csv(plotBioD2Dat, "intermediate-data/plot_biomass_2019_density_exp.csv")
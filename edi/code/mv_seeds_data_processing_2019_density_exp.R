##### outputs ####

# mv_plant_level_seeds_2019_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
bio_seeds <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
bag_seeds <- read_csv("data/mv_bag_seeds_2019_density_exp.csv")
bio_seeds_conv <- read_csv("data/mv_biomass_seeds_counted_2019_density_exp.csv")
bag_notes <- read_csv("data/focal_disease_sep_2019_density_exp.csv")


#### edit data #### 

# check notes
unique(bio_seeds$process_notes) # resolved
unique(bag_seeds$process_notes) # missing and extra bags
unique(bag_notes$field_notes) 
filter(bag_notes, field_notes == "2 seed bags") # all Elymus

# add field notes to bags
bag_seeds2 <- full_join(bag_seeds, bag_notes %>%
                          filter(sp == "Mv") %>%
                          mutate(plant = as.numeric(ID)) %>%
                          select(site, plot, treatment, plant, leaves_tot, flower_bags))

# make sure missing bags are NA
filter(bag_seeds2, seeds == 0) # yes

# number of bags
bag_number <- bag_seeds2 %>%
  group_by(site, plot, treatment, plant) %>%
  summarise(bags_counted = sum(!is.na(seeds)),
            bags_added = sum(two_clips == "yes", na.rm = T),
            bags_reported = unique(flower_bags)) %>%
  mutate(bags_added = replace_na(bags_added, 0),
         bags_from_plant = bags_counted - bags_added)
bag_number

# mismatches
filter(bag_number, bags_from_plant != bags_reported) %>% data.frame()

# add notes
# adjust for issues
bag_seeds3 <- bag_seeds2 %>%
  mutate(analysis_notes = case_when(site == "D1" & plot == 10 & treatment == "water" ~ "bags mixed up among plants",
                                    site == "D1" & plot == 1 & treatment == "fungicide" ~ "bags put in D1 2F paper bag, switched back based on field notes",
                                    site == "D1" & plot == 2 & treatment == "fungicide" ~ "bags put in D1 1F paper bag, switched back based on field notes",
                                    site == "D2" & plot == 4 & treatment == "fungicide" ~ "bags mixed up among plants",
                                    site == "D2" & plot == 9 & treatment == "fungicide" ~ "bags mixed up among plants",
                                    site == "D2" & plot == 9 & treatment == "water" ~ "bags mixed up among plants",
                                    site == "D2" & plot == 10 & treatment == "fungicide" ~ "bags mixed up among plants",
                                    site == "D3" & plot == 4 & treatment == "water" ~ "bags mixed up among plants",
                                    site == "D3" & plot == 6 & treatment == "fungicide" ~ "bags mixed up among plants",
                                    site == "D3" & plot == 10 & treatment == "water" & plant == 3 ~ "all bags lost",
                                    site == "D4" & plot == 4 & treatment == "fungicide" ~ "bags mixed up among plants",
                                    site == "D4" & plot == 4 & treatment == "water" ~ "bags mixed up among plants",
                                    site == "D4" & plot == 9 & treatment == "fungicide" & plant == 2 ~ "one bag lost",
                                    site == "D4" & plot == 10 & treatment == "fungicide" & plant == 1 ~ "removed bags collected in September",
                                    site == "D4" & plot == 10 & treatment == "fungicide" & plant == 3 ~ "removed bags collected in September",
                                    TRUE ~ NA_character_),
         plot = case_when(site == "D1" & plot == 1 & treatment == "fungicide" ~ 2,
                          site == "D1" & plot == 2 & treatment == "fungicide" ~ 1,
                          TRUE ~ plot)) %>%
  filter(!(site == "D4" & plot == 10 & treatment == "fungicide" & plant == 1 & !is.na(process_notes)) &
           !(site == "D4" & plot == 10 & treatment == "fungicide" & plant == 3 & !is.na(process_notes))) 
# filter removes extra bags collected in September that may have been mis-labelled based on bag numbers

# check notes for extra bags
filter(bag_notes, site == "D4" & plot == 10 & treatment == "fungicide")
filter(bag_seeds, site == "D4" & plot == 10 & treatment == "fungicide") %>% data.frame()
filter(bag_seeds3, site == "D4" & plot == 10 & treatment == "fungicide") %>% data.frame()

# average bag seeds by plant
bag_seeds_plant <- bag_seeds3 %>%
  group_by(site, plot, treatment, sp, plant) %>%
  summarise(mean_flower_seeds = mean(seeds, na.rm = T),
            analysis_notes = unique(analysis_notes)) %>%
  ungroup() %>%
  mutate(mean_flower_seeds = replace_na(mean_flower_seeds, 0),
         mean_flower_seeds = case_when(analysis_notes == "all bags lost" ~ NA_real_,
                                       TRUE ~ mean_flower_seeds))

# stem seed weight conversion
stem_seed_mod <- lm(seeds ~ stem_seed_weight.g, data = bio_seeds_conv)
summary(stem_seed_mod)

# convert stem seeds
bio_seeds2 <- bio_seeds %>%
  mutate(stem_seeds = predict(stem_seed_mod, newdata = .))

# missing data
filter(bio_seeds2, is.na(stem_seeds)) # D2 2F Mv 2
filter(bio_seeds2, is.na(flowers)) # D2 2F Mv 2 and D3 4W Mv 1
filter(bag_seeds_plant, is.na(mean_flower_seeds)) # D3 10W Mv 3
filter(bag_seeds_plant, site == "D2" & plot == 2 & treatment == "fungicide" & plant == 2) # 0 flower seeds
filter(bag_seeds_plant, site == "D3" & plot == 4 & treatment == "water" & plant == 1) # does have flower seeds
filter(bio_seeds2, site == "D3" & plot == 10 & treatment == "water" & plant == 3) # does have flowers

# use averages from other plants for the plants missing data
d22fmv2_stem_seeds <- filter(bio_seeds2, site == "D2" & plot == 2 & treatment == "fungicide") %>%
  summarise(stem_seeds = mean(stem_seeds, na.rm = T)) %>%
  as.numeric()

d34wmv1_flowers <- filter(bio_seeds2, site == "D3" & plot == 4 & treatment == "water") %>%
  summarise(flowers = mean(flowers, na.rm = T)) %>%
  as.numeric()

d310wmv3_flower_seeds <- filter(bag_seeds_plant, site == "D3" & plot == 10 & treatment == "water") %>%
  summarise(mean_flower_seeds = mean(mean_flower_seeds, na.rm = T)) %>%
  as.numeric()

# plant-level data
mv_plant_seed_dat <- bio_seeds2 %>%
  mutate(stem_seeds_adj = case_when(site == "D2" & plot == 2 & treatment == "fungicide" & plant == 2 ~ d22fmv2_stem_seeds,
                                TRUE ~ stem_seeds),
         flowers_adj = case_when(site == "D3" & plot == 4 & treatment == "water" & plant == 1 ~ d34wmv1_flowers,
                             site == "D2" & plot == 2 & treatment == "fungicide" & plant == 2 ~ 0,
                             TRUE ~ flowers)) %>%
  full_join(bag_seeds_plant %>%
              mutate(mean_flower_seeds_adj = case_when(site == "D3" & plot == 10 & treatment == "water" & plant == 3 ~ d310wmv3_flower_seeds,
                                                   TRUE ~ mean_flower_seeds))) %>%
  mutate(flower_seeds = flowers * mean_flower_seeds,
         seeds = flower_seeds + stem_seeds,
         flower_seeds_adj = flowers_adj * mean_flower_seeds_adj,
         seeds_adj = flower_seeds_adj + stem_seeds_adj)


#### visualize ####

# variation in flower bags among plants within a plot
ggplot(bag_seeds_plant, aes(x = plot, y = mean_flower_seeds)) +
  geom_point() +
  facet_grid(site ~ treatment)
# similar

# variation in flowers among plants within a plot
ggplot(bio_seeds, aes(x = plot, y = flowers)) +
  geom_point() +
  facet_grid(site ~ treatment)
# similar

# counted biomass seeds
bio_seeds_conv %>%
  filter(!is.na(seeds)) %>%
  ggplot(aes(x = stem_seed_weight.g, y = seeds)) +
  geom_point() +
  geom_smooth(method = "lm")

# converted biomass seeds
ggplot(bio_seeds2, aes(x = stem_seed_weight.g, y = stem_seeds)) +
  geom_point()

# plant and plot data
ggplot(mv_plant_seed_dat, aes(x = plot, y = seeds)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(site ~ treatment)


#### output ####
write_csv(mv_plant_seed_dat, "intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")

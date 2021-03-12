##### info ####

# file: CFY_relationships_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/3/21
# goal: explore relationships related to the theory of constant final yield (CFY)


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(car) # for logit
library(DescTools) # for Gini

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import focal fitness data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
# ev_seeds_data_processing_2018.R and ev_seeds_data_processing_2019.R
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") # mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") # ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")

# import growth data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") # bg_biomass_data_processing_2019_density_exp.R
mvBioD1Dat <- read_csv("intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1)) %>%
  select(plot, treatment, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal))

filter(survD1Dat2, is.na(survival))
filter(survD1Dat, site == "D3" & plot == 2 & sp == "Mv" & ID == "1" & treatment == "water") %>%
  data.frame()
# lost track of plant

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3))) %>%
  mutate(sp = "Mv",
         age = "seedling") %>%
  full_join(plotsD %>%
              select(plot, treatment) %>%
              expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
              expand_grid(tibble(ID = c("1", "2", "3", "A"),
                                 age = c(rep("seedling", 3), "adult"))) %>%
              mutate(sp = "Ev"))

# focal survival
survD2Dat2 <- survD2Dat %>%
  filter(focal == 1) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(survival = 1/(length(unique(replace_date)) + 1)) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  mutate(survival = replace_na(survival, 1))

# severity data
sevD1Dat2 <- sevD1Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  mutate(late_aug_severity_adj = case_when(is.na(late_aug_severity) & !is.na(jul_severity) & !is.na(sep_severity) ~ (jul_severity + sep_severity)/2,
                                           TRUE ~ late_aug_severity)) %>%
  ungroup() %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"))

sevD2Dat2 <- sevD2Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  filter(ID %in% c("1", "2", "3", "A")) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  mutate(early_aug_severity_adj = case_when(is.na(early_aug_severity) & !is.na(jul_severity) & !is.na(late_aug_severity) ~ (jul_severity + late_aug_severity)/2,
                                            TRUE ~ early_aug_severity)) %>%
  ungroup() %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"))

# severity change
sevChangeD2Dat <- sevD2Dat %>%
  filter(ID %in% c("1", "2", "3", "A") & !(month %in% c("may", "sep"))) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         severity_logit = logit(severity, adjust = 0.001),
         month_num = case_when(month == "jun" ~ 0,
                               month == "jul" ~ 1,
                               month == "early_aug" ~ 2,
                               month == "late_aug" ~ 3)) %>%
  filter(!is.na(severity)) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(severity_change = as.numeric(coef(lm(severity_logit ~ month_num))[2])) %>%
  ungroup() %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"))
  
# seeds
evSeedD1Dat2 <- evSeedD1Dat %>%
  filter(focal == 1) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds),
            ID_unclear = sum(ID_unclear)) %>%
  ungroup() %>% 
  right_join(survD1Dat2 %>%
               filter(sp == "Ev")) %>%
  filter(ID_unclear == 0 | is.na(ID_unclear)) %>%
  mutate(seeds = replace_na(seeds, 0)) %>%
  select(-ID_unclear)

mvSeedD1Dat2 <- growthD1Dat %>%
  filter(sp == "Mv") %>%
  select(site:ID, tillers_jul) %>%
  left_join(mvSeedD1Dat) %>%
  mutate(seeds = seeds_per_stem * tillers_jul) %>%
  select(site:ID, seeds)

evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling"))

# plot biomass data
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 2
filter(evBioD2Dat, is.na(weight)) # 1 seedling
filter(bgBioD2Dat, is.na(biomass.g)) # 0

# use average of other species in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         sp_age = paste(sp, age, sep = "_")) %>%
  pivot_wider(-c(sp, age),
              names_from = sp_age,
              names_glue = "{sp_age}_biomass",
              values_from = biomass.g) %>%
  mutate(Mv_seedling_biomass = replace_na(Mv_seedling_biomass, 0),
         Ev_seedling_biomass = replace_na(Ev_seedling_biomass, 0),
         Ev_adult_biomass = replace_na(Ev_adult_biomass, 0)) %>%
  select(-none_seedling_biomass) %>%
  left_join(mvBioD2Dat %>%
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment, ) %>%
              summarise(Mv_seedling_biomass_foc = sum(biomass_weight.g)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID %in% c("1", "2", "3")) %>%
              group_by(site, plot, treatment) %>%
              mutate(weight_adj = mean(weight, na.rm = T)) %>%
              ungroup() %>%
              mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                        TRUE ~ weight)) %>%
              group_by(site, plot, treatment) %>%
              summarise(Ev_seedling_biomass_foc = sum(weight)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID == "A") %>%
              select(site, plot, treatment, weight) %>%
              rename(Ev_adult_biomass_foc = weight)) %>%
  mutate(Mv_seedling_biomass = Mv_seedling_biomass + Mv_seedling_biomass_foc,
         Ev_seedling_biomass = Ev_seedling_biomass + Ev_seedling_biomass_foc,
         Ev_adult_biomass = Ev_adult_biomass + Ev_adult_biomass_foc) %>%
  select(-c(Mv_seedling_biomass_foc, Ev_seedling_biomass_foc, Ev_adult_biomass_foc))

# plot-level data
plotD1Dat <- plotDens %>%
  full_join(mvBioD1Dat %>%
              select(site, plot, treatment, bio.g))

plotD2Dat <- plotDens %>%
  full_join(plotBioD2Dat)

# combine Mv data
mvD1Dat <- mvSeedD1Dat2 %>%
  full_join(growthD1Dat %>%
              filter(sp == "Mv") %>%
              select(site:ID, tiller_growth)) %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Mv") %>%
              select(-age)) %>%
  full_join(sevD1Dat2 %>%
              select(site:ID, late_aug_severity_adj, jul_severity) %>%
              filter(sp == "Mv")) %>%
  full_join(plotDens) %>%
  mutate(age = "seedling")

mvD2Dat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  rename("ID" = "plant") %>%
  mutate(ID = as.character(ID)) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g, stem_seed_weight.g) %>%
              rename("ID" = "plant") %>%
              mutate(ID = as.character(ID))) %>%
  full_join(survD2Dat2 %>%
              filter(sp == "Mv")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Mv")) %>%
  full_join(sevChangeD2Dat %>%
              filter(sp == "Mv")) %>%
  full_join(plotDens) %>%
  full_join(plotBioD2Dat)

# combine Ev data
evD1Dat <- evSeedD1Dat2 %>%
  full_join(growthD1Dat %>%
              filter(sp == "Ev") %>%
              select(site:ID, tiller_growth) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(sevD1Dat2 %>%
              select(site:ID, age, late_aug_severity_adj, jul_severity) %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens)

evD2Dat <- evSeedD2Dat2 %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD2Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(sevChangeD2Dat %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens) %>%
  full_join(plotBioD2Dat) %>%
  rename("biomass_weight.g" = "weight")


#### Mv-Mv figures ####

# total biomass by density, treatments and sites separated
plotD1Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = bio.g, color = treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site)

plotD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = Mv_seedling_biomass, color = treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site)

# total biomass by density, treatments separated
plotD1Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = bio.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

plotD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = Mv_seedling_biomass, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# plant biomass by density, treatments and sites separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# plant biomass by density, treatments separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# variation in plant biomass by density, treatments and site separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "Gini") +
  stat_summary(geom = "point", fun = "Gini") +
  facet_wrap( ~ site)

# variation in plant biomass by density, treatments separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  group_by(Mv_seedling_density, treatment) %>%
  summarise(biomass_var = Gini(biomass_weight.g, na.rm = T),
            lower = Gini(biomass_weight.g, conf.level = 0.95, na.rm = T)[2],
            upper = Gini(biomass_weight.g, conf.level = 0.95, na.rm = T)[3]) %>%
  ungroup() %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_var, color = treatment)) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point()

# does plant biomass match plot biomass?
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  group_by(plot, Mv_seedling_density, treatment, site, Mv_seedling_biomass) %>%
  summarise(biomass_weight.g = mean(biomass_weight.g, na.rm = T) * Mv_seedling_density) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = Mv_seedling_biomass, y = biomass_weight.g, color = treatment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(shape = as.factor(plot))) +
  facet_wrap(~ site)

# survival by density, treatments and sites separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = survival, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# survival by density, treatments and sites separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = survival, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by density, treatments and sites separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# severity change by density, treatments separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity by density, treatments separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = late_aug_severity, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by biomass, treatments separated
mvD2Dat %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_biomass, y = severity_change, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm")

# severity by severity change, treatments separated
mvD2Dat %>%
  ggplot(aes(x = severity_change, y = late_aug_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ treatment)

# biomass by severity, treatments separated, all plots (looked at all of the months)
mvD2Dat %>%
  group_by(site, plot, treatment) %>%
  summarise(severity = mean(severity_change, na.rm = T),
            plant_biomass = mean(biomass_weight.g, na.rm = T)) %>%
  ungroup() %>%
  ggplot(aes(x = severity, y = plant_biomass, color = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ plot, scales = "free")
# need to look at this after accounting for density effects on severity and biomass

# seeds by biomass
ggplot(mvD2Dat, aes(biomass_weight.g, seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm")


#### EvS-EvS figures ####

# total biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = Ev_seedling_biomass, color = treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site)

# total biomass by density, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = Ev_seedling_biomass, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# plant biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# plant biomass by density, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# does plant biomass match plot biomass?
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  group_by(plot, Ev_seedling_density, treatment, site, Ev_seedling_biomass) %>%
  summarise(biomass_weight.g = mean(biomass_weight.g, na.rm = T) * Ev_seedling_density) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = Ev_seedling_biomass, y = biomass_weight.g, color = treatment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(shape = as.factor(plot))) +
  facet_wrap(~ site)

# severity change by density, treatments and sites separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# severity change by density, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = Ev_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")


#### EvA-EvA figures ####

# total biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = Ev_adult_biomass, color = treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ site)

# total biomass by density, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = Ev_adult_biomass, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# plant biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# plant biomass by density, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# does plant biomass match plot biomass?
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 5:7)) %>%
  group_by(plot, Ev_adult_density, treatment, site, Ev_adult_biomass) %>%
  summarise(biomass_weight.g = mean(biomass_weight.g, na.rm = T) * Ev_adult_density) %>%
  ungroup() %>%
  unique() %>%
  ggplot(aes(x = Ev_adult_biomass, y = biomass_weight.g, color = treatment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(shape = as.factor(plot))) +
  facet_wrap(~ site)

# severity change by density, treatments and sites separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# severity change by density, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = Ev_adult_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")


#### Mv-EvS figures ####

# plant biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# plant biomass by density, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by density, treatments and sites separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# severity change by density, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by biomass, treatments separated
evD2Dat %>%
  filter(age == "seedling" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_biomass, y = severity_change, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm")


#### Mv-EvA figures ####

# plant biomass by density, treatments and sites separated
evD2Dat %>%
  filter(age == "adult" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# plant biomass by density, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by density, treatments and sites separated
evD2Dat %>%
  filter(age == "adult" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ site)

# severity change by density, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_density, y = severity_change, color = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# severity change by biomass, treatments separated
evD2Dat %>%
  filter(age == "adult" & plot %in% 1:4) %>%
  ggplot(aes(x = Mv_seedling_biomass, y = severity_change, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm")

##### info ####

# file: piecewise_sem_2019_density_exp
# author: Amy Kendig
# date last edited: 2/9/21
# goal: fit piecewise SEM to focal plot-level data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(car)
library(tidyverse)
library(piecewiseSEM)
library(lme4)
library(GGally)


# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import focal fitness data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") # mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") # ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")

# import biomass data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# import severity data
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R

# import environmental variables
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") # temp_humidity_data_processing_2019_density_exp


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

filter(sevD2Dat2, is.na(early_aug_severity) & !is.na(early_aug_severity_adj)) %>%
  nrow() # replaced 8 missing values with the adjustment

range(sevD2Dat2$jul_severity, na.rm = T) # includes 0 and 1
range(sevD2Dat2$early_aug_severity_adj, na.rm = T) # includes 0 

# plot-scale severity
plotSevD2Dat <- sevD2Dat2 %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(early_aug_severity_adj = mean(early_aug_severity_adj, na.rm = T),
            early_aug_severity = mean(early_aug_severity, na.rm = T),
            jul_severity = mean(jul_severity, na.rm = T)) %>%
  ungroup()  %>%
  pivot_wider(names_from = sp,
              values_from = c(jul_severity, early_aug_severity, early_aug_severity_adj),
              names_glue = "{sp}_{.value}")

sum(is.na(plotSevD2Dat$Ev_jul_severity)) # 0
range(plotSevD2Dat$Ev_jul_severity)
sum(is.na(plotSevD2Dat$Ev_early_aug_severity)) # 3
range(plotSevD2Dat$Ev_early_aug_severity, na.rm = T) 
sum(is.na(plotSevD2Dat$Ev_early_aug_severity_adj)) # 1
range(plotSevD2Dat$Ev_early_aug_severity_adj, na.rm = T) 
sum(is.na(plotSevD2Dat$Mv_jul_severity)) # 0
range(plotSevD2Dat$Mv_jul_severity) # includes 0
sum(is.na(plotSevD2Dat$Mv_early_aug_severity)) # 0
range(plotSevD2Dat$Mv_early_aug_severity) # includes 0
sum(is.na(plotSevD2Dat$Mv_early_aug_severity_adj)) # 0
range(plotSevD2Dat$Mv_early_aug_severity_adj) # includes 0

# edge severity
edgeSevD2Dat <- mvEdgeSevD2Dat %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity) %>%
  filter(month == "jul")

sum(is.na(edgeSevD2Dat$edge_severity)) # 0
range(edgeSevD2Dat$edge_severity) 

# seeds
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling"))

# environmental data
envD2Dat2 <- envD2Dat %>%
  select(site, plot, treatment, month, dew_intensity2, hrs_10_35) %>%
  rename(dew_intensity = dew_intensity2) %>%
  filter(month == "early_aug")

# check focal notes
unique(mvSeedD2Dat$process_notes)
unique(mvSeedD2Dat$analysis_notes)

# missing data?
filter(mvSeedD2Dat, is.na(seeds)) %>%
  data.frame()
# 3 plants
filter(mvBioD2Dat, is.na(biomass_weight.g)) %>%
  data.frame()
# 2 plants
filter(evSeedD2Dat, is.na(seeds)) %>%
  data.frame()
# 0 plants
filter(evBioD2Dat, is.na(weight)) %>%
  data.frame()
# 1 plant

# combine data
disDat <- plotSevD2Dat %>%
  select(-c(Ev_early_aug_severity_adj, Mv_early_aug_severity_adj)) %>%
  left_join(edgeSevD2Dat %>%
              select(-month)) %>%
  left_join(plotDens) %>%
  filter(treatment == "water") %>%
  left_join(envD2Dat2 %>%
              select(-month)) %>%
  drop_na() %>%
  mutate(Ev_aug_severity_t = logit(Ev_early_aug_severity, adjust = 0.001),
         Mv_aug_severity_t = logit(Mv_early_aug_severity, adjust = 0.001),
         Ev_jul_severity_t = logit(Ev_jul_severity, adjust = 0.001),
         Mv_jul_severity_t = logit(Mv_jul_severity, adjust = 0.001),
         log_dew = log(dew_intensity),
         edge_severity_t = logit(edge_severity, adjust = 0.001))

mvDat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  rename("ID" = "plant") %>%
  mutate(ID = as.character(ID)) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g) %>%
              rename("ID" = "plant") %>%
              mutate(ID = as.character(ID))) %>%
  full_join(survD2Dat2 %>%
              filter(sp == "Mv")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Mv")) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_biomass = log(biomass_weight.g),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(early_aug_severity_adj, adjust = 0.001),
         survival_t = logit(survival, adjust = 0.001))

evDat <- evSeedD2Dat2 %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD2Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens) %>%
  rename("biomass_weight.g" = "weight") %>%
  mutate(log_seeds = log(seeds + 1),
         log_biomass = log(biomass_weight.g),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(early_aug_severity_adj, adjust = 0.001),
         survival_t = logit(survival, adjust = 0.001))

# combine
datJuly <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, jul_severity_t, survival_t) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, jul_severity_t, survival_t)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()

datAug <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, aug_severity_t, survival_t) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, aug_severity_t, survival_t)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()


#### visualizations ####

# correlations in disease data
disDat %>%
  select(Ev_jul_severity_t, Mv_jul_severity_t, Ev_aug_severity_t, Mv_aug_severity_t, log_dew, hrs_10_35, edge_severity_t) %>%
  ggpairs()
# Ev and Mv July severity correlated
# hrs 10 to 35 very limited in values
# Mv july severity very low

# disease random effects
ggplot(disDat, aes(x = site, y = Mv_aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# higher in D4

ggplot(disDat, aes(x = site, y = Ev_aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# lower in D4

# site effects
ggplot(datJuly, aes(x = site, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev produce more seeds at D3, Mv at D1

ggplot(datJuly, aes(x = site, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# same as seeds

ggplot(datJuly, aes(x = site, y = survival_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# lower at D4 for both seedlings

ggplot(datJuly, aes(x = site, y = jul_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# all tend to have higher at D4

ggplot(datAug, aes(x = site, y = aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Mv has higher at D4, Ev seedlings lower at D1

# plot effects
ggplot(datJuly, aes(x = plot_f, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datJuly, aes(x = plot_f, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datJuly, aes(x = plot_f, y = survival_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datJuly, aes(x = plot_f, y = jul_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datAug, aes(x = plot_f, y = aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

# distributions
ggplot(datJuly, aes(x = log_seeds)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")
# could make seeds a binary variable for Ev seedlings

ggplot(datJuly, aes(x = log_biomass)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")

ggplot(datJuly, aes(x = jul_severity_t)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")
# Mv has a lot of zeros

ggplot(datAug, aes(x = aug_severity_t)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")


#### fit disease model ####

# initial fit
dis_mod1 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disDat),
  lmer(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disDat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  lmer(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t
)

# summary
summary(dis_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values
# P = 0.042

# update model
dis_mod2 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disDat),
  lmer(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disDat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  lmer(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disDat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t,
  edge_severity_t %~~% Mv_jul_severity_t
)

# summary
summary(dis_mod2)
# no missing links
# good model fit (P = 0.541)
# dew intensity increase August Mv severity
# July Mv severity decreased August Ev severity, but edge Mv severity increased it
# august and edge correlations significant and positive


#### Mv July model ####

# Mv dat
mvDatJuly <- datJuly %>%
  filter(sp == "Mv")

# initial fit
mv_july_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatJuly),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatJuly),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatJuly),
  log_biomass%~~%jul_severity_t
)
# error from biomass model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
summary(mv_july_mod1)
# Mv seedlings affect seeds
# p = 0.011
# better R-squared with random effects

# update model
mv_july_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatJuly),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatJuly),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatJuly),
  log_biomass%~~%jul_severity_t
)

# model summary
summary(mv_july_mod2)
# no missing links
# P = 0.279
# biomass increased seeds
# Mv density reduced seeds
# severity reduced survival
# Mv density reduced severity
# severity and biomass positively correlated


#### Mv August model ####

# Mv dat
mvDatAug <- datAug %>%
  filter(sp == "Mv")

# initial fit
mv_aug_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatAug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatAug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatAug),
  log_biomass%~~%aug_severity_t
)
# error from biomass model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
summary(mv_aug_mod1)
# Mv seedlings affect seeds
# p = 0.02
# better R-squared with random effects

# update model
mv_aug_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatAug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatAug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatAug),
  log_biomass%~~%aug_severity_t
)

# model summary
summary(mv_aug_mod2)
# no missing links
# P = 0.331
# biomass increased seeds
# Mv density reduced seeds
# fungicide reduced severity
# severity and biomass positively correlated


#### Ev seedling July model ####

# Ev seedling dat
evSDatJuly <- datJuly %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_july_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatJuly),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatJuly),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_biomass%~~%jul_severity_t
)

# model summary
summary(evS_july_mod1)
# Ev seedling density affects survival
# fungicide affects biomass and survival
# seeds are related to survival
# p = 0.0
# better R-squared with random effects

# update model
evS_july_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatJuly),
  lmer(survival_t ~ log_biomass + jul_severity_t + Ev_seedling_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_biomass%~~%jul_severity_t,
  log_seeds%~~%survival_t
)

# model summary
summary(evS_july_mod2)
# no missing links
# P = 0.592
# biomass increased seeds
# severity decreased survival
# Ev seedling density increased survival
# fungicide increased survival and biomass
# Mv seedling density reduced biomass
# biomass and severity negatively correlated
# seeds and survival negatively correlated


#### Ev seedling August model ####

# Ev seedling dat
evSDatAug <- datAug %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_aug_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  log_biomass%~~%aug_severity_t
)

# model summary
summary(evS_aug_mod1)
# fungicide affects biomass and survival
# seeds are related to survival
# p = 0.0
# better R-squared with random effects

# update model
evS_aug_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  lmer(survival_t ~ log_biomass + aug_severity_t + fungicide + (1|site/plot_f), data = evSDatAug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  log_biomass%~~%aug_severity_t,
  log_seeds%~~%survival_t
)

# model summary
summary(evS_aug_mod2)
# no missing links
# P = 0.278
# biomass increased seeds
# severity increased seeds
# fungicide increased survival and biomass
# Mv seedling density reduced biomass
# fungicide reduced severity
# seeds and survival negatively correlated


#### Ev adult July model ####

# Ev adult dat
evADatJuly <- datJuly %>%
  filter(sp == "Ev" & age == "adult")
# only three individuals didn't survive

# initial fit
evA_july_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site), data = evADatJuly),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatJuly),
  log_biomass%~~%jul_severity_t
)

# model summary
summary(evA_july_mod1)
# no missing links
# p = 0.57
# better R-squared with random effects
# biomass increased seeds
# Ev adult density reduced severity
# biomass and severity correlated


#### Ev adult August model ####

# Ev adult data
evADatAug <- datAug %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_aug_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site), data = evADatAug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatAug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatAug),
  log_biomass%~~%aug_severity_t
)
# error from severity model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
summary(evA_aug_mod1)
# no missing links
# p = 0.371
# better R-squared with random effects
# biomass increased seeds
# biomass and severity correlated
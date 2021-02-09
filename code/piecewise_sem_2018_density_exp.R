##### info ####

# file: piecewise_sem_2018_density_exp
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
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
# ev_seeds_data_processing_2018.R and ev_seeds_data_processing_2019.R
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018

# import growth data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R

# import environmental variables
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") # covariate_data_processing_2018_density_exp


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

sum(is.na(sevD1Dat2$late_aug_severity)) # 116
sum(is.na(sevD1Dat2$late_aug_severity_adj)) # 91

# plot-scale severity
plotSevD1Dat <- sevD1Dat2 %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(late_aug_severity_adj = mean(late_aug_severity_adj, na.rm = T),
            late_aug_severity = mean(late_aug_severity, na.rm = T),
            jul_severity = mean(jul_severity, na.rm = T)) %>%
  ungroup()  %>%
  pivot_wider(names_from = sp,
              values_from = c(jul_severity, late_aug_severity, late_aug_severity_adj),
              names_glue = "{sp}_{.value}")

sum(is.na(plotSevD1Dat$Ev_late_aug_severity)) # 15
sum(is.na(plotSevD1Dat$Ev_late_aug_severity_adj)) # 12
sum(is.na(plotSevD1Dat$Mv_late_aug_severity)) # 3
sum(is.na(plotSevD1Dat$Mv_late_aug_severity_adj)) # 1

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

# combine data
disDat <- plotSevD1Dat %>%
  select(-c(Ev_late_aug_severity_adj, Mv_late_aug_severity_adj)) %>%
  left_join(plotDens) %>%
  left_join(envD1Dat %>%
              select(site, plot, treatment, mv_inf_jul.prop)) %>%
  drop_na() %>%
  mutate(Ev_aug_severity_t = logit(Ev_late_aug_severity, adjust = 0.001),
         Mv_aug_severity_t = logit(Mv_late_aug_severity, adjust = 0.001),
         Ev_jul_severity_t = logit(Ev_jul_severity, adjust = 0.001),
         Mv_jul_severity_t = logit(Mv_jul_severity, adjust = 0.001),
         edge_severity_t = logit(mv_inf_jul.prop, adjust = 0.001),
         fungicide = ifelse(treatment == "water", 0, 1))

mvDat <- mvSeedD1Dat2 %>%
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
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(late_aug_severity_adj, adjust = 0.001),
         age = "seedling")

evDat <- evSeedD1Dat2 %>%
  full_join(growthD1Dat %>%
              filter(sp == "Ev") %>%
              select(site:ID, tiller_growth) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(sevD1Dat2 %>%
              select(site:ID, age, late_aug_severity_adj, jul_severity) %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(late_aug_severity_adj, adjust = 0.001))

# combine
datJuly <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, jul_severity_t, survival) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, jul_severity_t, survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()

datAug <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, aug_severity_t, survival) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, aug_severity_t, survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()


#### visualizations ####

# correlations in disease data
disDat %>%
  select(Ev_aug_severity_t, Mv_aug_severity_t, Ev_jul_severity_t, Mv_jul_severity_t, edge_severity_t) %>%
  ggpairs()
# Mv edge severity and July plot severity correlated

# disease random effects
ggplot(disDat, aes(x = site, y = Ev_jul_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# D1 lower

ggplot(disDat, aes(x = site, y = Mv_jul_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# D1 higher

ggplot(disDat, aes(x = site, y = Ev_aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# no difference

ggplot(disDat, aes(x = site, y = Mv_aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)
# D3 slightly higher

# site effects
ggplot(datJuly, aes(x = site, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# pretty close

ggplot(datJuly, aes(x = site, y = log_growth)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# lower for Ev seedlings at D3, for Mv at D2 and D3

ggplot(datJuly, aes(x = site, y = survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# lower for all at D3

ggplot(datJuly, aes(x = site, y = jul_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev have lower at D1

ggplot(datAug, aes(x = site, y = aug_severity_t)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# pretty similar

# plot effects
ggplot(datJuly, aes(x = plot_f, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datJuly, aes(x = plot_f, y = log_growth)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datJuly, aes(x = plot_f, y = survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev adult only has one mortality

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

ggplot(datJuly, aes(x = log_growth)) +
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
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide + (1 | site), data = disDat),
  lm(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide, data = disDat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disDat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disDat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t
)
# error from Ev severity (boundary fit is singular) and site variance was very low - removed from random effects

# summary
summary(dis_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values except for Ev August severity
# poor fit (P = 0.006)

# update model
dis_mod2 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide + (1 | site), data = disDat),
  lm(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide, data = disDat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disDat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disDat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t,
  edge_severity_t %~~% Mv_jul_severity_t
)

# summary
summary(dis_mod2)
# no missing links
# P = 0.686
# edge severity positively correlated with Mv july severity
# fungicide reduced Mv Aug  and July severity
# Aug severity negatively correlated


#### Mv July model ####

# Mv dat
mvDatJuly <- datJuly %>%
  filter(sp == "Mv")

# initial fit
mv_july_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatJuly),
  glmer(survival ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(mv_july_mod1)
# Mv seedling affects seeds
# P = 0.037
# random effects improved R-squared

# update fit
mv_july_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatJuly),
  glmer(survival ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(mv_july_mod2)
# P = 0.406
# growth increased seeds
# severity increased seeds
# density decreased seeds


#### Mv August model ####

# Mv dat
mvDatAug <- datAug %>%
  filter(sp == "Mv")

# initial fit
mv_aug_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = mvDatAug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = mvDatAug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatAug),
  log_growth%~~%aug_severity_t
)
# error from survival model (boundary fit is singular) and plot variance was very low - removed from random effects

# model summary
summary(mv_aug_mod1)
# Mv density affects seeds
# P = 0.03
# random effects increase R-squared

# update fit
mv_aug_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatAug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = mvDatAug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatAug),
  log_growth%~~%aug_severity_t
)

# model summary
summary(mv_aug_mod2)
# P = 0.372
# growth increased seeds
# Mv density decreases seeds
# fungicide decreases severity
# severity and growth negatively correlated


#### Ev seedling July model ####

# Ev seedling dat
evSDatJuly <- datJuly %>%
  filter(sp == "Ev" & age == "seedling") %>%
  mutate(seeds_bin = ifelse(log_seeds > log(1), 1, 0))

# initial fit
evS_july_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site/plot_f), data = evSDatJuly),
  glmer(survival ~ log_growth + jul_severity_t + (1|plot_f), data = evSDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_growth%~~%jul_severity_t
)
# boundary is singular error, not sure where it came from, but tried removing site from survival random effects

# model summary
summary(evS_july_mod1)
# warnings
# Mv seedling density affects seeds and survival
# Ev adult density affects seeds
# fungicide affects growth
# P = 0.001
# R-squared higher with random effects

# update model
evS_july_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatJuly),
  glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(evS_july_mod2)
# warnings gone
# P = 0.628
# severity increased seeds
# MV density decreased seeds
# growth increased survival
# Mv density decreased survival
# fungicide increased growth
# Ev seedling and adult density increased severity

# see if Ev adult density link can be removed (not sig)
evS_july_mod3 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = evSDatJuly),
  glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(evS_july_mod3)
# warnings returned - decided to leave in


#### Ev seedling July model, binary seeds ####

# initial fit
evS_july_bin_mod1 <- psem(
  glmer(seeds_bin ~ log_growth + jul_severity_t + (1|site/plot_f), data = evSDatJuly, family = binomial),
  glmer(survival ~ log_growth + jul_severity_t + (1|plot_f), data = evSDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_growth%~~%jul_severity_t
)
# boundary is singular error, not sure where it came from, but tried removing site from survival random effects

# model summary
summary(evS_july_bin_mod1)
# warnings
# Mv seedlings density affects seeds and survival
# Ev adult density affects seeds
# fungicide affects growth

evS_july_bin_mod2 <- psem(
  glmer(seeds_bin ~ log_growth + jul_severity_t + Mv_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatJuly, family = binomial),
  glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatJuly, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(evS_july_bin_mod2)
# warnings gone
# P = 0.757
# severity increased seeds
# growth increased survival
# Mv density increased survival
# fungicide increased growth
# Ev seedling and adult dnesity affected severity
# all results the same as non-binary model


#### Ev seedling August model ####

# Ev seedling dat
evSDatAug <- datAug %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_aug_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = evSDatAug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = evSDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  log_growth%~~%aug_severity_t
)
# convergence error from survival model and plot random effects very small
# singular boundary error from growth model and site random effects very small

# model summary
summary(evS_aug_mod1)
# warnings
# fungicide affected growth
# P = 0.081
# random effects increased R-squared

# update model
evS_aug_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = evSDatAug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = evSDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  log_growth%~~%aug_severity_t
)

# model summary
summary(evS_aug_mod2)
# warnings

# try removing random effects from survival, don't increase R-squared much
evS_aug_mod3 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatAug),
  glm(survival ~ log_growth + aug_severity_t, data = evSDatAug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = evSDatAug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatAug),
  log_growth%~~%aug_severity_t
)

# model summary
summary(evS_aug_mod3)
# fungicide increased growth


#### Ev adult July model ####

# Ev adult dat
evADatJuly <- datJuly %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_july_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site), data = evADatJuly),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(evA_july_mod1)
# Ev seedling density affected seeds
# P = 0.102
# random effects increased R-squared

# update fit
evA_july_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Ev_seedling_density + (1|site), data = evADatJuly),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatJuly),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatJuly),
  log_growth%~~%jul_severity_t
)

# model summary
summary(evA_july_mod2)
# P = 0.305
# Ev seedling increased seeds
# Mv decreased growth
# severity and growth negatively correlated


#### Ev adult August model ####

# Ev adult data
evADatAug <- datAug %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_aug_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site), data = evADatAug),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatAug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatAug),
  log_growth%~~%aug_severity_t
)
# error from severity model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
summary(evA_aug_mod1)
# Ev seedlings affected seeds
# P = 0.37
# random effects only increased R-squared of seeds (nearly zero for growth)

# update fit
evA_aug_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + Ev_seedling_density + (1|site), data = evADatAug),
  lm(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatAug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatAug),
  log_growth%~~%aug_severity_t
)

# model summary
summary(evA_aug_mod2)
# P = 0.906
# severity increased seeds
# Ev density increased seeds
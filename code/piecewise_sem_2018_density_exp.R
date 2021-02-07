##### info ####

# file: piecewise_sem_2019_density_exp
# author: Amy Kendig
# date last edited: 2/4/21
# goal: fit piecewise SEM to focal plot-level data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(piecewiseSEM)
library(lme4)
library(nlme)
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
survD1Datb <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal))

filter(survD1Datb, is.na(survival))
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
evSeedD1Datb <- evSeedD1Dat %>%
  filter(focal == 1) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds),
            ID_unclear = sum(ID_unclear)) %>%
  ungroup() %>% 
  right_join(survD1Datb %>%
               filter(sp == "Ev")) %>%
  filter(ID_unclear == 0 | is.na(ID_unclear)) %>%
  mutate(seeds = replace_na(seeds, 0)) %>%
  select(-ID_unclear)

mvSeedD1Datb <- growthD1Dat %>%
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
              select(site, plot, treatment, soil_moisture_jun.prop, mv_inf_jul.prop)) %>%
  mutate(asr_Ev_aug_severity = asin(sqrt(Ev_late_aug_severity)),
         asr_Mv_aug_severity = asin(sqrt(Mv_late_aug_severity)),
         asr_Ev_jul_severity = asin(sqrt(Ev_jul_severity)),
         asr_Mv_jul_severity = asin(sqrt(Mv_jul_severity)),
         asr_edge_severity = asin(sqrt(mv_inf_jul.prop)),
         asr_soil_moisture = asin(sqrt(soil_moisture_jun.prop)),
         fungicide = ifelse(treatment == "water", 0, 1)) %>%
  drop_na()

mvDat <- mvSeedD1Datb %>%
  full_join(growthD1Dat %>%
              filter(sp == "Mv") %>%
              select(site:ID, tiller_growth)) %>%
  full_join(survD1Datb %>%
              filter(sp == "Mv") %>%
              select(-age)) %>%
  full_join(sevD1Dat2 %>%
              select(site:ID, late_aug_severity_adj) %>%
              filter(sp == "Mv")) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         asr_severity = asin(sqrt(late_aug_severity_adj)),
         age = "seedling")

evDat <- evSeedD1Datb %>%
  full_join(growthD1Dat %>%
              filter(sp == "Ev") %>%
              select(site:ID, tiller_growth) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD1Datb %>%
              filter(sp == "Ev")) %>%
  full_join(sevD1Dat2 %>%
              select(site:ID, age, late_aug_severity_adj) %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         asr_severity = asin(sqrt(late_aug_severity_adj)))

# combine
datSeeds <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, asr_severity, survival) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_growth, asr_severity, survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()


#### visualizations ####

# correlations in disease data
disDat %>%
  select(asr_Ev_jul_severity, asr_Mv_jul_severity, asr_Ev_aug_severity, asr_Mv_aug_severity, asr_soil_moisture, asr_edge_severity) %>%
  ggpairs()
# Mv edge severity and plot severity correlated

# disease histogram
ggplot(disDat, aes(x = asr_Mv_aug_severity)) +
  geom_histogram()

ggplot(disDat, aes(x = asr_Ev_aug_severity)) +
  geom_histogram()

# disease random effects
ggplot(disDat, aes(x = site, y = asr_Mv_aug_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)

ggplot(disDat, aes(x = site, y = asr_Ev_aug_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2)

# site effects
ggplot(datSeeds, aes(x = site, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev produce more seeds at D3 and D4, Mv slightly more at D1

ggplot(datSeeds, aes(x = site, y = log_growth)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Mv slightly more at D1 and D4

ggplot(datSeeds, aes(x = site, y = asr_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev seedlings have lower disease at D3 and Mv the highest at D3

ggplot(datSeeds, aes(x = site, y = survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev have no variation in survival
# mv has high survival at two sites

# plot effects
ggplot(datSeeds, aes(x = plot_f, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datSeeds, aes(x = plot_f, y = log_growth)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datSeeds, aes(x = plot_f, y = asr_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datSeeds, aes(x = plot_f, y = survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

# distributions
ggplot(datSeeds, aes(x = log_seeds)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")
# could make seeds a binary variable for Ev seedlings

ggplot(datSeeds, aes(x = log_growth)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")

ggplot(datSeeds, aes(x = asr_severity)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")


#### fit disease model ####

# check random effects
dis_lme1 <- lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, random = ~ 1 | site, data = disDat)
summary(dis_lme1)
# very small random effects
dis_lme2 <- lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, random = ~ 1 | site, data = disDat)
summary(dis_lme2)
# very small random effects
evsev_lme1 <- lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat)
summary(evsev_lme1)
mvsev_lme1 <- lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat)
summary(mvsev_lme1)

# initial fit
dis_mod1 <- psem(
  lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity
)

# summary
summary(dis_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values for July severities but not August
# poor fit (P = 0.004)

# update model
dis_mod2 <- psem(
  lm(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, data = disDat),
  lm(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity,
  asr_edge_severity%~~%asr_Mv_jul_severity
)

# summary
summary(dis_mod2)
# no missing links
# good model fit (P = 0.588)
# soil moisture reduced Mv august severity
# edge severity positively correlated with Mv july severity

# multigroup model
dis_mod3 <- multigroup(dis_mod2, group = "treatment")

# summary
dis_mod3
# P = 0.027, poor model fit

# remove non-significant correlations and re-try multigroup model
dis_mod4 <- psem(
  lm(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, data = disDat),
  lm(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  asr_edge_severity%~~%asr_Mv_jul_severity
)

dis_mod5 <- multigroup(dis_mod4, group = "treatment")
dis_mod5
# same model

# try fungicide in the model
# initial fit
dis_mod6 <- psem(
  lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity
)

summary(dis_mod6)
# no effectof random effects on August Ev severity
# edge effect on Mv severity: add correlation
# poor model fit (P = 0.007)

# update model
dis_mod7 <- psem(
  lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + asr_soil_moisture + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity,
  asr_edge_severity%~~%asr_Mv_jul_severity
)

summary(dis_mod7)
# no omitted links
# good model fit (P = 0.476)


#### Mv seeds models ####

# Mv dat
mvDatSeeds <- datSeeds %>%
  filter(sp == "Mv")

# check random effects
mv_seeds_lme1 <- lme(log_seeds ~ log_growth + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_seeds_lme1)
mv_surv_lm1 <- glm(survival ~ log_growth + asr_severity, family = binomial, data = mvDatSeeds)
summary(mv_surv_lm1)
mv_growth_lme1 <- lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_growth_lme1)
mv_severity_lme1 <- lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_severity_lme1)

# initial fit
mv_seeds_mod1 <- psem(
  lme(log_seeds ~ log_growth + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds),
  glm(survival ~ log_growth + asr_severity, data = mvDatSeeds),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds),
  log_growth%~~%asr_severity
)

# model summary
summary(mv_seeds_mod1)
# better R-squared with random effects
# d sep test: Mv density effect on seeds
# Fisher's C: P = 0.045, poor fit
# correlation not significant

# update model
mv_seeds_mod2 <- psem(
  lme(log_seeds ~ log_growth + asr_severity + Mv_seedling_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  glm(survival ~ log_growth + asr_severity, data = mvDatSeeds),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds),
  log_growth%~~%asr_severity
)

# summary
summary(mv_seeds_mod2)
# all omitted paths are supported
# better fit by Fisher's C (P = 0.373)
# random effects improve the R-squared values
# growth increases seeds
# Mv density reduces seeds
# fungicide reduces disease


#### Ev seedling seeds models ####

# Ev seedling dat
evSDatSeeds <- datSeeds %>%
  filter(sp == "Ev" & age == "seedling")

# check random effects
evS_seeds_lme1 <- lme(log_seeds ~ log_growth + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_seeds_lme1)
evS_growth_lme1 <- lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_growth_lme1)
evS_severity_lme1 <- lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_severity_lme1)

# initial fit
evS_seeds_mod1 <- psem(
  lme(log_seeds ~ log_growth + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  log_growth%~~%asr_severity
)

# model summary
summary(evS_seeds_mod1)
# better R-squared with random effects
# d sep test: fungicide directly affects growth
# Fisher's C: P = 0.009, poor fit

evS_seeds_mod2 <- psem(
  lme(log_seeds ~ log_growth + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  log_growth%~~%asr_severity
)

# model summary
summary(evS_seeds_mod2)
# better R-squared with random effects
# d sep test: all omitted paths supported
# Fisher's C: P = 0.141, good fit
# fungicide increases growth

# binomial model

# add column to data
evSDatSeeds2 <- evSDatSeeds %>%
  mutate(seeds_bin = as.numeric(log_seeds > 0))

# check model
evS_seeds_lme2 <- glmer(seeds_bin ~ log_growth + asr_severity + (1|plot_f), family = binomial, data = evSDatSeeds2)
summary(evS_seeds_lme2)

# initial fit
evS_bin_mod1 <- psem(
  glmer(seeds_bin ~ log_growth + asr_severity + (1|plot_f), family = binomial, data = evSDatSeeds2),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  log_growth%~~%asr_severity
)

# model summary
summary(evS_bin_mod1)
# better R-squared with random effects
# d sep test: fungicide directly affects growth and seeds
# Fisher's C: P = 0.004, poor fit

# update model
evS_bin_mod2 <- psem(
  glmer(seeds_bin ~ log_growth + asr_severity + fungicide + (1|plot_f), family = binomial, data = evSDatSeeds2),
  lme(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  log_growth%~~%asr_severity
)

# model summary
summary(evS_bin_mod2)
# better R-squared with random effects
# d sep test: all omitted links supported
# Fisher's C: P = 0.33, good fit
# fungicide increases growth and reduces seeds


#### Ev adult seeds models ####

# Ev adult dat
evADatSeeds <- datSeeds %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_seeds_mod1 <- psem(
  lm(log_seeds ~ log_growth + asr_severity, data = evADatSeeds),
  lm(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatSeeds),
  lm(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatSeeds),
  log_growth%~~%asr_severity
)

# model summary
summary(evA_seeds_mod1)
# d sep test: Ev seedling density affects seeds
# Fisher's C: P = 0.126, good fit

# update model
evA_seeds_mod2 <- psem(
  lm(log_seeds ~ log_growth + asr_severity + Ev_seedling_density, data = evADatSeeds),
  lm(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatSeeds),
  lm(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatSeeds),
  log_growth%~~%asr_severity
)

summary(evA_seeds_mod2)
# d sep test: support for all omitted links
# Fisher's C: P = 0.404, good fit
# severity increases seeds
# Ev seedling density increases seeds

#### start here: remove soil moisture from disease model, update graphs based on results ###

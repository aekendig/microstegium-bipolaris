##### info ####

# file: piecewise_sem_2018_2019_density_exp
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
library(cowplot)

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

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R

# import environmental variables
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") # covariate_data_processing_2018_density_exp
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") # temp_humidity_data_processing_2019_density_exp

# inverse-logit function (B. Bolker, https://stackoverflow.com/questions/23845283/logit-transformation-backwards)
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}


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

sum(is.na(sevD1Dat2$late_aug_severity)) # 116
sum(is.na(sevD1Dat2$late_aug_severity_adj)) # 91

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

# environmental data
envD2Dat2 <- envD2Dat %>%
  select(site, plot, treatment, month, dew_intensity2, hrs_10_35) %>%
  rename(dew_intensity = dew_intensity2) %>%
  filter(month == "early_aug")


# combine disease data
disD1Dat <- plotSevD1Dat %>%
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

disD2Dat <- plotSevD2Dat %>%
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
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(late_aug_severity_adj, adjust = 0.001),
         age = "seedling")

mvD2Dat <- mvSeedD2Dat %>%
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
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_growth = log(tiller_growth + 2),
         jul_severity_t = logit(jul_severity, adjust = 0.001),
         aug_severity_t = logit(late_aug_severity_adj, adjust = 0.001))

evD2Dat <- evSeedD2Dat2 %>%
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

# combine July data
datD1July <- mvD1Dat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_growth, jul_severity, jul_severity_t, survival) %>%
  full_join(evD1Dat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_growth, jul_severity, jul_severity_t, survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()

datD2July <- mvD2Dat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_biomass, jul_severity, jul_severity_t, survival, survival_t) %>%
  full_join(evD2Dat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_biomass, jul_severity, jul_severity_t, survival, survival_t)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()

# combine August data
datD1Aug <- mvD1Dat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_growth, late_aug_severity_adj, aug_severity_t, survival) %>%
  full_join(evD1Dat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_growth, late_aug_severity_adj, aug_severity_t, survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()

datD2Aug <- mvD2Dat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_biomass, early_aug_severity_adj,aug_severity_t, survival_t) %>%
  full_join(evD2Dat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, seeds, log_seeds, log_biomass, early_aug_severity_adj,aug_severity_t, survival_t)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()


#### visualizations ####

# see separate year R scripts


#### fit Y1 disease model ####

# initial fit
dis_y1_mod1 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide + (1 | site), data = disD1Dat),
  lm(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide, data = disD1Dat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disD1Dat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disD1Dat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t
)
# error from Ev severity (boundary fit is singular) and site variance was very low - removed from random effects

# summary
# summary(dis_y1_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values except for Ev August severity
# poor fit (P = 0.006)

# update model
dis_y1_mod2 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide + (1 | site), data = disD1Dat),
  lm(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide, data = disD1Dat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disD1Dat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1 | site), data = disD1Dat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t,
  edge_severity_t %~~% Mv_jul_severity_t
)

# summary
# summary(dis_y1_mod2)
# no missing links
# P = 0.686
# edge severity positively correlated with Mv july severity
# fungicide reduced Mv Aug  and July severity
# Aug severity negatively correlated


#### fit Y2 disease model ####

# initial fit
dis_y2_mod1 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat),
  lmer(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  lmer(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t
)

# summary
# summary(dis_y2_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values
# P = 0.042

# update model
dis_y2_mod2 <- psem(
  lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat),
  lmer(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat),
  lmer(Ev_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  lmer(Mv_jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  lmer(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1 | site), data = disD2Dat),
  Ev_jul_severity_t %~~% Mv_jul_severity_t,
  Ev_aug_severity_t %~~% Mv_aug_severity_t,
  edge_severity_t %~~% Mv_jul_severity_t
)

# summary
# summary(dis_y2_mod2)
# no missing links
# good model fit (P = 0.541)
# dew intensity increase August Mv severity
# July Mv severity decreased August Ev severity, but edge Mv severity increased it
# august and edge correlations significant and positive


#### Mv July Y1 model ####

# Mv dat
mvDatD1July <- datD1July %>%
  filter(sp == "Mv")

# initial fit
mv_july_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatD1July),
  glmer(survival ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatD1July, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(mv_july_y1_mod1)
# Mv seedling affects seeds
# P = 0.037
# random effects improved R-squared

# update fit
mv_july_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD1July),
  glmer(survival ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatD1July, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(mv_july_y1_mod2)
# P = 0.406
# growth increased seeds
# severity increased seeds
# density decreased seeds


#### Mv July Y2 model ####

# Mv dat
mvDatD2July <- datD2July %>%
  filter(sp == "Mv")

# initial fit
mv_july_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatD2July),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatD2July),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatD2July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD2July),
  log_biomass%~~%jul_severity_t
)
# error from biomass model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
# summary(mv_july_y2_mod1)
# Mv seedlings affect seeds
# p = 0.011
# better R-squared with random effects

# update model
mv_july_y2_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD2July),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatD2July),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatD2July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD2July),
  log_biomass%~~%jul_severity_t
)

# model summary
# summary(mv_july_y2_mod2)
# no missing links
# P = 0.279
# biomass increased seeds
# Mv density reduced seeds
# severity reduced survival
# Mv density reduced severity
# severity and biomass positively correlated


#### Mv August Y1 model ####

# Mv dat
mvDatD1Aug <- datD1Aug %>%
  filter(sp == "Mv")

# initial fit
mv_aug_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = mvDatD1Aug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = mvDatD1Aug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatD1Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD1Aug),
  log_growth%~~%aug_severity_t
)
# error from survival model (boundary fit is singular) and plot variance was very low - removed from random effects

# model summary
# summary(mv_aug_y1_mod1)
# Mv density affects seeds
# P = 0.03
# random effects increase R-squared

# update fit
mv_aug_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD1Aug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = mvDatD1Aug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = mvDatD1Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD1Aug),
  log_growth%~~%aug_severity_t
)

# model summary
# summary(mv_aug_y1_mod2)
# P = 0.372
# growth increased seeds
# Mv density decreases seeds
# fungicide decreases severity
# severity and growth negatively correlated


#### Mv August Y2 model ####

# Mv dat
mvDatD2Aug <- datD2Aug %>%
  filter(sp == "Mv")

# initial fit
mv_aug_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatD2Aug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatD2Aug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatD2Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD2Aug),
  log_biomass%~~%aug_severity_t
)
# error from biomass model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
# summary(mv_aug_y2_mod1)
# Mv seedlings affect seeds
# p = 0.02
# better R-squared with random effects

# update model
mv_aug_y2_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD2Aug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = mvDatD2Aug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDatD2Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = mvDatD2Aug),
  log_biomass%~~%aug_severity_t
)

# model summary
# summary(mv_aug_y2_mod2)
# no missing links
# P = 0.331
# biomass increased seeds
# Mv density reduced seeds
# fungicide reduced severity
# severity and biomass positively correlated


#### Ev seedling July Y1 model ####

# Ev seedling dat
evSDatD1July <- datD1July %>%
  filter(sp == "Ev" & age == "seedling") %>%
  mutate(seeds_bin = ifelse(log_seeds > log(1), 1, 0))

# initial fit
evS_july_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site/plot_f), data = evSDatD1July),
  glmer(survival ~ log_growth + jul_severity_t + (1|plot_f), data = evSDatD1July, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July),
  log_growth%~~%jul_severity_t
)
# boundary is singular error, not sure where it came from, but tried removing site from survival random effects

# model summary
# summary(evS_july_y1_mod1)
# warnings
# Mv seedling density affects seeds and survival
# Ev adult density affects seeds
# fungicide affects growth
# P = 0.001
# R-squared higher with random effects

# update model
evS_july_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatD1July),
  glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatD1July, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(evS_july_y1_mod2)
# warnings gone
# P = 0.628
# severity increased seeds
# MV density decreased seeds
# growth increased survival
# Mv density decreased survival
# fungicide increased growth
# Ev seedling and adult density increased severity

# see if Ev adult density link can be removed (not sig)
evS_july_y1_mod3 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = evSDatD1July),
  glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatD1July, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(evS_july_y1_mod3)
# warnings returned - decided to leave in


#### Ev seedling July Y2 model ####

# Ev seedling dat
evSDatD2July <- datD2July %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_july_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatD2July),
  lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatD2July),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatD2July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2July),
  log_biomass%~~%jul_severity_t
)

# model summary
# summary(evS_july_y2_mod1)
# Ev seedling density affects survival
# fungicide affects biomass and survival
# seeds are related to survival
# p = 0.0
# better R-squared with random effects

# update model
evS_july_y2_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site/plot_f), data = evSDatD2July),
  lmer(survival_t ~ log_biomass + jul_severity_t + Ev_seedling_density + fungicide + (1|site/plot_f), data = evSDatD2July),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2July),
  log_biomass%~~%jul_severity_t,
  log_seeds%~~%survival_t
)

# model summary
# summary(evS_july_y2_mod2)
# no missing links
# P = 0.592
# biomass increased seeds
# severity decreased survival
# Ev seedling density increased survival
# fungicide increased survival and biomass
# Mv seedling density reduced biomass
# biomass and severity negatively correlated
# seeds and survival negatively correlated


#### Ev seedling August Y1 model ####

# Ev seedling dat
evSDatD1Aug <- datD1Aug %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_aug_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatD1Aug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = evSDatD1Aug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = evSDatD1Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1Aug),
  log_growth%~~%aug_severity_t
)
# convergence error from survival model and plot random effects very small
# singular boundary error from growth model and site random effects very small

# model summary
# summary(evS_aug_y1_mod1)
# warnings
# fungicide affected growth
# P = 0.081
# random effects increased R-squared

# update model
evS_aug_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatD1Aug),
  glmer(survival ~ log_growth + aug_severity_t + (1|site), data = evSDatD1Aug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = evSDatD1Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1Aug),
  log_growth%~~%aug_severity_t
)

# model summary
# summary(evS_aug_y1_mod2)
# warnings

# try removing random effects from survival, don't increase R-squared much
evS_aug_y1_mod3 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site/plot_f), data = evSDatD1Aug),
  glm(survival ~ log_growth + aug_severity_t, data = evSDatD1Aug, family = binomial),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = evSDatD1Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1Aug),
  log_growth%~~%aug_severity_t
)

# model summary
# summary(evS_aug_y1_mod3)
# fungicide increased growth


#### Ev seedling August Y2 model ####

# Ev seedling dat
evSDatD2Aug <- datD2Aug %>%
  filter(sp == "Ev" & age == "seedling")

# initial fit
evS_aug_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatD2Aug),
  lmer(survival_t ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatD2Aug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site/plot_f), data = evSDatD2Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2Aug),
  log_biomass%~~%aug_severity_t
)

# model summary
# summary(evS_aug_y2_mod1)
# fungicide affects biomass and survival
# seeds are related to survival
# p = 0.0
# better R-squared with random effects

# update model
evS_aug_y2_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site/plot_f), data = evSDatD2Aug),
  lmer(survival_t ~ log_biomass + aug_severity_t + fungicide + (1|site/plot_f), data = evSDatD2Aug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2Aug),
  lmer(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD2Aug),
  log_biomass%~~%aug_severity_t,
  log_seeds%~~%survival_t
)

# model summary
# summary(evS_aug_y2_mod2)
# no missing links
# P = 0.278
# biomass increased seeds
# severity increased seeds
# fungicide increased survival and biomass
# Mv seedling density reduced biomass
# fungicide reduced severity
# seeds and survival negatively correlated


#### Ev adult July Y1 model ####

# Ev adult dat
evADatD1July <- datD1July %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_july_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + (1|site), data = evADatD1July),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(evA_july_y1_mod1)
# Ev seedling density affected seeds
# P = 0.102
# random effects increased R-squared

# update fit
evA_july_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + jul_severity_t + Ev_seedling_density + (1|site), data = evADatD1July),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD1July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatD1July),
  log_growth%~~%jul_severity_t
)

# model summary
# summary(evA_july_y1_mod2)
# P = 0.305
# Ev seedling increased seeds
# Mv decreased growth
# severity and growth negatively correlated


#### Ev adult July Y2 model ####

# Ev adult dat
evADatD2July <- datD2July %>%
  filter(sp == "Ev" & age == "adult")
# only three individuals didn't survive

# initial fit
evA_july_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + jul_severity_t + (1|site), data = evADatD2July),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD2July),
  lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatD2July),
  log_biomass%~~%jul_severity_t
)

# model summary
# summary(evA_july_y2_mod1)
# no missing links
# p = 0.57
# better R-squared with random effects
# biomass increased seeds
# Ev adult density reduced severity
# biomass and severity correlated


#### Ev adult August Y1 model ####

# Ev adult data
evADatD1Aug <- datD1Aug %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_aug_y1_mod1 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + (1|site), data = evADatD1Aug),
  lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD1Aug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatD1Aug),
  log_growth%~~%aug_severity_t
)
# error from severity model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
# summary(evA_aug_y1_mod1)
# Ev seedlings affected seeds
# P = 0.37
# random effects only increased R-squared of seeds (nearly zero for growth)

# update fit
evA_aug_y1_mod2 <- psem(
  lmer(log_seeds ~ log_growth + aug_severity_t + Ev_seedling_density + (1|site), data = evADatD1Aug),
  lm(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatD1Aug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatD1Aug),
  log_growth%~~%aug_severity_t
)

# model summary
# summary(evA_aug_y1_mod2)
# P = 0.906
# severity increased seeds
# Ev density increased seeds


#### Ev adult August Y2 model ####

# Ev adult data
evADatD2Aug <- datD2Aug %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_aug_y2_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + aug_severity_t + (1|site), data = evADatD2Aug),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD2Aug),
  lm(aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatD2Aug),
  log_biomass%~~%aug_severity_t
)
# error from severity model (boundary fit is singular) and site variance was very low - removed from random effects

# model summary
# summary(evA_aug_y2_mod1)
# no missing links
# p = 0.371
# better R-squared with random effects
# biomass increased seeds
# biomass and severity correlated


#### figure settings ####

# figure theme
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 12, hjust = -0.1),
        legend.background = element_blank())

# color palette
col_pal = c("black", "gray")
col_pal2 = c("#39568CFF", "3CBB75FF")
line_pal1 = c("solid", "dashed")
line_pal2 = c("dashed", "solid")
shape_pal = c(16, 17, 15)


#### fig: Mv density and Mv seeds, July models ####

# models
summary(mv_july_y1_mod2, intercepts = T)
mv_seeds_jul_y1_mod <- lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD1July)
summary(mv_seeds_jul_y1_mod)

summary(mv_july_y2_mod2, intercepts = T)
mv_seeds_jul_y2_mod <- lmer(log_seeds ~ log_biomass + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = mvDatD2July)
summary(mv_seeds_jul_y2_mod)

# predictions
mv_seeds_jul_y1_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(log_growth = mean(mvDatD1July$log_growth),
         jul_severity_t = mean(mvDatD1July$jul_severity_t)) %>%
  mutate(log_seeds = predict(mv_seeds_jul_y1_mod, newdata = ., re.form = NA),
         seeds = exp(log_seeds) - 1)

mv_seeds_jul_y2_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(log_biomass = mean(mvDatD2July$log_biomass),
         jul_severity_t = mean(mvDatD2July$jul_severity_t)) %>%
  mutate(log_seeds = predict(mv_seeds_jul_y2_mod, newdata = ., re.form = NA),
         seeds = exp(log_seeds) - 1)

mv_seeds_jul_dat <- mv_seeds_jul_y1_dat %>%
  mutate(year = "Year 1") %>%
  full_join(mv_seeds_jul_y2_dat %>%
              mutate(year = "Year 2"))

# figure
mv_seeds_dens_fig <- mvDatD1July %>% 
  mutate(year = "Year 1") %>%
  full_join(mvDatD2July %>%
              mutate(year = "Year 2")) %>%
  ggplot(aes(x = Mv_seedling_density, y = log_seeds, color = year)) +
  # geom_point(alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun.data = mean_se, position = position_dodge(1), size = 2, shape = shape_pal[1]) +
  geom_line(data = mv_seeds_jul_dat) +
  xlab("Mv density") +
  ylab("log(Mv seeds)") +
  ggtitle("C") +
  scale_color_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = c(0.78, 0.83),
        legend.title = element_blank())


#### fig: Mv density and Ev S seeds, July models ####

# models
summary(evS_july_y1_mod3, intercepts = T)
evS_seeds_jul_y1_mod <- lmer(log_seeds ~ log_growth + jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = evSDatD1July)
summary(evS_seeds_jul_y1_mod)

summary(evS_july_y2_mod2, intercepts = T)
# int: 1.5653
# compound path mediated by biomass (product of coefficients)
(evS_seeds_jul_y2_coef <- -0.0078*1.0463)
summary(lmer(log_seeds ~ jul_severity_t + Mv_seedling_density + (1|site/plot_f), data = evSDatD2July))
evS_seeds_jul_y2_int <- 1.8764

# predictions
evS_seeds_jul_y1_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(log_growth = mean(evSDatD1July$log_growth),
         jul_severity_t = mean(evSDatD1July$jul_severity_t)) %>%
  mutate(log_seeds = predict(evS_seeds_jul_y1_mod, newdata = ., re.form = NA))

evS_seeds_jul_y2_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(log_seeds = evS_seeds_jul_y2_int + evS_seeds_jul_y2_coef * Mv_seedling_density)

evS_seeds_jul_dat <- evS_seeds_jul_y1_dat %>%
  mutate(year = "Year 1") %>%
  full_join(evS_seeds_jul_y2_dat %>%
              mutate(year = "Year 2"))

# figure
ev_seeds_dens_fig <- evSDatD1July %>% 
  mutate(year = "Year 1") %>%
  full_join(evSDatD2July %>%
              mutate(year = "Year 2")) %>%
  ggplot(aes(x = Mv_seedling_density, y = log_seeds, color = year)) +
  # geom_point(alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "point", fun.data = mean_se, size = 2, shape = shape_pal[2]) +
  geom_line(data = evS_seeds_jul_dat) +
  xlab("Mv density") +
  ylab("log(Ev seedling seeds)") +
  ggtitle("D") +
  scale_color_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = "none",
        legend.title = element_blank())


#### fig: Mv density and Ev A growth, July models ####

# models
summary(evA_july_y1_mod2, intercepts = T) 
evA_growth_jul_y1_mod <- lmer(log_growth ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD1July)
summary(evA_growth_jul_y1_mod)

summary(evA_july_y2_mod1, intercepts = T) 
evA_growth_jul_y2_mod <- lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|site), data = evADatD2July)
summary(evA_growth_jul_y2_mod)

# predictions
evA_growth_jul_y1_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(Ev_seedling_density = mean(evSDatD1July$Ev_seedling_density),
         Ev_adult_density = mean(evSDatD1July$Ev_adult_density)) %>%
  mutate(log_growth = predict(evA_growth_jul_y1_mod, newdata = ., re.form = NA))

evA_growth_jul_y2_dat <- tibble(Mv_seedling_density = 3:67) %>%
  mutate(Ev_seedling_density = mean(evSDatD2July$Ev_seedling_density),
         Ev_adult_density = mean(evSDatD2July$Ev_adult_density)) %>%
  mutate(log_growth = predict(evA_growth_jul_y2_mod, newdata = ., re.form = NA))

evA_growth_jul_dat <- evA_growth_jul_y1_dat %>%
  mutate(year = "Year 1") %>%
  full_join(evA_growth_jul_y2_dat %>%
              mutate(year = "Year 2"))

# figure
ev_growth_dens_fig <- evADatD1July %>% 
  mutate(year = "Year 1") %>%
  full_join(evADatD2July %>%
              mutate(year = "Year 2",
                     log_growth = log_biomass)) %>%
  ggplot(aes(x = Mv_seedling_density, y = log_growth, color = year)) +
  # geom_point(alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "point", fun.data = mean_se, size = 2, shape = shape_pal[3]) +
  geom_line(data = evA_growth_jul_dat, aes(linetype = year)) +
  xlab("Mv density") +
  ylab("log(Ev adult growth)") +
  ggtitle("E") +
  scale_color_manual(values = col_pal) +
  scale_linetype_manual(values = line_pal1) +
  fig_theme +
  theme(legend.position = "none",
        legend.title = element_blank())


#### fig: severity and Mv survival, July models ####

# # models
# summary(mv_july_y1_mod2, intercepts = T)
# mv_surv_jul_y1_mod <- glmer(survival ~ log_growth + jul_severity_t + (1|site/plot_f), data = mvDatD1July, family = binomial)
# summary(mv_surv_jul_y1_mod)
# 
# summary(mv_july_y2_mod2, intercepts = T)
# mv_surv_jul_y2_mod <- lmer(survival_t ~ log_biomass + jul_severity_t + (1|site/plot_f), data = mvDatD2July)
# summary(mv_surv_jul_y2_mod)
# 
# # predictions
# mv_surv_jul_y1_dat <- tibble(jul_severity = seq(min(mvDatD1July$jul_severity), max(mvDatD1July$jul_severity), length.out = 100)) %>%
#   mutate(log_growth = mean(mvDatD1July$log_growth),
#          jul_severity_t = logit(jul_severity, adjust = 0.001)) %>%
#   mutate(survival = predict(mv_surv_jul_y1_mod, newdata = ., type = "response", re.form = NA))
# 
# mv_surv_jul_y2_dat <- tibble(jul_severity = seq(min(mvDatD2July$jul_severity), max(mvDatD2July$jul_severity), length.out = 100)) %>%
#   mutate(log_biomass = mean(mvDatD2July$log_biomass),
#          jul_severity_t = logit(jul_severity, adjust = 0.001)) %>%
#   mutate(survival = predict(mv_surv_jul_y2_mod, newdata = ., re.form = NA) %>%
#            inv.logit(., a = 0.001))
# 
# mv_surv_jul_dat <- mv_surv_jul_y1_dat %>%
#   mutate(year = "Year 1") %>%
#   full_join(mv_surv_jul_y2_dat %>%
#               mutate(year = "Year 2"))
# 
# # figure
# mvDatD1July %>% 
#   mutate(year = "Year 1",
#          jul_severity_cut = cut(jul_severity, breaks = seq(0, 1, by = 0.1), include.lowest = T) %>%
#            as.character()) %>%
#   group_by(jul_severity_cut) %>%
#   mutate(min_interval = parse_number(strsplit(jul_severity_cut, ",")[[1]])[1],
#          max_interval = parse_number(strsplit(jul_severity_cut, ",")[[1]])[2],
#          jul_severity_med = (max_interval + min_interval) / 2) %>%
#   full_join(mvDatD2July %>%
#               mutate(year = "Year 2",
#                      jul_severity_cut = cut(jul_severity, breaks = seq(0, 1, by = 0.1), include.lowest = T) %>%
#                        as.character()) %>%
#               group_by(jul_severity_cut) %>%
#               mutate(min_interval = parse_number(strsplit(jul_severity_cut, ",")[[1]])[1],
#                      max_interval = parse_number(strsplit(jul_severity_cut, ",")[[1]])[2],
#                      jul_severity_med = (max_interval + min_interval) / 2)) %>%
#   ggplot(aes(x = jul_severity_med, y = survival, color = year)) +
#   # geom_point(alpha = 0.5) +
#   stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
#   stat_summary(geom = "point", fun.data = mean_se, size = 2) +
#   geom_line(data = mv_surv_jul_dat, aes(x = jul_severity, linetype = year)) +
#   xlab("July disease severity") +
#   ylab("Mv survival") +
#   scale_color_manual(values = col_pal) +
#   scale_linetype_manual(values = line_pal2) +
#   fig_theme +
#   theme(legend.position = c(0.85, 0.85),
#         legend.title = element_blank())


#### fig: severity and Ev S survival, July models ####

# # models
# summary(evS_july_y1_mod3, intercepts = T)
# evS_surv_jul_y1_mod <- glmer(survival ~ log_growth + jul_severity_t + Mv_seedling_density + (1|plot_f), data = evSDatD1July, family = binomial)
# summary(evS_surv_jul_y1_mod)
# 
# summary(evS_july_y2_mod2, intercepts = T)
# evS_surv_jul_y2_mod <- lmer(survival_t ~ log_biomass + jul_severity_t + Ev_seedling_density + fungicide + (1|site/plot_f), data = evSDatD2July)
# summary(evS_surv_jul_y2_mod)
# 
# # predictions
# evS_surv_jul_y1_dat <- tibble(jul_severity = seq(min(evSDatD1July$jul_severity), max(evSDatD1July$jul_severity), length.out = 100)) %>%
#   mutate(log_growth = mean(evSDatD1July$log_growth),
#          Mv_seedling_density = mean(evSDatD1July$Mv_seedling_density),
#          jul_severity_t = logit(jul_severity, adjust = 0.001)) %>%
#   mutate(survival = predict(evS_surv_jul_y1_mod, newdata = ., type = "response", re.form = NA))
# 
# evS_surv_jul_y2_dat <- tibble(jul_severity = seq(min(evSDatD2July$jul_severity), max(evSDatD2July$jul_severity), length.out = 100)) %>%
#   mutate(log_biomass = mean(evSDatD2July$log_biomass),
#          Ev_seedling_density = mean(evSDatD2July$Ev_seedling_density),
#          fungicide = 0,
#          jul_severity_t = logit(jul_severity, adjust = 0.001)) %>%
#   mutate(survival = predict(evS_surv_jul_y2_mod, newdata = ., re.form = NA) %>%
#            inv.logit(., a = 0.001))
# 
# evS_surv_jul_dat <- evS_surv_jul_y1_dat %>%
#   mutate(year = "Year 1") %>%
#   full_join(evS_surv_jul_y2_dat %>%
#               mutate(year = "Year 2"))
# 
# # figure
# evSDatD1July %>% 
#   mutate(year = "Year 1") %>%
#   full_join(evSDatD2July %>%
#               mutate(year = "Year 2")) %>%
#   ggplot(aes(x = jul_severity, y = survival, color = year)) +
#   geom_point(alpha = 0.5) + 
#   geom_line(data = evS_surv_jul_dat, aes(linetype = year)) +
#   xlab("July disease severity") +
#   ylab("Ev seedling survival") +
#   scale_color_manual(values = col_pal) +
#   scale_linetype_manual(values = line_pal1) +
#   fig_theme +
#   theme(legend.position = c(0.85, 0.35),
#         legend.title = element_blank())


#### fig: fungicide and severity, August models ####

# models
summary(mv_aug_y1_mod2)
summary(mv_aug_y2_mod2)
summary(evS_aug_y1_mod3)
summary(evS_aug_y2_mod2)
summary(evA_aug_y1_mod2)
summary(evA_aug_y2_mod1)

# significance table
sev_fung_sig <- tibble(plant_type = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 2),
                       year = rep(c("Year 1", "Year 2"), 3),
                       sig = c("y", "y", "n", "y", "n", "n")) %>%
  mutate(plant_type =  fct_relevel(plant_type, "Mv seedling", "Ev seedling"))

# figure
sev_fung_fig <- datD1Aug %>%
  mutate(year = "Year 1",
         aug_severity = late_aug_severity_adj) %>%
  full_join(datD2Aug %>%
              mutate(year = "Year 2",
                     aug_severity = early_aug_severity_adj)) %>%
  mutate(plant_type = paste(sp, age, sep = " ") %>%
           fct_relevel("Mv seedling", "Ev seedling"),
         treatment = recode(treatment, water = "water (control)") %>%
           fct_rev()) %>%
  left_join(sev_fung_sig) %>%
  ggplot(aes(x = treatment, y = aug_severity, color = year, group = interaction(year, plant_type))) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.1), width = 0) +
  stat_summary(geom = "line", fun.data = mean_se, position = position_dodge(0.1), aes(linetype = sig)) +
  stat_summary(geom = "point", fun.data = mean_se, position = position_dodge(0.1), size = 2, aes(shape = plant_type)) +
  # annotate(geom = "text",
  #          label = expression(paste("Mv seedling Year 1: ", italic(P), " < 0.001", sep = "")),
  #          x = 0.7, y = 0.62, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Mv seedling Year 2: ", italic(P), " = 0.002", sep = "")),
  #          x = 0.7, y = 0.6, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Ev seedling Year 1: ", italic(P), " = 0.488", sep = "")),
  #          x = 0.697, y = 0.58, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Ev seedling Year 1: ", italic(P), " = 0.030", sep = "")),
  #          x = 0.697, y = 0.56, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Ev adult Year 1: ", italic(P), " = 0.204", sep = "")),
  #          x = 0.67, y = 0.54, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Ev adult Year 2: ", italic(P), " = 0.055", sep = "")),
  #          x = 0.67, y = 0.52, size = 2.5) +
  xlab("Treatment") +
  ylab("August disease severity") +
  ggtitle("F") +
  scale_color_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  scale_linetype_manual(values = line_pal2) +
  fig_theme +
  theme(legend.position = "none",
        legend.title = element_blank())


#### fig: fungicide and Ev S survival, July models ####

# models
summary(evS_july_y1_mod3, intercepts = T)
summary(evS_july_y2_mod2, intercepts = T)

# figure
evS_surv_fung_fig <- evSDatD1July %>% 
  mutate(year = "Year 1",
         treatment = recode(treatment, water = "water (control)") %>%
           fct_rev()) %>%
  full_join(evSDatD2July %>%
              mutate(year = "Year 2",
                     treatment = recode(treatment, water = "water (control)") %>%
                       fct_rev())) %>%
  ggplot(aes(x = treatment, y = survival, color = year, group = year)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun.data = mean_se, aes(linetype = year)) +
  stat_summary(geom = "point", fun.data = mean_se, size = 2, shape = shape_pal[2]) +
  # annotate(geom = "text",
  #          label = expression(paste("Year 1: ", italic(P), " = 0.468", sep = "")),
  #          x = 0.6, y = 1, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Year 2: ", italic(P), " = 0.001", sep = "")),
  #          x = 0.6, y = 0.99, size = 2.5) +
  xlab("Treatment") +
  ylab("Ev seedling survival") +
  ggtitle("G") +
  scale_color_manual(values = col_pal) +
  scale_linetype_manual(values = line_pal2) +
  fig_theme +
  theme(legend.position = "none",
        legend.title = element_blank())


#### fig: fungicide and Ev S growth, July models ####

# models
summary(evS_july_y1_mod3, intercepts = T)
summary(evS_july_y2_mod2, intercepts = T)

# figure
evS_growth_fung_fig <- evSDatD1July %>% 
  mutate(year = "Year 1",
         treatment = recode(treatment, water = "water (control)") %>%
           fct_rev()) %>%
  full_join(evSDatD2July %>%
              mutate(year = "Year 2",
                     treatment = recode(treatment, water = "water (control)") %>%
                       fct_rev(),
                     log_growth = log_biomass)) %>%
  ggplot(aes(x = treatment, y = log_growth, color = year, group = year)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.1), width = 0) +
  stat_summary(geom = "line", fun.data = mean_se, position = position_dodge(0.1)) +
  stat_summary(geom = "point", fun.data = mean_se, position = position_dodge(0.1), size = 2, shape = shape_pal[2]) +
  # annotate(geom = "text",
  #          label = expression(paste("Year 1: ", italic(P), " = 0.008", sep = "")),
  #          x = 0.6, y = 0.57, size = 2.5) +
  # annotate(geom = "text",
  #          label = expression(paste("Year 2: ", italic(P), " = 0.020", sep = "")),
  #          x = 0.6, y = 0.555, size = 2.5) +
  xlab("Treatment") +
  ylab("log(Ev seedling growth)") +
  ggtitle("H") +
  scale_color_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = "none",
        legend.title = element_blank())


#### fig: severity and Ev A growth, July models ####

# # the correlations are pretty weak and may be driven by a few outliers
# 
# # models
# summary(evA_july_y1_mod2, intercepts = T)
# evADatD1July$growth_resid <- resid(evA_growth_jul_y1_mod) %>%
#   scale()
# # correlation coefficient: -0.2444
# evADatD1July$sev_resid <- resid(lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatD1July)) %>%
#   scale()
# summary(lm(growth_resid ~ sev_resid, data = evADatD1July))
# 
# summary(evA_july_y2_mod1, intercepts = T)
# # correlation coefficient: 0.2107
# evADatD2July$growth_resid <- resid(evA_growth_jul_y2_mod) %>%
#   scale()
# evADatD2July$sev_resid <- resid(lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site), data = evADatD2July)) %>%
#   scale()
# summary(lm(growth_resid ~ sev_resid, data = evADatD2July))
# 
# # figure
# evADatD1July %>% 
#   mutate(year = "Year 1") %>%
#   full_join(evADatD2July %>%
#               mutate(year = "Year 2")) %>%
#   ggplot(aes(x = sev_resid, y = growth_resid, color = year)) +
#   geom_point(alpha = 0.5) + 
#   geom_smooth(method = "lm", se = F) +
#   xlab("July disease severity residuals") +
#   ylab("log(Ev adult growth) residuals") +
#   scale_color_manual(values = col_pal) +
#   fig_theme +
#   theme(legend.position = "none",
#         legend.title = element_blank())


#### fig: Ev density and Ev S severity, July Year 1 model ####

# model
summary(evS_july_y1_mod3, intercepts = T)
evS_sev_jul_y1_mod <- lmer(jul_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|site/plot_f), data = evSDatD1July)
summary(evS_sev_jul_y1_mod)

# predictions
evS_sev_jul_y1_evS_dat <- tibble(Ev_seedling_density = seq(min(evSDatD1July$Ev_seedling_density), max(evSDatD1July$Ev_seedling_density), length.out = 100)) %>%
  mutate(fungicide = 0,
         Mv_seedling_density = mean(evSDatD1July$Mv_seedling_density),
         Ev_adult_density = mean(evSDatD1July$Ev_adult_density)) %>%
  mutate(jul_severity = predict(evS_sev_jul_y1_mod, newdata = ., re.form = NA) %>%
           inv.logit(., a = 0.001))

evS_sev_jul_y1_evA_dat <- tibble(Ev_adult_density = seq(min(evSDatD1July$Ev_adult_density), max(evSDatD1July$Ev_adult_density), length.out = 100)) %>%
  mutate(fungicide = 0,
         Mv_seedling_density = mean(evSDatD1July$Mv_seedling_density),
         Ev_seedling_density = mean(evSDatD1July$Ev_seedling_density)) %>%
  mutate(jul_severity = predict(evS_sev_jul_y1_mod, newdata = ., re.form = NA) %>%
           inv.logit(., a = 0.001))

evS_sev_jul_y1_dat <- evS_sev_jul_y1_evS_dat %>%
  mutate(Density = "seedling",
         Ev_density = Ev_seedling_density) %>%
  full_join(evS_sev_jul_y1_evA_dat %>%
              mutate(Density = "adult",
                     Ev_density = Ev_adult_density))

# figure
evS_sev_dens_fig <- evSDatD1July %>% 
  mutate(Density = "seedling",
         Ev_density = Ev_seedling_density) %>%
  full_join(evSDatD1July %>%
              mutate(Density = "adult",
                     Ev_density = Ev_adult_density)) %>%
  mutate(Density = fct_relevel(Density, "seedling")) %>%
  ggplot(aes(x = Ev_density, y = jul_severity, shape = Density, color = Density)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.1), width = 0) +
  stat_summary(geom = "point", fun.data = mean_se, position = position_dodge(0.1), size = 2) +
  geom_line(data = evS_sev_jul_y1_dat) +
  xlab("Ev density") +
  ylab("Ev seedling July disease severity") +
  ggtitle("I") +
  scale_color_manual(values = col_pal2) +
  scale_shape_manual(values = shape_pal[2:3]) +
  fig_theme +
  theme(legend.position = c(0.77, 0.85),
        legend.title = element_blank())


#### fig: plot edge and Ev August severity, disease ####

# did with Mv July severity correlation, but the severity values are so low, it's difficult to see anything in the figure

# model
summary(dis_y1_mod2, intercepts = T)
ev_dis_aug_y1_mod <- lm(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + fungicide, data = disD1Dat)
summary(ev_dis_aug_y1_mod)
summary(dis_y2_mod2, intercepts = T)
ev_dis_aug_y2_mod <- lmer(Ev_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat)
summary(ev_dis_aug_y2_mod)
# cor.test(~ Mv_jul_severity_t + edge_severity_t, data = disD2Dat)
# mv_dis_jul_y2_int <- -5.2667
# mv_dis_jul_y2_coef <- 0.4499 # from SEM

# predictions
ev_dis_aug_y1_dat <- tibble(edge_severity = seq(min(disD1Dat$mv_inf_jul.prop), max(disD1Dat$mv_inf_jul.prop), length.out = 100)) %>%
  mutate(edge_severity_t = logit(edge_severity, adjust = 0.001),
         Mv_seedling_density = mean(disD1Dat$Mv_seedling_density),
         Ev_seedling_density = mean(disD1Dat$Ev_seedling_density),
         Ev_adult_density = mean(disD1Dat$Ev_adult_density),
         Ev_jul_severity_t = mean(disD1Dat$Ev_jul_severity_t),
         Mv_jul_severity_t = mean(disD1Dat$Mv_jul_severity_t),
         fungicide = 0,
         year = "Year 1") %>%
  mutate(severity = predict(ev_dis_aug_y1_mod, newdata = ., re.form = NA) %>%
           inv.logit(., a = 0.001))

ev_dis_aug_y2_dat <- tibble(edge_severity = seq(min(disD2Dat$edge_severity), max(disD2Dat$edge_severity), length.out = 100)) %>%
  mutate(edge_severity_t = logit(edge_severity, adjust = 0.001),
         Mv_seedling_density = mean(disD2Dat$Mv_seedling_density),
         Ev_seedling_density = mean(disD2Dat$Ev_seedling_density),
         Ev_adult_density = mean(disD2Dat$Ev_adult_density),
         Ev_jul_severity_t = mean(disD2Dat$Ev_jul_severity_t),
         Mv_jul_severity_t = mean(disD2Dat$Mv_jul_severity_t),
         log_dew = mean(disD2Dat$log_dew),
         year = "Year 2") %>%
  mutate(severity = predict(ev_dis_aug_y2_mod, newdata = ., re.form = NA) %>%
           inv.logit(., a = 0.001))

# mv_dis_jul_y2_dat <- tibble(edge_severity = seq(min(disD2Dat$edge_severity), max(disD2Dat$edge_severity), length.out = 100)) %>%
#   mutate(edge_severity_t = logit(edge_severity, adjust = 0.001),
#          Mv_jul_severity_t = mv_dis_jul_y2_int + mv_dis_jul_y2_coef * edge_severity_t,
#          severity = inv.logit(Mv_jul_severity_t, a = 0.001),
#          Species_month = "Mv July")

ev_dis_aug_dat <- full_join(ev_dis_aug_y1_dat, ev_dis_aug_y2_dat)

# figure
ev_sev_edge_fig <- disD1Dat %>% 
  mutate(year = "Year 1",
         severity = Ev_late_aug_severity,
         edge_severity = mv_inf_jul.prop) %>%
  full_join(disD2Dat %>% 
              mutate(year = "Year 2",
                     severity = Ev_early_aug_severity)) %>%
  ggplot(aes(x = edge_severity, y = severity, color = year)) +
  geom_point(alpha = 0.5, shape = shape_pal[2]) + 
  geom_line(data = ev_dis_aug_dat) +
  xlab("Edge Mv July disease severity") +
  ylab("Ev August disease severity") +
  ggtitle("J") +
  scale_color_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = c(0.78, 0.83),
        legend.title = element_blank())


#### fig: dew and Mv August severity, disease ####

# model
summary(dis_y2_mod2, intercepts = T)
mv_dis_aug_y2_mod <- lmer(Mv_aug_severity_t ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + Ev_jul_severity_t + Mv_jul_severity_t + edge_severity_t + log_dew + (1 | site), data = disD2Dat)
summary(mv_dis_aug_y2_mod)

# predictions
mv_dis_aug_y2_dat <- tibble(dew_intensity = seq(min(disD2Dat$dew_intensity), max(disD2Dat$dew_intensity), length.out = 100)) %>%
  mutate(log_dew = log(dew_intensity),
         edge_severity_t = mean(disD2Dat$edge_severity_t),
         Mv_seedling_density = mean(disD2Dat$Mv_seedling_density),
         Ev_seedling_density = mean(disD2Dat$Ev_seedling_density),
         Ev_adult_density = mean(disD2Dat$Ev_adult_density),
         Ev_jul_severity_t = mean(disD2Dat$Ev_jul_severity_t),
         Mv_jul_severity_t = mean(disD2Dat$Mv_jul_severity_t)) %>%
  mutate(severity = predict(mv_dis_aug_y2_mod, newdata = ., re.form = NA) %>%
           inv.logit(., a = 0.001))

# figure
mv_sev_dew_fig <- disD2Dat %>% 
  mutate(severity = Mv_early_aug_severity) %>%
  ggplot(aes(x = dew_intensity, y = severity)) +
  geom_point(color = col_pal[2]) + 
  geom_line(data = mv_dis_aug_y2_dat, color = col_pal[2]) +
  xlab("Dew intensity") +
  ylab("Mv August disease severity") +
  ggtitle("K") +
  fig_theme


#### combine figures ####

sem_fig <- plot_grid(mv_seeds_dens_fig,
                     ev_seeds_dens_fig,
                     ev_growth_dens_fig,
                     sev_fung_fig,
                     evS_surv_fung_fig,
                     evS_growth_fung_fig,
                     evS_sev_dens_fig,
                     ev_sev_edge_fig,
                     mv_sev_dew_fig,
                     nrow = 3)

pdf("output/piecewise_sem_2018_2019_fig.pdf", width = 7, height = 7)
sem_fig
dev.off()

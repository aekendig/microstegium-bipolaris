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
survD2Dat <- survD2Dat %>%
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
  nrow()
# replaced 8 missing values

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

sum(is.na(plotSevD2Dat$Ev_early_aug_severity)) # 3
sum(is.na(plotSevD2Dat$Ev_early_aug_severity_adj)) # 1
sum(is.na(plotSevD2Dat$Mv_early_aug_severity)) # 0
sum(is.na(plotSevD2Dat$Mv_early_aug_severity_adj)) # 0

# edge severity
edgeSevD2Dat <- mvEdgeSevD2Dat %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity) %>%
  filter(month == "jul")

# seeds
evSeedD2Datb <- evSeedD2Dat %>%
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
  mutate(asr_Ev_aug_severity = asin(sqrt(Ev_early_aug_severity)),
         asr_Mv_aug_severity = asin(sqrt(Mv_early_aug_severity)),
         asr_Ev_jul_severity = asin(sqrt(Ev_jul_severity)),
         asr_Mv_jul_severity = asin(sqrt(Mv_jul_severity)),
         log_dew = log(dew_intensity),
         asr_edge_severity = asin(sqrt(edge_severity))) %>%
  drop_na()

mvDat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  rename("ID" = "plant") %>%
  mutate(ID = as.character(ID)) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g) %>%
              rename("ID" = "plant") %>%
              mutate(ID = as.character(ID))) %>%
  full_join(survD2Dat %>%
              filter(sp == "Mv")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Mv")) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         log_biomass = log(biomass_weight.g),
         asr_severity = asin(sqrt(early_aug_severity_adj)),
         asr_survival = asin(sqrt(survival)))

evDat <- evSeedD2Datb %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD2Dat %>%
              filter(sp == "Ev")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Ev")) %>%
  full_join(plotDens) %>%
  rename("biomass_weight.g" = "weight") %>%
  mutate(log_seeds = log(seeds + 1),
         log_biomass = log(biomass_weight.g),
         asr_severity = asin(sqrt(early_aug_severity_adj)),
         asr_survival = asin(sqrt(survival)))

# combine
datSeeds <- mvDat %>%
  select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, asr_severity, asr_survival) %>%
  full_join(evDat %>%
              select(site, plot, treatment, sp, ID, age, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, log_seeds, log_biomass, asr_severity, asr_survival)) %>%
  mutate(plant_type = paste(sp, age, sep = "_"),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = "")) %>%
  drop_na()


#### visualizations ####

# correlations in disease data
disDat %>%
  select(asr_Ev_jul_severity, asr_Mv_jul_severity, asr_Ev_aug_severity, asr_Mv_aug_severity, log_dew, hrs_10_35, asr_edge_severity) %>%
  ggpairs()
# Ev and Mv July severity correlated
# hrs 10 to 35 very limited in values

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
# Ev produce more seeds at D3, Mv at D1

ggplot(datSeeds, aes(x = site, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# same as seeds

ggplot(datSeeds, aes(x = site, y = asr_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev seedlings have higher disease at D3 and Mv the highest at D4

ggplot(datSeeds, aes(x = site, y = asr_survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

# plot effects
ggplot(datSeeds, aes(x = plot_f, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev produce more seeds at D3, Mv at D1

ggplot(datSeeds, aes(x = plot_f, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# same as seeds

ggplot(datSeeds, aes(x = plot_f, y = asr_severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

ggplot(datSeeds, aes(x = plot_f, y = asr_survival)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# adults have very little variation

# distributions
ggplot(datSeeds, aes(x = log_seeds)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")
# could make seeds a binary variable for Ev seedlings

ggplot(datSeeds, aes(x = log_biomass)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")

ggplot(datSeeds, aes(x = asr_severity)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")


#### fit disease model ####

# check random effects
dis_lme1 <- lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat)
summary(dis_lme1)
dis_lme2 <- lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat)
summary(dis_lme2)
evsev_lme1 <- lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat)
summary(evsev_lme1)
mvsev_lme1 <- lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat)
summary(mvsev_lme1)
dew_lme1 <- lme(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat)
summary(dew_lme1)

# initial fit
dis_mod1 <- psem(
  lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity
)

# summary
summary(dis_mod1)
# edge severity predicts Mv july severity (add correlation)
# random effects improve R-squared values
# good fit (P = 0.216)

# update model
dis_mod2 <- psem(
  lme(asr_Mv_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_aug_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + asr_Ev_jul_severity + asr_Mv_jul_severity + asr_edge_severity + log_dew, random = ~ 1 | site, data = disDat),
  lme(asr_Ev_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(asr_Mv_jul_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  lme(log_dew ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | site, data = disDat),
  asr_Ev_jul_severity %~~% asr_Mv_jul_severity,
  asr_Ev_aug_severity %~~% asr_Mv_aug_severity,
  asr_edge_severity %~~% asr_Mv_jul_severity
)

# summary
summary(dis_mod2)
# no missing links
# good model fit (P = 0.653)
# dew intensity increase August Mv severity
# July Mv severity decreased August Ev severity, but edge Mv severity increased it
# correlations significant and positive


#### fit seeds model ####

# abandoned this approach - poor model fit

# check random effects
seeds_lme1 <- lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = datSeeds)
summary(seeds_lme1)
surv_lme1 <- lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = datSeeds)
summary(seeds_lme1)
biomass_lme1 <- lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = datSeeds)
summary(biomass_lme1)
# very small random effect sd
severity_lme1 <- lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = datSeeds)
summary(severity_lme1)

# initial fit
seeds_mod1 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = datSeeds),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = datSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = datSeeds),
  log_biomass%~~%asr_severity
)

# model summary
summary(seeds_mod1)
# small increase in R-squared with log_seeds and none with biomass
# d sep test: fungicide direct effect on seeds
# Fisher's C: P = 0.03, not great fit

# update model
seeds_mod2 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity + fungicide, random = ~ 1 | plot_f, data = datSeeds),
  lm(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = datSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = datSeeds),
  log_biomass%~~%asr_severity
)

# summary
summary(seeds_mod2)
# all omitted paths are supported
# better fit by Fisher's C
# random effects improve the R-squared values

# multi-group model
seeds_mod3 <- multigroup(seeds_mod2, group = "plant_type")

# summary
seeds_mod3
# poor model fit
# the negative effect of severity on seeds could be because Ev have higher severity in August, but produce fewer seeds


#### Mv seeds models ####

# Mv dat
mvDatSeeds <- datSeeds %>%
  filter(sp == "Mv")

# check random effects
mv_seeds_lme1 <- lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_seeds_lme1)
mv_surv_lme1 <- lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_surv_lme1)
mv_biomass_lme1 <- lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_biomass_lme1)
mv_severity_lme1 <- lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds)
summary(mv_severity_lme1)

# initial fit
mv_seeds_mod1 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds),
  log_biomass%~~%asr_severity
)

# model summary
summary(mv_seeds_mod1)
# better R-squared with random effects
# d sep test: Mv density effect on seeds
# Fisher's C: P = 0.079, okay fit
# correlation not significant

# update model
mv_seeds_mod2 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity + Mv_seedling_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = mvDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = mvDatSeeds),
  log_biomass%~~%asr_severity
)

# summary
summary(mv_seeds_mod2)
# all omitted paths are supported
# better fit by Fisher's C (P = 0.462)
# random effects improve the R-squared values
# biomass increases seeds
# Mv density reduces seeds
# fungicide reduces disease
# severity reduces survival


#### Ev seedling seeds models ####

# Ev seedling dat
evSDatSeeds <- datSeeds %>%
  filter(sp == "Ev" & age == "seedling")

# check random effects
evS_seeds_lme1 <- lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_seeds_lme1)
evS_surv_lme1 <- lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_surv_lme1)
evS_biomass_lme1 <- lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_biomass_lme1)
evS_severity_lme1 <- lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds)
summary(evS_severity_lme1)

# initial fit
evS_seeds_mod1 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  log_biomass%~~%asr_severity
)

# model summary
summary(evS_seeds_mod1)
# better R-squared with random effects
# d sep test: fungicide directly affects biomass and survival, seeds related to survival
# Fisher's C: P = 0.003, poor fit

evS_seeds_mod2 <- psem(
  lme(log_seeds ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_survival ~ log_biomass + asr_severity + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds),
  log_biomass%~~%asr_severity,
  log_seeds%~~%asr_survival
)

# model summary
summary(evS_seeds_mod2)
# better R-squared with random effects
# d sep test: all omitted paths supported
# Fisher's C: P = 0.3, good fit
# biomass increases seeds
# severity increases seeds
# fungicide increases survival and biomass
# Mv density decreases biomass
# fungicide decreases severity
# survival and seeds negatively correlated

# binomial model

# add column to data
evSDatSeeds2 <- evSDatSeeds %>%
  mutate(seeds_bin = as.numeric(log_seeds > 0))

# check model
evS_seeds_lme2 <- glmer(seeds_bin ~ log_biomass + asr_severity + (1|plot_f), family = binomial, data = evSDatSeeds2)
summary(evS_seeds_lme2)

# initial fit
evS_bin_mod1 <- psem(
  glmer(seeds_bin ~ log_biomass + asr_severity + (1|plot_f), family = binomial, data = evSDatSeeds2),
  lme(asr_survival ~ log_biomass + asr_severity, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  log_biomass%~~%asr_severity
)

# model summary
summary(evS_bin_mod1)
# better R-squared with random effects
# d sep test: fungicide directly affects biomass and survival, seeds related to survival
# Fisher's C: P = 0.01, poor fit

# update model
evS_bin_mod2 <- psem(
  glmer(seeds_bin ~ log_biomass + asr_severity + (1|plot_f), family = binomial, data = evSDatSeeds2),
  lme(asr_survival ~ log_biomass + asr_severity + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  lme(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, random = ~ 1 | plot_f, data = evSDatSeeds2),
  log_biomass%~~%asr_severity,
  asr_survival%~~%seeds_bin
)

# model summary
summary(evS_bin_mod2)
# better R-squared with random effects
# d sep test: fungicide directly affects biomass and survival, seeds related to survival
# Fisher's C: P = 0.404, good fit
# biomass increases seeds
# fungicide increases survival and biomass
# Mv density decreases biomass
# fungicide decreases severity
# survival and seeds negatively correlated
# different from continuous seeds: severity doesn't affect seeds


#### Ev adult seeds models ####

# Ev adult dat
evADatSeeds <- datSeeds %>%
  filter(sp == "Ev" & age == "adult")

# initial fit
evA_seeds_mod1 <- psem(
  lm(log_seeds ~ log_biomass + asr_severity, data = evADatSeeds),
  lm(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatSeeds),
  lm(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatSeeds),
  log_biomass%~~%asr_severity
)

# model summary
summary(evA_seeds_mod1)
# d sep test: no missing links
# Fisher's C: P = 0.394, good fit
# biomass increases seeds
# biomass and severity are positively correlated

# binomial model

# add column to data
evADatSeeds2 <- evADatSeeds %>%
  mutate(seeds_bin = as.numeric(log_seeds > 0))

# refit model
evA_seeds_mod2 <- psem(
  glm(seeds_bin ~ log_biomass + asr_severity, family = binomial, data = evADatSeeds2),
  lm(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density, data = evADatSeeds2),
  lm(asr_severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide, data = evADatSeeds2),
  log_biomass%~~%asr_severity
)

summary(evA_seeds_mod2)
# d sep test: no missing links
# Fisher's C: P = 0.59, good fit
# biomass increases seeds
# biomass and severity are positively correlated
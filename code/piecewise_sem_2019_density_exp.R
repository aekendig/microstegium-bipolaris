##### info ####

# file: piecewise_sem_2019_density_exp
# author: Amy Kendig
# date last edited: 2/2/21
# goal: fit piecewise SEM to focal plot-level data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(piecewiseSEM)
library(lme4)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import focal fitness data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") # mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") # ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")

# import biomass data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") # bg_biomass_data_processing_2019_density_exp.R

# import severity data
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R

# import environmental variables
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") # covariate_data_processing_2018_density_exp
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

# Mv biomass data for plot biomass
unique(mvBioD2Dat$process_notes) # all issues addressed
filter(mvBioD2Dat, is.na(biomass_weight.g))

# Ev biomass data for plot biomass
unique(evBioD2Dat$processing_notes) # all issues addressed
filter(evBioD2Dat, is.na(weight)) # one seedling missing

# background biomass data
filter(bgBioD2Dat, is.na(biomass.g)) # none missing

# plot biomass data
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
  select(-none_seedling_biomass)

# severity data
sevD2Dat2 <- sevD2Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  filter(ID %in% c("1", "2", "3", "A")) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity)

# edge severity
edgeSevD2Dat <- mvEdgeSevD2Dat %>%
              filter(month != "sep") %>%
              mutate(edge_severity = lesion_area.pix / leaf_area.pix,
                     edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
              select(month, site, plot, treatment, edge_severity) %>%
              pivot_wider(names_from = month,
                          names_glue = "{month}_edge_severity",
                          values_from = edge_severity)
  
# germination
# average across trials
mvGermD1Dat <- mvGermD1Dat1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final)) %>%
  select(site_plot, trial, seeds, germination_final) %>%
  full_join(mvGermD1Dat2 %>%
              select(site_plot, trial, seeds, germination_final)) %>%
  mutate(site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           recode("F" = "fungicide", "W" = "water"),
         site = ifelse(site == "P1", "D1", site)) %>%
  group_by(site, plot, treatment) %>%
  summarise(germination = mean(germination_final/seeds)) %>%
  ungroup()

# select data from year 2
# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
evGermD2Dat <- evGermDat %>%
  filter(seeds_planted > 0 & year == 2019) %>%
  mutate(week_4_emerg = case_when(week_4_emerg < week_3_emerg ~ week_3_emerg,
                                  TRUE ~ week_4_emerg),
         week_3_cut_tops = case_when(week_4_cut_tops > week_3_cut_tops & week_4_cut_tops <= week_2_emerg ~ week_4_cut_tops,
                                     TRUE ~ week_3_cut_tops),
         week_4_cut_tops = case_when(week_4_cut_tops > week_2_emerg ~ week_3_cut_tops,
                                     week_4_cut_tops < week_3_cut_tops ~ week_3_cut_tops,
                                     TRUE ~ week_4_cut_tops),
         week_3_new_emerg = week_3_emerg - week_3_cut_tops,
         week_4_new_emerg = week_4_emerg - week_4_cut_tops,
         emerg = week_2_emerg + week_4_new_emerg + week_4_soil_germ,
         germination = emerg/seeds_planted) %>%
  group_by(site, treatment, age) %>%
  summarise(germination = mean(germination)) %>%
  ungroup()

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
envDDat <- envD2Dat %>%
  select(site, plot, treatment, month, dew_intensity2, temp_avg) %>%
  rename(dew_intensity = dew_intensity2) %>%
  pivot_wider(values_from = c(dew_intensity, temp_avg),
              names_from = month,
              names_glue = "{month}_{.value}") %>%
  select(site, plot, treatment, early_aug_dew_intensity, late_aug_dew_intensity, early_aug_temp_avg, late_aug_temp_avg) %>%
  mutate_at(c("early_aug_dew_intensity", "late_aug_dew_intensity"), log) %>%
  full_join(envD1Dat %>%
              mutate(canopy_open.prop = 1 - canopy_cover.prop,
                     canopy_open.prop = asin(sqrt(canopy_open.prop))) %>%
              select(site, plot, treatment, soil_moisture_jun.prop, soil_moisture_oct.prop, canopy_open.prop))


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
  full_join(mvGermD1Dat) %>%
  full_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         fitness = log(seeds * germination * survival + 1),
         log_biomass = log(biomass_weight.g),
         severity = asin(sqrt(early_aug_severity)),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = ""))

evDat <- evSeedD2Datb %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD2Dat %>%
              filter(sp == "Ev")) %>%
  full_join(sevD2Dat2 %>%
              filter(sp == "Ev") %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(evGermD2Dat) %>%
  full_join(plotDens) %>%
  rename("biomass_weight.g" = "weight") %>%
  mutate(log_seeds = log(seeds + 1),
         fitness = log(seeds * germination * survival + 1),
         log_biomass = log(biomass_weight.g),
         severity = asin(sqrt(early_aug_severity)),
         fungicide = ifelse(treatment == "water", 0, 1),
         plot_f = paste(site, plot, substring(treatment, 1, 1), sep = ""))

# make wide
dat <- mvDat %>%
  full_join(evDat) %>%
  full_join(envDDat) %>%
  full_join(plotBioD2Dat) %>%
  full_join(edgeSevD2Dat) %>%
  mutate(plant_type = paste(sp, age, sep = "_"))


#### visualizations ####

# site effects
ggplot(dat, aes(x = site, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev produce more seeds at D3, Mv at D1

ggplot(dat, aes(x = site, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# same as seeds

ggplot(dat, aes(x = site, y = severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev seedlings have higher disease at D3 and Mv the highest at D4

# plot effects
ggplot(dat, aes(x = plot_f, y = log_seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# Ev produce more seeds at D3, Mv at D1

ggplot(dat, aes(x = plot_f, y = log_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)
# same as seeds

ggplot(dat, aes(x = plot_f, y = severity)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun.data = "mean_cl_boot", size = 2) +
  facet_wrap(~ plant_type)

# distributions
ggplot(dat, aes(x = log_seeds)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")
# could make seeds a binary variable for Ev seedlings

ggplot(dat, aes(x = log_biomass)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")

ggplot(dat, aes(x = severity)) +
  geom_histogram() +
  facet_wrap(~ plant_type, scales = "free")


#### fit Mv model ####

#### start here: put in more variables to explain disease? if not, fit Ev SEMs ####

# remove rows with NAs
mvDat2 <- mvDat %>%
  select(site, plot, treatment, plot_f, sp, ID, Mv_seedling_density, Ev_seedling_density, Ev_adult_density, fungicide, log_seeds, log_biomass, severity) %>% 
  drop_na()

# initial fit
mv_mod1 <- psem(
  lmer(log_seeds ~ log_biomass + severity + (1|plot_f), data = mvDat2),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDat2),
  lmer(severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = mvDat2),
  log_biomass%~~%severity
)

# basis set
basisSet(mv_mod1)
# densities may directly affect seeds
# fungicide may affect biomass
# fungicide may affect seeds

# model summary
summary(mv_mod1)
# d sep test: P > 0.05: fail to reject the hypothesis that the variables are conditionally independent
# Mv density on seeds is significant, others are not
# Fisher's C: P = 0.054, not great fit

# update model with direct density effect
mv_mod2 <- psem(
  lmer(log_seeds ~ log_biomass + severity + Mv_seedling_density + (1|plot_f), data = mvDat2),
  lmer(log_biomass ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + (1|plot_f), data = mvDat2),
  lmer(severity ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density + fungicide + (1|plot_f), data = mvDat2),
  log_biomass%~~%severity
)

# summary
summary(mv_mod2)
# all omitted paths are supported
# better fit by Fisher's C
# random effects improve the R-squared values

# see how removing variables affects the fit
mv_mod3 <- psem(
  lmer(log_seeds ~ log_biomass + severity + Mv_seedling_density + (1|plot_f), data = mvDat2),
  lmer(log_biomass ~ Mv_seedling_density + (1|plot_f), data = mvDat2),
  lmer(severity ~ Mv_seedling_density + fungicide + (1|plot_f), data = mvDat2),
  log_biomass%~~%severity
)

# summary
summary(mv_mod3)
# Mv density becomes a significant predictor of biomass
# was marginally significant before
# estimate is similar
# R-squared the same
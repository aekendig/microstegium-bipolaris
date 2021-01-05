##### info ####

# file: microstegium_sem_2019_density_exp
# author: Amy Kendig
# date last edited: 1/4/21
# goal: fit SEM to focal Mv data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lavaan)
library(GGally)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import focal fitness data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
# adj values substitute averages of other two plants in the plot when a plant is missing data
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")

# import biomass data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")

# import severity data
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
mvEdgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")

# import environmental variables
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv")
envD2Dat <- read_csv("intermediate-data/average_humidity_total_biomass_2019_density_exp.csv")
# from temp_humidity_analysis_2019_density_exp.R

# knowns function
known_fun <- function(knowns){
  out = knowns * (knowns + 1) / 2
  return(out)
}


#### edit data ####

# plant group densities
# add focals (-1 for Mv)
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 2, 
                                         TRUE ~ 2),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1)) %>%
  select(plot, treatment, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
mvFocD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3)))

# focal survival
mvGsSurvD2Dat <- survD2Dat %>%
  filter(focal == 1 & sp == "Mv") %>%
  group_by(site, plot, treatment, ID) %>%
  summarise(survival = 1/(length(unique(replace_date)) + 1)) %>%
  ungroup() %>%
  full_join(mvFocD2Dat) %>%
  mutate(survival = replace_na(survival, 1))

# Mv biomass data for plot biomass
unique(mvBioD2Dat$process_notes) # all issues addressed
filter(mvBioD2Dat, is.na(biomass_weight.g))
# because two individuals are missing, multiply the average individual weight by 3

# Ev biomass data for plot biomass
unique(evBioD2Dat$processing_notes) # all issues addressed
filter(evBioD2Dat, is.na(weight)) # one seedling missing
# because one individual is missing, multiply the average individual weight by 3

# background biomass data
filter(bgBioD2Dat, is.na(biomass.g)) # none missing

# plot biomass data
# only include two Mv to remove the focal
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         focal = 0) %>%
  filter(sp != "none") %>%
  full_join(mvBioD2Dat %>%
              group_by(site, plot, treatment, sp) %>%
              summarise(biomass.g = mean(biomass_weight.g, na.rm = T) * 2) %>%
              ungroup() %>%
              mutate(age = "seedling",
                     focal = 1)) %>%
  full_join(evBioD2Dat %>%
              mutate(age = case_when(ID == "A" ~ "adult",
                                     TRUE ~ "seedling")) %>%
              group_by(site, plot, treatment, sp, age) %>%
              summarise(biomass.g = mean(weight, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass.g = case_when(age == "seedling" ~ biomass.g * 3,
                                           TRUE ~ biomass.g),
                     focal = 1)) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(plot_biomass.g = sum(biomass.g)) %>%
  ungroup() %>%
  mutate(sp_age = paste(sp, age, sep = "_")) %>%
  pivot_wider(-c(sp, age),
              names_from = sp_age,
              names_glue = "{sp_age}_biomass",
              values_from = plot_biomass.g)

# severity data
# average severity by sp/age combo
# edge severity
sevD2Dat2 <- sevD2Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  group_by(month, site, plot, treatment, sp) %>%
  summarise(severity = mean(severity, na.rm = T)) %>%
  ungroup() %>%
  mutate(severity = ifelse(severity > 1, 1, severity),
         severity = asin(sqrt(severity))) %>%
  pivot_wider(names_from = c(month, sp),
              names_glue = "{month}_{sp}_severity",
              values_from = severity) %>%
  full_join(mvEdgeSevD2Dat %>%
              filter(month != "sep") %>%
              mutate(edge_severity = lesion_area.pix / leaf_area.pix,
                     edge_severity = ifelse(edge_severity > 1, 1, edge_severity),
                     edge_severity = asin(sqrt(edge_severity))) %>%
              select(month, site, plot, treatment, edge_severity) %>%
              pivot_wider(names_from = month,
                          names_glue = "{month}_edge_severity",
                          values_from = edge_severity))
  
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
         site = ifelse(site == "P1", "D1", site),
         germination = germination_final) %>%
  group_by(site, plot, treatment) %>%
  summarise(germination = mean(germination_final/seeds)) %>%
  ungroup()

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
unique(mvBioD2Dat$process_notes)

# missing data?
filter(mvSeedD2Dat, is.na(seeds)) %>%
  data.frame()
# 3 plants
filter(mvBioD2Dat, is.na(biomass_weight.g)) %>%
  data.frame()
# 2 plants

# combine data
# remove focal biomass from plot biomass and calculate ratio
mvDat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, flower_seeds, stem_seeds, seeds) %>%
  rename("ID" = "plant") %>%
  mutate(ID = as.character(ID)) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g) %>%
              rename("ID" = "plant") %>%
              mutate(ID = as.character(ID))) %>%
  full_join(mvGsSurvD2Dat) %>%
  full_join(mvGermD1Dat) %>%
  full_join(plotBioD2Dat) %>%
  full_join(sevD2Dat2) %>%
  full_join(envDDat) %>%
  full_join(plotDens) %>%
  mutate(log_flower_seeds = log(flower_seeds + 1),
         log_stem_seeds =  log(stem_seeds),
         log_seeds = log(seeds),
         log_biomass = log(biomass_weight.g),
         log_Mv_biomass = log(Mv_seedling_biomass), 
         log_Ev_seedling_biomass = log(Ev_seedling_biomass),
         log_Ev_adult_biomass = log(Ev_adult_biomass),
         log_Ev_biomass = log(Ev_seedling_biomass + Ev_adult_biomass),
         log_plot_biomass = log(Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass),
         fungicide = ifelse(treatment == "water", 0, 1),
         log_Mv_density = log(Mv_seedling_density),
         log_Ev_seedling_density = log(Ev_seedling_density),
         log_Ev_adult_density = log(Ev_adult_density))


#### visualizations ####

# fitness metrics
ggpairs(mvDat %>%
          select(log_stem_seeds, log_flower_seeds, log_biomass, survival))
# biomass and seeds are correlated, but not with survival

# plot biomass
ggpairs(mvDat %>%
          select(log_Mv_biomass, log_Ev_seedling_biomass, log_Ev_adult_biomass, log_Ev_biomass, log_plot_biomass))
# not good indicators of the same latent variable

# biomass and density
ggpairs(mvDat %>%
          select(log_Mv_biomass, log_Mv_density, Mv_seedling_density))

# severity
ggpairs(mvDat %>%
          select(jun_Mv_severity, jul_Mv_severity, early_aug_Mv_severity, late_aug_Mv_severity, may_edge_severity, jun_edge_severity, jul_edge_severity, early_aug_edge_severity))

ggpairs(mvDat %>%
          select(jun_Mv_severity, jul_Mv_severity, early_aug_Mv_severity, late_aug_Mv_severity, may_Ev_severity, jun_Ev_severity, jul_Ev_severity, early_aug_Ev_severity, late_aug_Ev_severity))


# environment
ggpairs(envDDat %>%
          select(canopy_open.prop, soil_moisture_jun.prop, soil_moisture_oct.prop, early_aug_dew_intensity, early_aug_temp_avg))


#### fit model ####

# plant fitness as a latent variable: Maddox and Antonovics 1983

# define model
mvMod1 <- '# latent variables
          fitness =~ NA*log_stem_seeds + log_flower_seeds + log_biomass + survival
          competition =~ log_plot_biomass
          disease =~ NA*jun_Mv_severity + jul_Mv_severity + early_aug_Mv_severity + late_aug_Mv_severity
          inter_disease =~ NA*may_Ev_severity + jun_Ev_severity + jul_Ev_severity + early_aug_Ev_severity + late_aug_Ev_severity
          light =~ canopy_open.prop
          water =~ w*soil_moisture_jun.prop + w*soil_moisture_oct.prop
          
          # regressions
          fitness ~ competition + disease + light + water
          competition ~ log_Mv_density + log_Ev_seedling_density + log_Ev_adult_density + disease + inter_disease + light + water
          disease ~ log_Mv_density + log_Ev_seedling_density + log_Ev_adult_density + fungicide + light + water
          inter_disease ~ log_Mv_density + log_Ev_seedling_density + log_Ev_adult_density + fungicide + light + water
          jun_Mv_severity ~ may_edge_severity
          jul_Mv_severity ~ jun_edge_severity
          early_aug_Mv_severity ~ jul_edge_severity
          late_aug_Mv_severity ~ early_aug_edge_severity
          jun_Ev_severity ~ may_edge_severity
          jul_Ev_severity ~ jun_edge_severity
          early_aug_Ev_severity ~ jul_edge_severity
          late_aug_Ev_severity ~ early_aug_edge_severity

          # correlations
          log_stem_seeds ~~ log_flower_seeds + log_biomass
          log_flower_seeds ~~ log_biomass
          disease ~~ inter_disease

          # constraints
          fitness ~~ 1*fitness
          disease ~~ 1*disease
          inter_disease ~~ 1*inter_disease'

# fit model
mvFit1 <- sem(mvMod1, data = mvDat,  
              missing = "fiml")

# examine model
varTable(mvFit1)
lavInspect(mvFit1, "theta")
summary(mvFit1, fit.measures = T, standardized = T)
modificationIndices(mvFit1, sort. = TRUE, minimum.value = 3)

#### start here: 
# figure out how to improve model fit

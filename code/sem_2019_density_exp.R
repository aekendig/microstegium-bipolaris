##### info ####

# file: sem_2019_density_exp
# author: Amy Kendig
# date last edited: 1/20/21
# goal: fit SEM to focal plot-level data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lavaan)
library(GGally)
library(lavaan.survey)
library(semPower)

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
  mutate(mv_seeds = log(seeds * germination * survival),
         mv_biomass = log(biomass_weight.g),
         mv_severity = asin(sqrt(early_aug_severity))) %>%
  select(site, plot, treatment, sp, ID, mv_seeds, mv_biomass, mv_severity)

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
  mutate(ev_seeds = log(seeds * germination * survival + 1),
         ev_biomass = log(weight),
         ev_severity = asin(sqrt(early_aug_severity))) %>%
  select(site, plot, treatment, sp, ID, ev_seeds, ev_biomass, ev_severity)

# make wide
dat <- mvDat %>%
  select(-sp) %>%
  pivot_wider(names_from = ID, 
              values_from = c(mv_seeds, mv_biomass, mv_severity),
              names_glue = "{.value}_{ID}") %>%
  full_join(evDat %>%
              select(-sp) %>%
              pivot_wider(names_from = ID, 
                          values_from = c(ev_seeds, ev_biomass, ev_severity),
                          names_glue = "{.value}_{ID}")) %>%
  full_join(envDDat) %>%
  full_join(plotDens) %>%
  full_join(plotBioD2Dat) %>%
  full_join(edgeSevD2Dat) %>%
  mutate(fungicide = ifelse(treatment == "water", 0, 1),
         log_mv_biomass = log(Mv_seedling_biomass + 1),
         log_evS_biomass = log(Ev_seedling_biomass + 1),
         log_evA_biomass = log(Ev_adult_biomass + 1))


#### visualizations ####

# histograms
ggplot(mvDat, aes(x = mv_seeds)) +
  geom_histogram()
ggplot(mvDat, aes(x = mv_biomass)) +
  geom_histogram()
ggplot(mvDat, aes(x = mv_severity)) +
  geom_histogram()

ggplot(evDat, aes(x = ev_seeds, fill = ID)) +
  geom_histogram()
ggplot(evDat, aes(x = ev_biomass)) +
  geom_histogram()
ggplot(evDat, aes(x = ev_severity)) +
  geom_histogram()

# correlations
ggpairs(dat %>%
          select(mv_seeds_1, mv_seeds_2, mv_seeds_3))
ggpairs(dat %>%
          select(mv_biomass_1, mv_biomass_2, mv_biomass_3))
ggpairs(dat %>%
          select(mv_severity_1, mv_severity_2, mv_severity_3))
ggpairs(dat %>%
          select(mv_seeds_1, mv_biomass_1, mv_severity_1))
ggpairs(dat %>%
          select(mv_seeds_2, mv_biomass_2, mv_severity_2))
ggpairs(dat %>%
          select(mv_seeds_3, mv_biomass_3, mv_severity_3))

ggpairs(dat %>%
          select(ev_seeds_1, ev_seeds_2, ev_seeds_3, ev_seeds_A))
ggpairs(dat %>%
          select(ev_biomass_1, ev_biomass_2, ev_biomass_3, ev_biomass_A))
ggpairs(dat %>%
          select(ev_severity_1, ev_severity_2, ev_severity_3, ev_severity_A))
ggpairs(dat %>%
          select(ev_seeds_1, ev_biomass_1, ev_severity_1))
ggpairs(dat %>%
          select(ev_seeds_2, ev_biomass_2, ev_severity_2))
ggpairs(dat %>%
          select(ev_seeds_3, ev_biomass_3, ev_severity_3))
ggpairs(dat %>%
          select(ev_seeds_A, ev_biomass_A, ev_severity_A))


#### fit density model ####

# define model
mod1 <- '# latent variables
          mv_resources =~ mb*mv_biomass_1 + mb*mv_biomass_2 + mb*mv_biomass_3
          mv_fitness =~ mp*mv_seeds_1 + mp*mv_seeds_2 + mp*mv_seeds_3
          mv_disease =~ md*mv_severity_1 + md*mv_severity_2 + md*mv_severity_3
          
          ev_resources =~ ebs*ev_biomass_1 + ebs*ev_biomass_2 + ebs*ev_biomass_3 + eba*ev_biomass_A
          ev_fitness =~ eps*ev_seeds_1 + eps*ev_seeds_2 + eps*ev_seeds_3 + epa*ev_seeds_A
          ev_disease =~ eds*ev_severity_1 + eds*ev_severity_2 + eds*ev_severity_3 + eda*ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          ev_biomass_A ~~ ev_seeds_A
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease
          '

# fit model
fit1 <- sem(mod1, data = dat,  
              missing = "fiml")

fit1b <- sem(mod1, data = dat)

# add site random effect
design <- svydesign(ids = ~site, nest = T,
                    data = dat)
fit.me1 <- lavaan.survey(lavaan.fit = fit1b, 
                         survey.design = design)
# error
summary(fit1b) 
# by removing missing data, there are only 53 replicates instead of 80

# summary
summary(fit1, fit.measures = T, standardized = T)

# parameters to free
modificationIndices(fit1, sort. = TRUE, minimum.value = 3)

# correlations among residuals
residuals(fit1, type = "cor")
# high correlations mean there's an issue with the model fit to the data, which is captured in the RMSEA score

# power analysis
pow1 <- semPower.postHoc(effect = 0.05, effect.measure = 'RMSEA',
                         alpha = 0.05, N = 80, df = fit1@test[[1]]$df)
summary(pow1)

# update model: remove indicator variable constraints, forces correlations among them to 1
mod2 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3 + ev_biomass_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          ev_biomass_A ~~ ev_seeds_A
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease
          '

# fit model
fit2 <- sem(mod2, data = dat,  
            missing = "fiml")
summary(fit2, fit.measures = T, standardized = T)
anova(fit1, fit2)
# improved fit, but still not great
# CFI and TLI indicate okay fit (>0.9) and RMSEA indicates okay fit (P > 0.05)

# power analysis
pow2 <- semPower.postHoc(effect = 0.05, effect.measure = 'RMSEA',
                         alpha = 0.05, N = 80, df = fit2@test[[1]]$df)
summary(pow2)

# parameters to free
modificationIndices(fit2, sort. = TRUE, minimum.value = 3)
# I don't think any make sense to add in

# links to remove?
summary(fit2, fit.measures = T, standardized = T)
# correlations:
# mv and ev fitness
# mv and ev disease (interested in reporting that though, suggests limited transmission)
# resources and disease for each (interested in reporting, part of theoretical model)
# all others significant and within-individual
# regressions
# all are part of main question

#### start here: remove non-significant links and test with anova (Grace et al. 2020) ####

#### fit biomass model ####
modb <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3 + ev_biomass_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass
          mv_disease ~ fungicide + Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass
          ev_disease ~ fungicide + Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          ev_biomass_A ~~ ev_seeds_A
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease
          '

# fit model
fitb <- sem(modb, data = dat,  
            missing = "fiml")
summary(fitb, fit.measures = T, standardized = T)
anova(fit2, fitb)
# slightly higher chi-square
# some differences in interpretation, e.g. higher ev plot biomass -> higher ev biomass
##### info ####

# file: sem_2019_density_exp
# author: Amy Kendig
# date last edited: 1/21/21
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
library(lavaanPlot)

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
  mutate(fungicide = ifelse(treatment == "water", 0, 1))


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
          ev_disease ~~ ev_resources + mv_disease'

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

# visualize
lavaanPlot(model = fit2,
           node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE,
           covs = TRUE)


#### remove links in model ####

# remove interspecific disease correlation (P = 0.997, Std.all = -0.001)
mod3 <- '# latent variables
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
          ev_disease ~~ ev_resources'
fit3 <- sem(mod3, data = dat,  
            missing = "fiml")
anova(fit2, fit3) # not sig diff
summary(fit3, fit.measures = T, standardized = T)

# remove intraspecific Ev adult competition (P = 0.937, Std.all = 0.01)
mod4 <- '# latent variables
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
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
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
          ev_disease ~~ ev_resources'
fit4 <- sem(mod4, data = dat,  
            missing = "fiml")
anova(fit3, fit4) # not sig diff
summary(fit4, fit.measures = T, standardized = T)

# remove Ev adult density from Mv disease (P = 0.799, Std.all = -0.033)
mod5 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3 + ev_biomass_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
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
          ev_disease ~~ ev_resources'
fit5 <- sem(mod5, data = dat,  
            missing = "fiml")
anova(fit4, fit5) # not sig diff
summary(fit5, fit.measures = T, standardized = T)

# remove Ev seedling density from Mv disease (P = 0.781, Std.all = -0.034)
mod6 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3 + ev_biomass_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
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
          ev_disease ~~ ev_resources'
fit6 <- sem(mod6, data = dat,  
            missing = "fiml")
anova(fit5, fit6) # not sig diff
summary(fit6, fit.measures = T, standardized = T)

# remove Ev adult density from Mv resources (P = 0.678, Std.all = 0.040)
mod7 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3 + ev_biomass_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
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
          ev_disease ~~ ev_resources'
fit7 <- sem(mod7, data = dat,  
            missing = "fiml")
anova(fit6, fit7) # not sig diff
summary(fit7, fit.measures = T, standardized = T)

# remove Ev adult biomass as an indicator for Ev resources (P = 0.670, Std.all = 0.047) and its correlation with seeds 
mod8 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit8 <- sem(mod8, data = dat,  
            missing = "fiml")
anova(fit7, fit8) # warning about comparison, but there was a big decrease in AIC, BIC, and Chisq
summary(fit8, fit.measures = T, standardized = T)

# remove Ev seedling density from Ev disease (P = 0.646, Std.all = 0.065)
mod9 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit9 <- sem(mod9, data = dat,  
            missing = "fiml")
anova(fit8, fit9) # not significantly different
summary(fit9, fit.measures = T, standardized = T)

# remove Ev adult density from Ev disease (P = 0.548, Std.all = 0.080)
mod10 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit10 <- sem(mod10, data = dat,  
            missing = "fiml")
anova(fit9, fit10) # warning because Ev adult density has been completely removed, lower AIC, BIC, and Chisq
summary(fit10, fit.measures = T, standardized = T)

# remove Ev seedling density from Mv resources (P = 0.547, Std.all = 0.069)
mod11 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide + Mv_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit11 <- sem(mod11, data = dat,  
             missing = "fiml")
anova(fit10, fit11) # no significant difference
summary(fit11, fit.measures = T, standardized = T)

# remove Mv seedling density from Mv disease (P = 0.529, Std.all = -0.073)
mod12 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit12 <- sem(mod12, data = dat,  
             missing = "fiml")
anova(fit11, fit12) # no significant difference
summary(fit12, fit.measures = T, standardized = T)

# remove Ev seedling density from Ev resources (P = 0.469, Std.all = 0.065)
mod13 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources'
fit13 <- sem(mod13, data = dat,  
             missing = "fiml")
anova(fit12, fit13) # warning because Ev seedling density removed from the model, AIC, BIC, and Chisq went down
summary(fit13, fit.measures = T, standardized = T)

# remove correlation between Mv resources and disease (P = 0.203, Std.all = 0.133)
mod14 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources'
fit14 <- sem(mod14, data = dat,  
             missing = "fiml")
anova(fit13, fit14) # no significant difference
summary(fit14, fit.measures = T, standardized = T)

# remove mv disease from mv fitness (P = 0.149, Std.all = -0.099)
mod15 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources'
fit15 <- sem(mod15, data = dat,  
             missing = "fiml")
anova(fit14, fit15) # no significant difference
summary(fit15, fit.measures = T, standardized = T)

lavaanPlot(model = fit15,
           node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE,
           covs = TRUE)

# fix the correlation between Mv disease and Ev fitness (P = 0.7, Std.all = -0.057)
mod16 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources
          mv_disease ~~ 0*ev_fitness'
fit16 <- sem(mod16, data = dat,  
             missing = "fiml")
anova(fit15, fit16) # no significant difference
summary(fit16, fit.measures = T, standardized = T)

# fix the correlation between Mv and Ev fitness (P = 0.149, Std.all = -0.23)
mod17 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources
          
          # constraints
          mv_disease ~~ 0*ev_fitness
          mv_fitness ~~ 0*ev_fitness'
fit17 <- sem(mod17, data = dat,  
             missing = "fiml")
anova(fit16, fit17) # no significant difference
summary(fit17, fit.measures = T, standardized = T)

# fix the correlation between Mv fitness and Mv disease (P = 0.15, Std.all = -0.21)
mod18 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources
          
          # constraints
          mv_disease ~~ 0*ev_fitness
          mv_fitness ~~ 0*ev_fitness + 0*mv_disease'
fit18 <- sem(mod18, data = dat,  
             missing = "fiml")
anova(fit17, fit18) # no significant difference
summary(fit18, fit.measures = T, standardized = T)

# remove the effect of Mv density of Ev resources (P = 0.135, Std.all = -0.177)
mod19 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_disease ~ fungicide + Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources
          
          # constraints
          mv_disease ~~ 0*ev_fitness
          mv_fitness ~~ 0*ev_fitness + 0*mv_disease'
fit19 <- sem(mod19, data = dat,  
             missing = "fiml")
anova(fit18, fit19) # no significant difference
summary(fit19, fit.measures = T, standardized = T)

# remove the effect of Mv density of Ev disease (P = 0.214, Std.all = -0.155)
mod20 <- '# latent variables
          mv_resources =~ mv_biomass_1 + mv_biomass_2 + mv_biomass_3
          mv_fitness =~ mv_seeds_1 + mv_seeds_2 + mv_seeds_3
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_biomass_1 + ev_biomass_2 + ev_biomass_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ Mv_seedling_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_disease ~ fungicide
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_biomass_1 ~~ mv_seeds_1
          mv_biomass_2 ~~ mv_seeds_2
          mv_biomass_3 ~~ mv_seeds_3
          
          ev_biomass_1 ~~ ev_seeds_1
          ev_biomass_2 ~~ ev_seeds_2
          ev_biomass_3 ~~ ev_seeds_3
          
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources
          
          # constraints
          mv_disease ~~ 0*ev_fitness
          mv_fitness ~~ 0*ev_fitness + 0*mv_disease'
fit20 <- sem(mod20, data = dat,  
             missing = "fiml")
anova(fit19, fit20) # no significant difference
summary(fit20, fit.measures = T, standardized = T)

# visualize
lavaanPlot(model = fit20,
           node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE,
           covs = TRUE)


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
##### info ####

# file: sem_2018_density_exp
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
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv") 
# ev_seeds_data_processing_2018.R and ev_seeds_data_processing_2019.R
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")

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
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  filter(!is.na(survival)) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(survival = mean(survival))
# 233 entries, the one missing is the Ev adult that was not virginicus

# severity data
sevD1Datb <- sevD1Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity)) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
  pivot_wider(names_from = month,
              names_glue = "{month}_severity",
              values_from = severity)

sum(!is.na(sevD1Datb$jul_severity))
sum(!is.na(sevD1Datb$late_aug_severity))
sum(!is.na(sevD1Datb$sep_severity))
  
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

# select data from year 1
# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
# use average of three sites to add in missing value
evGermD1Dat <- evGermDat %>%
  filter(seeds_planted > 0 & year == 2018) %>%
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
  ungroup() %>%
  group_by(treatment, age) %>%
  mutate(germination2 = mean(germination)) %>%
  ungroup %>%
  full_join(tibble(site = "D3",
                   treatment = "fungicide",
                   age = "seedling",
                   germination = 0.161)) %>%
  select(-germination2)

# seeds
evSeedD1Datb <- evSeedD1Dat %>%
  filter(focal == 1 & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup()

# combine data
mvDat <- mvSeedD1Dat %>%
  full_join(growthD1Dat %>%
              filter(sp == "Mv") %>%
              select(site:ID, height_growth, tiller_growth)) %>%
  full_join(survD1Datb %>%
              filter(sp == "Mv")) %>%
  full_join(sevD1Datb %>%
              filter(sp == "Mv")) %>%
  full_join(mvGermD1Dat) %>%
  rowwise() %>%
  mutate(late_aug_severity = ifelse(is.na(late_aug_severity), 
                                    mean(c(jul_severity, sep_severity), na.rm = T), 
                                    late_aug_severity)) %>%
  ungroup() %>%
  mutate(mv_seeds = log(total_seeds * germination * survival + 1),
         mv_height = height_growth,
         mv_tiller = tiller_growth,
         mv_severity = asin(sqrt(late_aug_severity))) %>%
  select(site, plot, treatment, sp, ID, mv_seeds, mv_height, mv_tiller, mv_severity)

# one plant is not Ev (remove from all: D2 7W Ev A)
# could not fit SEM with individual-level Ev severity (too many missing)
evDat <- evSeedD1Datb %>%
  full_join(growthD1Dat %>%
              filter(sp == "Ev") %>%
              select(site:ID, height_growth:basal_growth) %>%
              mutate(age = ifelse(ID == "A", "adult", "seedling"))) %>%
  full_join(survD1Datb %>%
              filter(sp == "Ev")) %>%
  full_join(sevD1Datb %>%
              filter(sp == "Ev" &
                       !(site == "D2" & plot == 7 & treatment == "water" & ID == "A")) %>%
              rowwise() %>%
              mutate(late_aug_severity = ifelse(is.na(late_aug_severity), 
                                                mean(c(jul_severity, sep_severity), na.rm = T), 
                                                late_aug_severity)) %>%
              ungroup() %>%
              group_by(site, plot, treatment, sp) %>%
              summarise(late_aug_severity = mean(late_aug_severity, na.rm = T)) %>%
              ungroup()) %>%
  full_join(evGermD1Dat) %>%
  filter(!(site == "D2" & plot == 7 & treatment == "water" & age == "adult")) %>%
  mutate(seeds = replace_na(seeds, 0),
         ev_seeds = log(seeds * germination * survival + 1),
         ev_height = height_growth,
         ev_tiller = tiller_growth,
         ev_basal = basal_growth,
         ev_severity = asin(sqrt(late_aug_severity))) %>%
  select(site, plot, treatment, sp, ID, ev_seeds, ev_height, ev_tiller, ev_basal, ev_severity)

# make wide
dat <- mvDat %>%
  select(-sp) %>%
  pivot_wider(names_from = ID, 
              values_from = c(mv_seeds, mv_height, mv_tiller, mv_severity),
              names_glue = "{.value}_{ID}") %>%
  full_join(evDat %>%
              select(-sp) %>%
              pivot_wider(names_from = ID, 
                          values_from = c(ev_seeds, ev_height, ev_tiller, ev_basal, ev_severity),
                          names_glue = "{.value}_{ID}")) %>%
  full_join(envD1Dat) %>%
  full_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "water", 0, 1),
         log_mv_density = log(Mv_seedling_density),
         log_evS_density = log(Ev_seedling_density),
         log_evA_density = log(Ev_adult_density)) %>%
  filter(!(site == "D4" & plot %in% c(8, 10) & treatment == "fungicide"))


#### visualizations ####

# histograms
ggplot(dat, aes(x = mv_seeds_1)) +
  geom_histogram()
ggplot(mvDat, aes(x = mv_height)) +
  geom_histogram()
ggplot(mvDat, aes(x = mv_tiller)) +
  geom_histogram()
ggplot(mvDat, aes(x = mv_severity)) +
  geom_histogram()

ggplot(evDat, aes(x = ev_seeds, fill = ID)) +
  geom_histogram()
ggplot(evDat, aes(x = ev_height)) +
  geom_histogram()
ggplot(evDat, aes(x = ev_tiller)) +
  geom_histogram()
ggplot(evDat, aes(x = ev_basal)) +
  geom_histogram()
ggplot(dat, aes(x = ev_severity_1)) +
  geom_histogram()

# correlations
ggpairs(sevD1Datb %>%
          select(jul_severity, late_aug_severity, sep_severity))
ggpairs(dat %>%
          select(mv_seeds_1, mv_seeds_2, mv_seeds_3))
ggpairs(dat %>%
          select(mv_height_1, mv_height_2, mv_height_3))
ggpairs(dat %>%
          select(mv_tiller_1, mv_tiller_2, mv_tiller_3))
ggpairs(dat %>%
          select(mv_severity_1, mv_severity_2, mv_severity_3))
ggpairs(dat %>%
          select(mv_height_1, mv_tiller_1, mv_severity_1)) # tiller and height 0.24
ggpairs(dat %>%
          select(mv_height_2, mv_tiller_2, mv_severity_2))
ggpairs(dat %>%
          select(mv_height_3, mv_tiller_3, mv_severity_3))

ggpairs(dat %>%
          select(ev_seeds_1, ev_seeds_2, ev_seeds_3, ev_seeds_A))
ggpairs(dat %>%
          select(ev_height_1, ev_height_2, ev_height_3, ev_height_A))
ggpairs(dat %>%
          select(ev_tiller_1, ev_tiller_2, ev_tiller_3, ev_tiller_A))
ggpairs(dat %>%
          select(ev_basal_1, ev_basal_2, ev_basal_3, ev_basal_A))
ggpairs(dat %>%
          select(ev_severity_1, ev_severity_2, ev_severity_3, ev_severity_A))
ggpairs(dat %>%
          select(ev_seeds_1, ev_height_1, ev_tiller_1, ev_basal_1, ev_severity_1)) # height and severity 0.68, correlations among tiller, basal, and height, seeds and basal 0.35
ggpairs(dat %>%
          select(ev_seeds_2, ev_height_2, ev_tiller_2, ev_basal_2, ev_severity_2)) # correlations among tiller, basal, and height, seeds and basal 0.38
ggpairs(dat %>%
          select(ev_seeds_3, ev_height_3, ev_tiller_3, ev_basal_3, ev_severity_3)) # correlations among tiller, basal, and height, seeds and basal 0.41
ggpairs(dat %>%
          select(ev_seeds_A, ev_height_A, ev_tiller_A, ev_basal_A, ev_severity_A)) # correlations among tiller, basal, and height, seeds and basal 0.23


#### fit density model ####

# define model
mod1 <- '# latent variables
          mv_resources =~ mv_height_1 + mv_height_2 + mv_height_3 + mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_height_1 + ev_height_2 + ev_height_3 + ev_height_A + ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A + ev_basal_1 + ev_basal_2 + ev_basal_3 + ev_basal_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1 + ev_severity_2 + ev_severity_3 + ev_severity_A
          
          # regressions
          mv_resources ~ log_mv_density + log_evS_density + log_evA_density
          mv_disease ~ fungicide + log_mv_density + log_evS_density + log_evA_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ log_mv_density + log_evS_density + log_evA_density
          ev_disease ~ fungicide + log_mv_density + log_evS_density + log_evA_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease'

# fit model
fit1 <- sem(mod1, data = dat,  
            missing = "fiml")
summary(fit1, fit.measures = T, standardized = T)
# does not converge

# refit model with "size" instead of separate size measurements
# use plot-level ev severity (too much missing data)
# okay to use non-log-transformed density
mod2 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease'

# fit model
fit2 <- sem(mod2, data = dat,  
            missing = "fiml")
summary(fit2, fit.measures = T, standardized = T)

# power analysis
pow2 <- semPower.postHoc(effect = 0.05, effect.measure = 'RMSEA',
                         alpha = 0.05, N = 78, df = fit2@test[[1]]$df)
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

# remove Mv seedling density from disease (P = 0.994, Std.all = -0.001)
mod3 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Ev_seedling_density + Ev_adult_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease'

# fit model
fit3 <- sem(mod3, data = dat,  
            missing = "fiml")
anova(fit2, fit3) # not sig diff
summary(fit3, fit.measures = T, standardized = T)

# constrain correlation between Mv and Ev fitness (P = 0.968, Std.all = -0.006)
mod4 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Ev_seedling_density + Ev_adult_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          mv_fitness ~~ 0*ev_fitness'

# fit model
fit4 <- sem(mod4, data = dat,  
            missing = "fiml")
anova(fit3, fit4) # not sig diff
summary(fit4, fit.measures = T, standardized = T)

# remove Ev adult density from Mv disease (P = 0.959, Std.all = -0.005)
mod5 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Ev_seedling_density
          mv_fitness ~ mv_resources + mv_disease
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          mv_fitness ~~ 0*ev_fitness'

# fit model
fit5 <- sem(mod5, data = dat,  
            missing = "fiml")
anova(fit4, fit5) # not sig diff
summary(fit5, fit.measures = T, standardized = T)

# remove Mv disease from Mv fitness (P = 0.918, Std.all = 0.012)
mod6 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Ev_seedling_density
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          mv_fitness ~~ 0*ev_fitness'

# fit model
fit6 <- sem(mod6, data = dat,  
            missing = "fiml")
anova(fit5, fit6) # not sig diff
summary(fit6, fit.measures = T, standardized = T)

# constrain correlation between Mv disease and Ev fitness (P = 0.984, Std.all = 0.003)
mod7 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide + Ev_seedling_density
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit7 <- sem(mod7, data = dat,  
            missing = "fiml")
anova(fit6, fit7) # not sig diff
summary(fit7, fit.measures = T, standardized = T)

# remove Ev seedling density from Mv disease (P = 0.898, Std.all = 0.012)
mod8 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit8 <- sem(mod8, data = dat,  
            missing = "fiml")
anova(fit7, fit8) # not sig diff
summary(fit8, fit.measures = T, standardized = T)

# remove Mv seedling density from Mv resources (P = 0.863, Std.all = 0.028)
mod9 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_seedling_density + Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit9 <- sem(mod9, data = dat,  
            missing = "fiml")
anova(fit8, fit9) # not sig diff
summary(fit9, fit.measures = T, standardized = T)

# remove Ev seedling density from Mv resources (P = 0.761, Std.all = -0.045)
mod10 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ fungicide + Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit10 <- sem(mod10, data = dat,  
            missing = "fiml")
anova(fit9, fit10) # not sig diff
summary(fit10, fit.measures = T, standardized = T)

# remove fungicide from Ev disease (P = 0.703, Std.all = -0.042)
mod11 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ mv_disease + ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit11 <- sem(mod11, data = dat,  
             missing = "fiml")
anova(fit10, fit11) # not sig diff
summary(fit11, fit.measures = T, standardized = T)

# remove correlation between Mv resources and disease (P = 0.613, Std.all = -0.089)
mod12 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ ev_resources
          ev_disease ~~ ev_resources + mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit12 <- sem(mod12, data = dat,  
             missing = "fiml")
anova(fit11, fit12) # not sig diff
summary(fit12, fit.measures = T, standardized = T)

# remove correlation between Ev resources and disease (P = 0.470, Std.all = 0.130)
mod13 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ ev_resources
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit13 <- sem(mod13, data = dat,  
             missing = "fiml")
anova(fit12, fit13) # not sig diff
summary(fit13, fit.measures = T, standardized = T)

# remove Ev seedling density from Ev resources (P = 0.482, Std.all = 0.121)
mod14 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ ev_resources
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease'

# fit model
fit14 <- sem(mod14, data = dat,  
             missing = "fiml")
anova(fit13, fit14) # not sig diff
summary(fit14, fit.measures = T, standardized = T)

# constrain correlation between Mv fitness and disease (P = 0.444, Std.all = -0.105)
mod15 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          mv_resources ~~ ev_resources
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit15 <- sem(mod15, data = dat,  
             missing = "fiml")
anova(fit14, fit15) # not sig diff
summary(fit15, fit.measures = T, standardized = T)

# remove the corelation between Mv and Ev resources (P = 0.336, Std.all = 0.188)
mod16 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density + Ev_adult_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit16 <- sem(mod16, data = dat,  
             missing = "fiml")
anova(fit15, fit16) # not sig diff
summary(fit16, fit.measures = T, standardized = T)

# remove Ev adult density from Ev resources (P = 0.192, Std.all = 0.201)
mod17 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density + Ev_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit17 <- sem(mod17, data = dat,  
             missing = "fiml")
anova(fit16, fit17) # not sig diff
summary(fit17, fit.measures = T, standardized = T)

# remove Ev seedling density from Ev disease (P = 0.171, Std.all = 0.165)
mod18 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density + Ev_adult_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit18 <- sem(mod18, data = dat,  
             missing = "fiml")
anova(fit17, fit18) # warning because Ev seedling density was removed from the model, chisq and BIC decreased, AIC didn't
summary(fit18, fit.measures = T, standardized = T)

# remove Ev adult density from Ev disease (P = 0.170, Std.all = 0.153)
mod19 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3 + ev_seeds_A
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit19 <- sem(mod19, data = dat,  
             missing = "fiml")
anova(fit18, fit19) # not sig diff
summary(fit19, fit.measures = T, standardized = T)

# remove Ev A seeds from Ev fitness (P = 0.163, Std.all = 0.189)
mod20 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease
          
          # correlations
          ev_disease ~~ mv_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit20 <- sem(mod20, data = dat,  
             missing = "fiml")
anova(fit19, fit20) # warning about using different variables, but all three fit measures decreased
summary(fit20, fit.measures = T, standardized = T)

# remove Mv and Ev disease correlation (P = 0.132, Std.all = -0.203)
mod21 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3 + ev_tiller_A
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit21 <- sem(mod21, data = dat,  
             missing = "fiml")
anova(fit20, fit21) # not sig diff
summary(fit21, fit.measures = T, standardized = T)

# remove Ev A tiller from Ev resources (P = 0.121, Std.all = 0.331)
mod22 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_disease ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit22 <- sem(mod22, data = dat,  
             missing = "fiml")
anova(fit21, fit22) # warning about different number of variables, but all fit measures decreased
summary(fit22, fit.measures = T, standardized = T)

# remove Mv seedling density from Ev resource (P = 0.080, Std.all = -0.283)
mod23 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_disease ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit23 <- sem(mod23, data = dat,  
             missing = "fiml")
anova(fit22, fit23) # significantly worse - keep in

# remove Mv seedling density from Ev disease (P = 0.074, Std.all = 0.199)
mod24 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_fitness ~ ev_resources + ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit24 <- sem(mod24, data = dat,  
             missing = "fiml")
anova(fit22, fit24) # not sig diff
summary(fit24, fit.measures = T, standardized = T)

# remove Ev resources from Ev fitness (P = 0.058, Std.all = 0.337)
mod25 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_fitness ~ ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit25 <- sem(mod25, data = dat,  
             missing = "fiml")
anova(fit24, fit25) # warning about comparison, fit metrics all increased
summary(fit25, fit.measures = T, standardized = T)

# remove Ev resources from Ev fitness and constrain correlations that appear because of this (P = 0.058, Std.all = 0.337)
mod26 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_resources ~ Mv_seedling_density
          ev_fitness ~ ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease
          ev_resources ~~ 0*mv_fitness + 0*mv_disease + 0*ev_fitness'

# fit model
fit26 <- sem(mod26, data = dat,  
             missing = "fiml")
anova(fit24, fit26) # model with Ev resources significantly better

# check again to make sure Mv seedling density should be included (P = 0.079, Std.all = -0.283)
mod27 <- '# latent variables
          mv_resources =~ mv_tiller_1 + mv_tiller_2 + mv_tiller_3
          mv_fitness =~ mv_seeds_1
          mv_disease =~ mv_severity_1 + mv_severity_2 + mv_severity_3
          
          ev_resources =~ ev_tiller_1 + ev_tiller_2 + ev_tiller_3
          ev_fitness =~ ev_seeds_1 + ev_seeds_2 + ev_seeds_3
          ev_disease =~ ev_severity_1
          
          # regressions
          mv_resources ~ Ev_adult_density
          mv_disease ~ fungicide
          mv_fitness ~ mv_resources
          
          ev_fitness ~ ev_resources + ev_disease

          # constraints
          ev_fitness ~~ 0*mv_fitness + 0*mv_disease
          mv_fitness ~~ 0*mv_disease'

# fit model
fit27 <- sem(mod27, data = dat,  
             missing = "fiml")
anova(fit24, fit27) # warning about variables, AIC and BIC increased

# visualize
lavaanPlot(model = fit24,
           node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE,
           covs = TRUE)

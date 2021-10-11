##### info ####

# file: focal_seed_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 9/29/21
# goal: analyses of seeds as a function of density and severity

# Ricker model for density analysis:
# log(Nt+1/Nt) = r - alphaNN x Nt - alphaNM x Mt
# haven't updated this to full model, with all background and focal in same
# do the update if using seed ~ density results
# use focal_growth_2018_2019_density_exp.R for guide
# use hurdle model (as in severity) for fit


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)
library(GGally)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import seed data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
# ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R

#import severity data
sevD1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R

# import survival data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")

# import growth data
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD1Dat <- read_csv("intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")
# mv_biomass_data_processing_2018_density_exp.R
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# severity data
sevD1Dat2 <- sevD1Dat %>%
  select(month, site, plot, treatment, sp, age, severity) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")

sevD2Dat2 <- sevD2Dat %>%
  select(month, site, plot, treatment, sp, age, severity) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")

# plant group densities
plotDens <- plotsD %>%
  mutate(density = case_when(plot %in% 2:4 ~ background_density + 3,
                             plot %in% 5:7 ~ background_density + 3,
                             plot%in% 8:10 ~ background_density + 1,
                             plot == 1 ~ 0)) %>%
  select(plot, treatment, background, density_level, density)

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal)) %>%
  filter(!is.na(survival))

# 2019 list of all plants
# all dead plants were replaced
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

# Ev seeds 2018
evSeedD1Dat2 <- evSeedD1Dat %>%
  filter(focal == 1 & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev" & survival == 1) %>%
              select(-survival)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# Ev seeds 2019
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(focD2Dat %>%
              filter(sp == "Ev")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# Mv seeds 2018
# quadrat was 0.49 x 0.25 m
mvSeedD1Dat2 <- mvSeedD1Dat %>% # none missing 
  left_join(tillerD1Dat %>%
              filter(sp == "Mv") %>%
              group_by(site, plot, treatment, sp) %>%
              summarise(tillers_jul = mean(tillers_jul, na.rm = T)) %>%
              ungroup()) %>%
  left_join(plotDens) %>%
  mutate(seeds_per_plant = seeds_per_stem * tillers_jul,
         seeds_per_m2 = (seeds_bio + seeds_soil) / (0.49 * 0.25),
         age = "seedling")

# Mv seeds 2019
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# combine by year
seedD1Dat <- evSeedD1Dat2 %>%
  full_join(mvSeedD1Dat2 %>%
              rename(seeds = seeds_per_plant)) %>%
  left_join(sevD1Dat2) %>%
  left_join(tillerD1Dat %>%
              filter(sp == "Ev") %>%
              mutate(plant_growth = log(tillers_jul/tillers_jun)) %>%
              select(site, plot, treatment, sp, ID, plant_growth) %>%
              full_join(mvBioD1Dat %>%
                          mutate(sp = "Mv",
                                 plant_growth = log(bio.g * 0.25 * 0.49)) %>%
                          select(site, plot, treatment, sp, plant_growth))) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         seeds1 = seeds + 1,
         log_seeds = log(seeds1),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1)))

seedD2Dat <- evSeedD2Dat2 %>%
  full_join(mvSeedD2Dat2) %>%
  left_join(sevD2Dat2) %>%
  left_join(mvBioD2Dat %>%
              mutate(ID = as.character(plant)) %>%
              select(site, plot, treatment, sp, ID, biomass_weight.g) %>%
              full_join(evBioD2Dat %>%
                          rename(biomass_weight.g = weight) %>%
                          select(site, plot, treatment, sp, ID, biomass_weight.g))) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         seeds1 = seeds + 1,
         log_seeds = log(seeds1),
         plant_growth = log(biomass_weight.g),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1)))


#### initial severity visualizations ####

seedD1Dat %>%
  filter(foc == "s") %>%
  select(jul_severity, late_aug_severity, sep_severity, log_seeds) %>%
  ggpairs()

seedD1Dat %>%
  filter(foc == "a") %>%
  select(jul_severity, late_aug_severity, sep_severity, log_seeds) %>%
  ggpairs()

seedD1Dat %>%
  filter(foc == "m") %>%
  select(jul_severity, late_aug_severity, sep_severity, log_seeds) %>%
  ggpairs()

seedD2Dat %>%
  filter(foc == "s") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, log_seeds) %>%
  ggpairs()

seedD2Dat %>%
  filter(foc == "a") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, log_seeds) %>%
  ggpairs()

seedD2Dat %>%
  filter(foc == "m") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, log_seeds) %>%
  ggpairs()


#### fit severity models ####

# remove missing data
seedD1Dat2 <- seedD1Dat %>%
  filter(!is.na(sep_severity) & !is.na(plant_growth))

seedD2Dat2 <- seedD2Dat %>%
  filter(!is.na(late_aug_severity) & !is.na(plant_growth))

# fit models
seedSevD1Mod <- brm(data = seedD1Dat2, family = hurdle_lognormal,
                    bf(seeds ~ foc * (sep_severity + plant_growth) + (1|plotf),
                       hu ~ foc * (sep_severity + plant_growth) + (1|plotf)),
                    prior <- c(prior(normal(1, 1), class = "Intercept"),
                               prior(normal(0, 10), class = "Intercept", dpar = "hu"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedSevD1Mod)

seedSevD2Mod <- brm(data = seedD2Dat2, family = hurdle_lognormal,
                    bf(seeds ~ foc * (late_aug_severity + plant_growth) + (1|plotf),
                       hu ~ foc * (late_aug_severity + plant_growth) + (1|plotf)),
                    prior <- c(prior(normal(1, 1), class = "Intercept"),
                               prior(normal(0, 10), class = "Intercept", dpar = "hu"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedSevD2Mod)

# save models
save(seedSevD1Mod, file = "output/seed_severity_model_2018_density_exp.rda")
save(seedSevD2Mod, file = "output/seed_severity_model_2019_density_exp.rda")


#### severity coefficients ####

# 2018
evS_2018_sev = "sep_severity = 0"
mv_2018_sev = "sep_severity + focm:sep_severity = 0"
evA_2018_sev = "sep_severity + foca:sep_severity = 0"
evS_bin_2018_sev = "hu_sep_severity = 0"
mv_bin_2018_sev = "hu_sep_severity + hu_focm:sep_severity = 0"
evA_bin_2018_sev = "hu_sep_severity + hu_foca:sep_severity = 0"
evS_growth = "plant_growth = 0"
mv_growth = "plant_growth + focm:plant_growth = 0"
evA_growth = "plant_growth + foca:plant_growth = 0"
evS_bin_growth = "hu_plant_growth = 0"
mv_bin_growth = "hu_plant_growth + hu_focm:plant_growth = 0"
evA_bin_growth = "hu_plant_growth + hu_foca:plant_growth = 0"

seedSevD1Coef <- hypothesis(seedSevD1Mod, c(mv_2018_sev, evS_2018_sev, evA_2018_sev,
                                            mv_bin_2018_sev, evS_bin_2018_sev, evA_bin_2018_sev,
                                            mv_growth, evS_growth, evA_growth,
                                            mv_bin_growth, evS_bin_growth, evA_bin_growth))
# no significant effects

# 2019
evS_2019_sev = "late_aug_severity = 0"
mv_2019_sev = "late_aug_severity + focm:late_aug_severity = 0"
evA_2019_sev = "late_aug_severity + foca:late_aug_severity = 0"
evS_bin_2019_sev = "hu_late_aug_severity = 0"
mv_bin_2019_sev = "hu_late_aug_severity + hu_focm:late_aug_severity = 0"
evA_bin_2019_sev = "hu_late_aug_severity + hu_foca:late_aug_severity = 0"

seedSevD2Coef <- hypothesis(seedSevD2Mod, c(mv_2019_sev, evS_2019_sev, evA_2019_sev,
                                            mv_bin_2019_sev, evS_bin_2019_sev, evA_bin_2019_sev,
                                            mv_growth, evS_growth, evA_growth,
                                            mv_bin_growth, evS_bin_growth, evA_bin_growth))
# no significant severity effects
# strong growth effects

# add columns
seedSevD1Coef2 <- seedSevD1Coef[[1]] %>%
  mutate(year = 2018,
         focal = rep(c("Mv", "Ev seedling", "Ev adult"), 4),
         response = rep(rep(c("continuous", "binary"), each = 3), 2),
         variable = rep(c("severity", "growth"), each = 6)) %>%
  select(year, focal, variable, response, Estimate:CI.Upper)

seedSevD2Coef2 <- seedSevD2Coef[[1]] %>%
  mutate(year = 2019,
         focal = rep(c("Mv", "Ev seedling", "Ev adult"), 4),
         response = rep(rep(c("continuous", "binary"), each = 3), 2),
         variable = rep(c("severity", "growth"), each = 6)) %>%
  select(year, focal, variable, response, Estimate:CI.Upper)

# combine
seedSevCoef <- seedSevD1Coef2 %>%
  full_join(seedSevD2Coef2) %>%
  arrange(year, focal, variable, response)

# save
write_csv(seedSevCoef, "output/seed_severity_growth_coefficients_2018_2019_density_exp.csv")


#### density models ####

# remove plot 1
seedD1Dat3 <- seedD1Dat %>%
  filter(plot > 1) %>%
  mutate(foc = fct_relevel(foc, "a"),
         bg = fct_relevel(bg, "a"))

seedD2Dat3 <- seedD2Dat %>%
  filter(plot > 1) %>%
  mutate(foc = fct_relevel(foc, "a"),
         bg = fct_relevel(bg, "a"))

# initial visualizations
ggplot(seedD1Dat3 %>% filter(seeds > 0), aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(1)) +
  facet_grid(focal ~ background, scales = "free")
# non non-zero seeds for Ev seedlings with low Ev seedling density

# fit models
# seedD1Mod <- brm(data = seedD1Dat3, family = hurdle_lognormal,
#                  bf(seeds ~ density * foc * bg * fungicide + (1|site),
#                     hu ~ density * foc * bg * fungicide + (1|site)),
#                  prior <- c(prior(normal(3.5, 1), class = "Intercept"),
#                             prior(normal(0, 10), class = "Intercept", dpar = "hu"),
#                             prior(normal(0, 10), class = "b")), # use default for sigma
#                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
# cannot fit this model - a lot of divergent transitions and transitions that exceed max treedepth
# I think because there's too much missing data when it's split

seedD1Mod <- brm(data = seedD1Dat3, family = lognormal,
                 seeds1 ~ density * foc * bg * fungicide + (1|site),
                 prior <- c(prior(normal(3.5, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedD1Mod)

seedD2Mod <- brm(data = seedD2Dat3, family = lognormal,
                 seeds1 ~ density * foc * bg * fungicide + (1|site),
                 prior <- c(prior(normal(3.5, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                 control = list(max_treedepth = 15)) 
mod_check_fun(seedD2Mod)

# save models
save(seedD1Mod, file = "output/focal_seed_density_model_2018_density_exp.rda")
save(seedD2Mod, file = "output/focal_seed_density_model_2019_density_exp.rda")


#### fungicide-only models ####

# average fungicide effect across background plot types
# keep plot 1

# initial visualizations
ggplot(seedD1Dat, aes(x = foc, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean")

ggplot(seedD2Dat, aes(x = foc, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean")

# fit models
seedFungD1Mod <- brm(data = seedD1Dat, family = hurdle_lognormal,
                    bf(seeds ~ foc * fungicide + (1|plotf),
                       hu ~ foc * fungicide + (1|plotf)),
                    prior <- c(prior(normal(1, 1), class = "Intercept"),
                               prior(normal(0, 10), class = "Intercept", dpar = "hu"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedFungD1Mod)

seedFungD2Mod <- brm(data = seedD2Dat, family = hurdle_lognormal,
                    bf(seeds ~ foc * fungicide + (1|plotf),
                       hu ~ foc * fungicide + (1|plotf)),
                    prior <- c(prior(normal(1, 1), class = "Intercept"),
                               prior(normal(0, 10), class = "Intercept", dpar = "hu"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedFungD2Mod)

# save models
save(seedFungD1Mod, file = "output/seed_fungicide_model_2018_density_exp.rda")
save(seedFungD2Mod, file = "output/seed_fungicide_model_2019_density_exp.rda")


#### fungicide effects ####

# fungicide - control when density = 0
mv_cont_fung_eff = "fungicide + focm:fungicide = 0"
evS_cont_fung_eff = "fungicide = 0"
evA_cont_fung_eff = "fungicide + foca:fungicide = 0"

mv_bin_fung_eff = "hu_fungicide + hu_focm:fungicide = 0"
evS_bin_fung_eff = "hu_fungicide = 0"
evA_bin_fung_eff = "hu_fungicide + hu_foca:fungicide = 0"

seedFungEff <- hypothesis(seedFungD1Mod,
                            c(mv_cont_fung_eff,
                              evS_cont_fung_eff,
                              evA_cont_fung_eff,
                              mv_bin_fung_eff,
                              evS_bin_fung_eff,
                              evA_bin_fung_eff))[[1]] %>%
  mutate(Year = 2018) %>%
  full_join(hypothesis(seedFungD2Mod,
                       c(mv_cont_fung_eff,
                         evS_cont_fung_eff,
                         evA_cont_fung_eff,
                         mv_bin_fung_eff,
                         evS_bin_fung_eff,
                         evA_bin_fung_eff))[[1]] %>%
              mutate(Year = 2019)) %>%
  mutate(Focal = rep(c("Mv", "Ev first-year", "Ev adult"), 4), 
         Distribution = rep(rep(c("continuous", "binary"), each = 3), 2)) %>%
  select(Year, Focal, Distribution, Estimate:CI.Upper) %>%
  arrange(Year, Focal, Distribution)

write_csv(seedFungEff, "output/focal_seed_fungicide_effect_2018_2019_density_exp.csv")


#### interaction coefficients (alphas) ####

# are alphas different than 0? (_alpha)
# does fungicide treatment affect alphas? (_trt_eff, fung_alpha - ctrl_alpha)

# EvA background
evA_evA_ctrl_alpha = "density = 0"
evA_evA_fung_alpha = "density + density:fungicide = 0"
evA_evA_trt_eff = "density:fungicide = 0"
evS_evA_ctrl_alpha = "density + density:focs = 0"
evS_evA_fung_alpha = "density + density:fungicide + density:focs + density:focs:fungicide = 0"
evS_evA_trt_eff = "density:fungicide + density:focs:fungicide = 0"
mv_evA_ctrl_alpha = "density + density:focm = 0"
mv_evA_fung_alpha = "density + density:fungicide + density:focm + density:focm:fungicide = 0"
mv_evA_trt_eff = "density:fungicide + density:focm:fungicide = 0"

# EvS background
evS_evS_ctrl_alpha = "density + density:focs + density:bgs + density:focs:bgs = 0"
evS_evS_fung_alpha = "density + density:focs + density:bgs + density:focs:bgs + density:fungicide + density:focs:fungicide + density:bgs:fungicide + density:focs:bgs:fungicide = 0"
evS_evS_trt_eff = "density:fungicide + density:focs:fungicide + density:bgs:fungicide + density:focs:bgs:fungicide = 0"
evA_evS_ctrl_alpha = "density +  density:bgs = 0"
evA_evS_fung_alpha = "density +  density:bgs + density:fungicide + density:bgs:fungicide = 0"
evA_evS_trt_eff = "density:fungicide + density:bgs:fungicide = 0"
mv_evS_ctrl_alpha = "density +  density:bgs + density:focm + density:focm:bgs = 0"
mv_evS_fung_alpha = "density +  density:bgs + density:focm + density:focm:bgs + density:fungicide +  density:bgs:fungicide + density:focm:fungicide + density:focm:bgs:fungicide = 0"
mv_evS_trt_eff = "density:fungicide +  density:bgs:fungicide + density:focm:fungicide + density:focm:bgs:fungicide = 0"

# Mv background
mv_mv_ctrl_alpha = "density + density:focm + density:bgm + density:focm:bgm = 0"
mv_mv_fung_alpha = "density + density:focm + density:bgm + density:focm:bgm + density:fungicide + density:focm:fungicide + density:bgm:fungicide + density:focm:bgm:fungicide = 0"
mv_mv_trt_eff = "density:fungicide + density:focm:fungicide + density:bgm:fungicide + density:focm:bgm:fungicide = 0"
evA_mv_ctrl_alpha = "density +  density:bgm = 0"
evA_mv_fung_alpha = "density +  density:bgm + density:fungicide + density:bgm:fungicide = 0"
evA_mv_trt_eff = "density:fungicide + density:bgm:fungicide = 0"
evS_mv_ctrl_alpha = "density + density:bgm + density:focs + density:focs:bgm = 0"
evS_mv_fung_alpha = "density + density:bgm + density:focs + density:focs:bgm + density:fungicide + density:bgm:fungicide + density:focs:fungicide + density:focs:bgm:fungicide = 0"
evS_mv_trt_eff = "density:fungicide + density:bgm:fungicide + density:focs:fungicide + density:focs:bgm:fungicide = 0"

seedD1alphas <- hypothesis(seedD1Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

seedD2alphas <- hypothesis(seedD2Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

seedD1TrtEff <- hypothesis(seedD1Mod,
                             c(mv_mv_trt_eff, 
                               evS_mv_trt_eff, 
                               evA_mv_trt_eff,
                               evS_evS_trt_eff,
                               mv_evS_trt_eff, 
                               evA_evS_trt_eff, 
                               evA_evA_trt_eff,
                               mv_evA_trt_eff,
                               evS_evA_trt_eff))

seedD2TrtEff <- hypothesis(seedD2Mod,
                           c(mv_mv_trt_eff, 
                             evS_mv_trt_eff, 
                             evA_mv_trt_eff,
                             evS_evS_trt_eff,
                             mv_evS_trt_eff, 
                             evA_evS_trt_eff, 
                             evA_evA_trt_eff,
                             mv_evA_trt_eff,
                             evS_evA_trt_eff))

# combine alphas
alphaDat <- seedD1alphas[[1]] %>%
  mutate(year = "2018",
         foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                        "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                        "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                        "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                        "s_a_ctrl", "s_a_fung")) %>%
  full_join(seedD2alphas[[1]] %>%
              mutate(year = "2019",
                     foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                                    "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                                    "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                                    "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                                    "s_a_ctrl", "s_a_fung"))) %>%
  select(-Hypothesis) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg_trt, "_")[[1]][1],
         bg = str_split(foc_bg_trt, "_")[[1]][2],
         trt = str_split(foc_bg_trt, "_")[[1]][3]) %>%
  ungroup() %>%
  mutate(treatment = fct_recode(trt, "control (water)" = "ctrl",
                                "fungicide" = "fung"),
         sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))

# edit to save
alphaDatSave <- alphaDat %>%
  left_join(seedD1Dat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)"),
         focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  select(year, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper, sig) %>%
  arrange(year, background, focal, treatment)

# combine treatment effects
trtEffDatSave <- seedD1TrtEff[[1]] %>%
  mutate(Year = 2018,
         Response = "seeds",
         Gradient = "density") %>%
  full_join(seedD2TrtEff[[1]] %>%
              mutate(Year = 2019,
                     Response = "seeds",
                     Gradient = "density")) %>%
  mutate(foc_bg = rep(c("m_m", "s_m", "a_m", "s_s", "m_s", "a_s", "a_a", "m_a", "s_a"), 2)) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg, "_")[[1]][1],
         bg = str_split(foc_bg, "_")[[1]][2]) %>%
  ungroup() %>%
  left_join(seedD1Dat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob, Star, foc_bg, foc, bg)) %>%
  relocate(Year, Response, Gradient, focal, background) %>%
  arrange(Year, background, focal)

# save
write_csv(alphaDatSave, "output/focal_seed_competition_coefficients_2018_2019_density_exp.csv")
write_csv(trtEffDatSave, "output/focal_seed_competition_treatment_effect_2018_2019_density_exp.csv")

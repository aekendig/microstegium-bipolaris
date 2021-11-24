##### info ####

# file: focal_growth_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 9/29/21
# goal: analyses of plant growth as a function of density and severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)

# import plot information
plotsD <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# import severity data
sevD1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R

# import plot biomass data
plotBioD2Dat <- read_csv("intermediate-data/plot_biomass_2019_density_exp.csv")

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
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1,
                             TRUE ~ background_density)) %>%
  select(plot, treatment, background, density)

plotDens2 <- plotsD2 %>%
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1,
                             TRUE ~ background_density)) %>%
  select(plot, treatment, background, density)

# missing data
filter(tillerD1Dat, is.na(tillers_jul)) # 0
filter(tillerD1Dat, is.na(tillers_jun)) # 0
filter(tillerD1Dat, tillers_jul == 0 | tillers_jun == 0) # these plants are dead
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling
filter(plotBioD2Dat, is.na(biomass.g_m2)) # all no bg plots

# combine data
growthD1Dat <- tillerD1Dat %>%
  # left_join(plotDens) %>% # only need for initial viz in next section
  left_join(plotDens2) %>%
  mutate(plant_growth = log(tillers_jul/tillers_jun),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1))) %>%
  filter(tillers_jun > 0 & tillers_jul > 0) %>%
  left_join(sevD1Dat2)

growthD2Dat <- mvBioD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  # left_join(plotDens) %>% # only need for initial viz in next section
  left_join(plotDens2) %>% 
  mutate(plant_growth = log(biomass_weight.g),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1))) %>%
  filter(!is.na(biomass_weight.g)) %>%
  left_join(sevD2Dat2) %>%
  left_join(plotBioD2Dat %>%
              select(site, treatment, plot, biomass.g_m2) %>%
              rename(plot_biomass = biomass.g_m2)) %>%
  mutate(plot_biomass = if_else(as.character(focal) == as.character(background), 
                                plot_biomass - biomass_weight.g, plot_biomass))


#### competition models ####

# initial visualization
ggplot(growthD1Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

ggplot(growthD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

# remove plot 1
growthD1Dat2 <- growthD1Dat %>%
  filter(plot > 1)

growthD2Dat2 <- growthD2Dat %>%
  filter(plot > 1)

# biomass visualization
ggplot(growthD2Dat2, aes(plot_biomass, plant_growth, color = treatment)) +
  geom_point() +
  facet_grid(focal ~ background, scales = "free")

# fit models
growthD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (density + density:bg) + (1|plotf),
                   prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.9999, max_treedepth = 15)) 
# 36 divergent transitions, 864 transitions exceeded max_treedepth
mod_check_fun(growthD1Mod)

growthD2Mod <- brm(data = growthD2Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (density + density:bg) + (1|plotf),
                   prior <- c(prior(normal(3, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99999, max_treedepth = 15)) 
mod_check_fun(growthD2Mod)

# biomass model
growthD2Mod2 <- brm(data = growthD2Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (plot_biomass + plot_biomass:bg) + (1|plotf),
                   prior <- c(prior(normal(3, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99999, max_treedepth = 15)) 
mod_check_fun(growthD2Mod2)

# subset data to remove high EvA biomass
growthD2Dat2b <- growthD2Dat2 %>%
  filter(!(bg == "a" & plot_biomass > 150))

growthD2Mod2b <- update(growthD2Mod2, newdata = growthD2Dat2b)
summary(growthD2Mod2b)

# save models and data
save(growthD1Mod, file = "output/focal_growth_density_model_2018_density_exp.rda")
save(growthD2Mod, file = "output/focal_growth_density_model_2019_density_exp.rda")
save(growthD2Mod2, file = "output/focal_growth_biomass_model_2019_density_exp.rda")
save(growthD2Mod2b, file = "output/focal_growth_biomass_model_no_high_EvA_2019_density_exp.rda")
write_csv(growthD1Dat2, "intermediate-data/focal_growth_density_data_2018_density_exp.csv")
write_csv(growthD2Dat2, "intermediate-data/focal_growth_density_data_2019_density_exp.csv")


#### interaction coefficients (alphas) ####

# are alphas different than 0? (_alpha)
# does fungicide treatment affect alphas? (_trt_eff, fung_alpha - ctrl_alpha)

# Mv background
mv_mv_ctrl_alpha = "density = 0"
mv_mv_fung_alpha = "density + fungicide:density = 0"
mv_mv_trt_eff = "fungicide:density = 0"
evS_mv_ctrl_alpha = "density + focs:density = 0"
evS_mv_fung_alpha = "density + fungicide:density + focs:density + focs:fungicide:density = 0"
evS_mv_trt_eff = "fungicide:density + focs:fungicide:density = 0"
evA_mv_ctrl_alpha = "density + foca:density = 0"
evA_mv_fung_alpha = "density + fungicide:density + foca:density + foca:fungicide:density = 0"
evA_mv_trt_eff = "fungicide:density + foca:fungicide:density = 0"

# EvS background
evS_evS_ctrl_alpha = "density + focs:density + density:bgs + focs:density:bgs = 0"
evS_evS_fung_alpha = "density + focs:density + density:bgs + focs:density:bgs + fungicide:density + focs:fungicide:density + fungicide:density:bgs + focs:fungicide:density:bgs = 0"
evS_evS_trt_eff = "fungicide:density + focs:fungicide:density + fungicide:density:bgs + focs:fungicide:density:bgs = 0"
mv_evS_ctrl_alpha = "density +  density:bgs = 0"
mv_evS_fung_alpha = "density +  density:bgs + fungicide:density + fungicide:density:bgs = 0"
mv_evS_trt_eff = "fungicide:density + fungicide:density:bgs = 0"
evA_evS_ctrl_alpha = "density +  density:bgs + foca:density + foca:density:bgs = 0"
evA_evS_fung_alpha = "density +  density:bgs + foca:density + foca:density:bgs + fungicide:density +  fungicide:density:bgs + foca:fungicide:density + foca:fungicide:density:bgs = 0"
evA_evS_trt_eff = "fungicide:density +  fungicide:density:bgs + foca:fungicide:density + foca:fungicide:density:bgs = 0"

# EvA background
evA_evA_ctrl_alpha = "density + foca:density + density:bga + foca:density:bga = 0"
evA_evA_fung_alpha = "density + foca:density + density:bga + foca:density:bga + fungicide:density + foca:fungicide:density + fungicide:density:bga + foca:fungicide:density:bga = 0"
evA_evA_trt_eff = "fungicide:density + foca:fungicide:density + fungicide:density:bga + foca:fungicide:density:bga = 0"
mv_evA_ctrl_alpha = "density +  density:bga = 0"
mv_evA_fung_alpha = "density +  density:bga + fungicide:density + fungicide:density:bga = 0"
mv_evA_trt_eff = "fungicide:density + fungicide:density:bga = 0"
evS_evA_ctrl_alpha = "density + density:bga + focs:density + focs:density:bga = 0"
evS_evA_fung_alpha = "density + density:bga + focs:density + focs:density:bga + fungicide:density + fungicide:density:bga + focs:fungicide:density + focs:fungicide:density:bga = 0"
evS_evA_trt_eff = "fungicide:density + fungicide:density:bga + focs:fungicide:density + focs:fungicide:density:bga = 0"

densityD1alphas <- hypothesis(growthD1Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

densityD1TrtEff <- hypothesis(growthD1Mod,
                          c(mv_mv_trt_eff, 
                            evS_mv_trt_eff, 
                            evA_mv_trt_eff,
                            evS_evS_trt_eff,
                            mv_evS_trt_eff, 
                            evA_evS_trt_eff, 
                            evA_evA_trt_eff,
                            mv_evA_trt_eff,
                            evS_evA_trt_eff))

densityD2alphas <- hypothesis(growthD2Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

densityD2TrtEff <- hypothesis(growthD2Mod,
                              c(mv_mv_trt_eff, 
                                evS_mv_trt_eff, 
                                evA_mv_trt_eff,
                                evS_evS_trt_eff,
                                mv_evS_trt_eff, 
                                evA_evS_trt_eff, 
                                evA_evA_trt_eff,
                                mv_evA_trt_eff,
                                evS_evA_trt_eff))


# Mv background
mv_mv_ctrl_alpha = "plot_biomass = 0"
mv_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass = 0"
mv_mv_trt_eff = "fungicide:plot_biomass = 0"
evS_mv_ctrl_alpha = "plot_biomass + focs:plot_biomass = 0"
evS_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + focs:plot_biomass + focs:fungicide:plot_biomass = 0"
evS_mv_trt_eff = "fungicide:plot_biomass + focs:fungicide:plot_biomass = 0"
evA_mv_ctrl_alpha = "plot_biomass + foca:plot_biomass = 0"
evA_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + foca:plot_biomass + foca:fungicide:plot_biomass = 0"
evA_mv_trt_eff = "fungicide:plot_biomass + foca:fungicide:plot_biomass = 0"

# EvS background
evS_evS_ctrl_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs = 0"
evS_evS_fung_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs + fungicide:plot_biomass + focs:fungicide:plot_biomass + fungicide:plot_biomass:bgs + focs:fungicide:plot_biomass:bgs = 0"
evS_evS_trt_eff = "fungicide:plot_biomass + focs:fungicide:plot_biomass + fungicide:plot_biomass:bgs + focs:fungicide:plot_biomass:bgs = 0"
mv_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs = 0"
mv_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + fungicide:plot_biomass + fungicide:plot_biomass:bgs = 0"
mv_evS_trt_eff = "fungicide:plot_biomass + fungicide:plot_biomass:bgs = 0"
evA_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs = 0"
evA_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs + fungicide:plot_biomass +  fungicide:plot_biomass:bgs + foca:fungicide:plot_biomass + foca:fungicide:plot_biomass:bgs = 0"
evA_evS_trt_eff = "fungicide:plot_biomass +  fungicide:plot_biomass:bgs + foca:fungicide:plot_biomass + foca:fungicide:plot_biomass:bgs = 0"

# EvA background
evA_evA_ctrl_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga = 0"
evA_evA_fung_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga + fungicide:plot_biomass + foca:fungicide:plot_biomass + fungicide:plot_biomass:bga + foca:fungicide:plot_biomass:bga = 0"
evA_evA_trt_eff = "fungicide:plot_biomass + foca:fungicide:plot_biomass + fungicide:plot_biomass:bga + foca:fungicide:plot_biomass:bga = 0"
mv_evA_ctrl_alpha = "plot_biomass +  plot_biomass:bga = 0"
mv_evA_fung_alpha = "plot_biomass +  plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga = 0"
mv_evA_trt_eff = "fungicide:plot_biomass + fungicide:plot_biomass:bga = 0"
evS_evA_ctrl_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga = 0"
evS_evA_fung_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga + focs:fungicide:plot_biomass + focs:fungicide:plot_biomass:bga = 0"
evS_evA_trt_eff = "fungicide:plot_biomass + fungicide:plot_biomass:bga + focs:fungicide:plot_biomass + focs:fungicide:plot_biomass:bga = 0"

biomassD2alphas <- hypothesis(growthD2Mod2, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

biomassD2TrtEff <- hypothesis(growthD2Mod2,
                             c(mv_mv_trt_eff, 
                               evS_mv_trt_eff, 
                               evA_mv_trt_eff,
                               evS_evS_trt_eff,
                               mv_evS_trt_eff, 
                               evA_evS_trt_eff, 
                               evA_evA_trt_eff,
                               mv_evA_trt_eff,
                               evS_evA_trt_eff))

biomassD2balphas <- hypothesis(growthD2Mod2b, 
                              c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                                evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                                evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                                evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                                mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                                evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                                evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                                mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                                evS_evA_ctrl_alpha, evS_evA_fung_alpha))

biomassD2bTrtEff <- hypothesis(growthD2Mod2b,
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
alphaDat <- densityD1alphas[[1]] %>%
  mutate(year = 2018,
         response = "tillers",
         gradient = "density",
         data = "full") %>%
  full_join(densityD2alphas[[1]] %>%
              mutate(year = 2019,
                     response = "biomass",
                     gradient = "density",
                     data = "full")) %>%
  full_join(biomassD2alphas[[1]] %>%
              mutate(year = 2019,
                     response = "biomass",
                     gradient = "biomass",
                     data = "full")) %>%
  full_join(biomassD2balphas[[1]] %>%
              mutate(year = 2019,
                     response = "biomass",
                     gradient = "biomass",
                     data = "no high Ev A")) %>%
  mutate(foc_bg_trt = rep(c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                            "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                            "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                            "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                            "s_a_ctrl", "s_a_fung"), 4)) %>%
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
  left_join(growthD2Dat2 %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)"),
         focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  select(year, response, gradient, data, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper, sig) %>%
  arrange(year, response, gradient, data, background, focal, treatment)

# combine treatment effects
trtEffDatSave <- densityD1TrtEff[[1]] %>%
  mutate(Year = 2018,
         Data = "full",
         Response = "tillers",
         Gradient = "density") %>%
  full_join(densityD2TrtEff[[1]] %>%
              mutate(Year = 2019,
                     Data = "full",
                     Response = "biomass",
                     Gradient = "density")) %>%
  full_join(biomassD2TrtEff[[1]] %>%
              mutate(Year = 2019,
                     Data = "full",
                     Response = "biomass",
                     Gradient = "biomass")) %>%
  full_join(biomassD2bTrtEff[[1]] %>%
              mutate(Year = 2019,
                     Data = "no high Ev A",
                     Response = "biomass",
                     Gradient = "biomass")) %>%
  mutate(foc_bg = rep(c("m_m", "s_m", "a_m", "s_s", "m_s", "a_s", "a_a", "m_a", "s_a"), 4)) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg, "_")[[1]][1],
         bg = str_split(foc_bg, "_")[[1]][2]) %>%
  ungroup() %>%
  left_join(growthD2Dat2 %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling"),
         sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0")) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob, Star, foc_bg, foc, bg)) %>%
  relocate(Year, Response, Gradient, Data, focal, background) %>%
  arrange(Year, Response, Gradient, Data, background, background, focal)

# save
write_csv(alphaDatSave, "output/focal_growth_competition_coefficients_2018_2019_density_exp.csv")
write_csv(trtEffDatSave, "output/focal_growth_competition_treatment_effect_2018_2019_density_exp.csv")


#### growth rate ####

mv_trt_eff = "fungicide = 0"
evS_trt_eff = "fungicide + focs:fungicide = 0"
evA_trt_eff = "fungicide + foca:fungicide = 0"

growthDensityD1TrtEff <- hypothesis(growthD1Mod, c(mv_trt_eff, evS_trt_eff, evA_trt_eff))
growthDensityD2TrtEff <- hypothesis(growthD2Mod, c(mv_trt_eff, evS_trt_eff, evA_trt_eff))
growthBiomassD2TrtEff <- hypothesis(growthD2Mod2, c(mv_trt_eff, evS_trt_eff, evA_trt_eff))

growthTrtEff <- growthDensityD1TrtEff[[1]] %>%
  mutate(year = 2018,
         gradient = "density") %>%
  full_join(growthDensityD2TrtEff[[1]] %>%
              mutate(year = 2019,
                     gradient = "density")) %>%
  full_join(growthBiomassD2TrtEff[[1]] %>%
              mutate(year = 2019,
                     gradient = "biomass")) %>%
  mutate(focal = rep(c("Mv", "Ev first-year", "Ev adult"), 3)) %>%
  select(-Hypothesis) %>%
  relocate(year, gradient, focal) %>%
  arrange(year, gradient, focal)

write_csv(growthTrtEff, "output/focal_growth_rate_treatment_effect_2018_2019_density_exp.csv")


#### figure ####

# density function
dens_fun <- function(f, b, trt, yr){
  
  if(yr == "2018"){
    dat <- growthD1Dat2 %>% filter(foc == f & bg == b & treatment == trt)
    density <- seq(min(dat$density), max(dat$density), length.out = 100)
  }else{
    dat <- growthD2Dat2 %>% filter(foc == f & bg == b & treatment == trt)
    density <- seq(min(dat$plot_biomass), max(dat$plot_biomass), length.out = 100)
  }
  
  return(density)
}

# predicted data
predD1Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  mutate(plant_growth = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predD2Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2019",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  rename(plot_biomass = density) %>%
  mutate(plant_growth = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, Mv = "m", "Ev first-year" = "s", "Ev adult" = "a") %>%
           fct_relevel("Ev adult", "Ev first-year"),
         background = fct_recode(bg, Mv = "m", "Ev first-year" = "s", "Ev adult" = "a") %>%
           fct_relevel("Ev adult", "Ev first-year"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

# combine with preddat
predDat2 <- predDat %>%
  left_join(alphaDat) %>%
  mutate(intra = case_when(foc == bg ~ "yes",
                           foc %in% c("a", "s") & bg %in% c("a", "s") ~ "yes",
                           TRUE ~ "no"))

# raw data
figDat <- growthD1Dat2 %>%
  select(site, plot, treatment, sp, age, ID, focal, background, density, plant_growth) %>%
  mutate(year = "2018") %>%
  full_join(growthD2Dat2 %>%
              select(site, plot, treatment, sp, age, ID, focal, background, plot_biomass, plant_growth) %>%
              mutate(year = "2019")) %>%
  left_join(plotsD %>%
              select(plot, density_level) %>%
              unique()) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling"),
         density_level = fct_relevel(density_level, "low", "medium"),
         intra = case_when(focal ==  background ~ "yes",
                           focal %in% c("Ev adult", "Ev first-year") &  background %in% c("Ev adult", "Ev first-year") ~ "yes",
                           TRUE ~ "no"))

# alphas for figure
# sig beta values
alphaDat2 <- alphaDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, foc, focal, bg, background) %>%
              summarise(density = max(density, na.rm = T),
                        plot_biomass = max(plot_biomass, na.rm = T))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(lower = min(CI.Lower)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(growth = min(plant_growth)) %>%
                          ungroup()) %>%
              rowwise() %>%
              mutate(plant_growth = min(c(growth, lower))) %>%
              ungroup() %>%
              select(-c(growth, lower))) %>%
  mutate(comp = as.character(sprintf("%.3f",round(Estimate, 3))))


# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")
textSize = 3
box_shade = "gray92"

# yearText <- tibble(year = c("2018", "2019"),
#                    density = c(3.1, 3.1),
#                    plant_growth = c(0.87, 3.09),
#                    background = "Ev adult",
#                    focal = "Ev adult",
#                    treatment = "fungicide")



# legend
legFig <- ggplot(predDat2, aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = figDat, aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.5) +
  facet_grid(rows = vars(focal),
             cols = vars(background)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density treatment") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical") +
  guides(shape = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

leg <- get_legend(legFig)


# 2018
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = plant_growth)) +
  geom_rect(aes(fill = intra), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c(col_pal, "white", box_shade)) +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme

# print
tiff("output/focal_growth_pairwise_figure_2018_density_exp.tiff", width = 9, height = 10, units = "cm", res = 300, compression = "lzw")
plot_grid(pairD1Fig, leg,
          nrow = 2, 
          rel_heights = c(1, 0.15))
dev.off()

# 2019
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = log(plot_biomass), y = plant_growth)) +
  geom_rect(aes(fill = intra), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat2, year == "2019"), 
            aes(label = paste("alpha", "==", comp, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 0.2, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c(col_pal, "white", box_shade)) +
  scale_x_continuous(breaks = c(2, 3, 4, 5, 6)) +
  xlab(expression(paste("Competitor biomass (ln[g/", m^2, "])", sep = ""))) +
  ylab("Focal biomass (ln[g])") +
  fig_theme

# compD2Fig <- ggplot(growthD2hyps, aes(x = Competitor, y = estimate, color = treatment)) +
#   geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(0.4)) +
#   geom_point(size = 2, position = position_dodge(0.4)) +
#   scale_color_manual(values = col_pal, name = "Treatment") +
#   ylab(expression(paste("Inter - intra competition ( ", 10^-3, ")", sep = ""))) +
#   fig_theme +
#   theme(axis.title.y = element_text(size = 10, hjust = 0))

# legend
# leg <- get_legend(pairD2Fig)
# 
# # right column
# rightFig <- plot_grid(compD2Fig, leg,
#                       ncol = 1,
#                       labels = c("B", NA),
#                       rel_heights = c(1, 0.4))

# combine
tiff("output/focal_growth_pairwise_figure_2019_density_exp.tiff", width = 9, height = 10, units = "cm", res = 300, compression = "lzw")
# plot_grid(pairD2Fig + theme(legend.position = "none"), rightFig,
#           nrow = 1, ncol = 2,
#           labels = c("A", NA),
#           rel_widths = c(1, 0.4))
plot_grid(pairD2Fig, leg,
          nrow = 2, 
          rel_heights = c(1, 0.15))
dev.off()


#### supplementary figure ####

yearText <- tibble(year = c("2018", "2019"),
                   density = c(-0.1, -0.1),
                   plant_growth = c(0.45, 2.3),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

# 2018
rawD1Fig <- growthD1Dat %>%
  mutate(focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  ggplot(aes(x = density, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_color_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# 2019
rawD2Fig <- growthD2Dat %>%
  mutate(focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  ggplot(aes(x = density, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg2 <- get_legend(rawD1Fig)

# combine plots
combFig2 <- plot_grid(rawD1Fig + theme(legend.position = "none"), rawD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9),
                     label_x = c(0, -0.01))

# combine
tiff("output/focal_raw_pairwise_figure_2018_2019_density_exp.tiff", width = 6.5, height = 3.5, units = "in", res = 300, compression = "lzw")
plot_grid(combFig2, leg2,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()


#### severity initial visualizations ####

# remove repeat plot 1's (will work whether or not they're present)
growthD1Dat3 <- growthD1Dat %>%
  filter(!(plot == 1 & background %in% c("Ev seedling", "Ev adult"))) %>%
  mutate(background = case_when(plot == 1 ~ "none",
                                TRUE ~ as.character(background)))

growthD2Dat3 <- growthD2Dat %>%
  filter(!(plot == 1 & background %in% c("Ev seedling", "Ev adult"))) %>%
  mutate(background = case_when(plot == 1 ~ "none",
                                TRUE ~ as.character(background)))

# initial visualizations
growthD1Dat3 %>%
  filter(sp == "Mv") %>%
  select(jul_severity, late_aug_severity, sep_severity, plant_growth) %>%
  ggpairs()

growthD1Dat3 %>%
  filter(focal == "Ev seedling") %>%
  select(jul_severity, late_aug_severity, sep_severity, plant_growth) %>%
  ggpairs()

growthD1Dat3 %>%
  filter(focal == "Ev adult") %>%
  select(jul_severity, late_aug_severity, sep_severity, plant_growth) %>%
  ggpairs()
# really only makes sense to use July severity because growth was measured
# june to july, and those severity-growth relationships are week

growthD2Dat3 %>%
  filter(sp == "Mv") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD2Dat3 %>%
  filter(focal == "Ev seedling") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD2Dat3 %>%
  filter(focal == "Ev adult") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()
# mostly positive relationships between growth and severity
# interesting to get this result after no sig transmission this year
# maybe driven by 1 plots?

# visualizations excluding 1 plots
growthD1Dat2 %>%
  filter(sp == "Mv") %>%
  select(jul_severity, sep_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD1Dat3 %>%
  filter(focal == "Ev seedling") %>%
  select(jul_severity, sep_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD1Dat3 %>%
  filter(focal == "Ev adult") %>%
  select(jul_severity, sep_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD2Dat2 %>%
  filter(sp == "Mv") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD2Dat3 %>%
  filter(focal == "Ev seedling") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()

growthD2Dat3 %>%
  filter(focal == "Ev adult") %>%
  select(jul_severity, early_aug_severity, late_aug_severity, plant_growth) %>%
  ggpairs()
# still positive
# may be because plots that are good for growth are good for disease

# severity months
sum(!is.na(growthD1Dat2$jul_severity)) # most complete, broadest for Ev
sum(!is.na(growthD1Dat2$late_aug_severity))
sum(!is.na(growthD1Dat2$sep_severity)) # broadest for Mv

sum(!is.na(growthD2Dat2$jul_severity)) # most complete
sum(!is.na(growthD2Dat2$early_aug_severity))
sum(!is.na(growthD2Dat2$late_aug_severity)) # broadest distribution

# only use July for 2018 because growth measurements were taken in July
# doesn't make sense for severity later in the season to affect early growth


#### fit severity models ####

# July
growthJulSevD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                      plant_growth ~ foc * (jul_severity + density + density:bg) + (1|plotf),
                      prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                                 prior(normal(0, 1), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(growthJulSevD1Mod)

growthJulSevD2Mod <- brm(data = growthD2Dat2, family = gaussian,
                      plant_growth ~ foc * (jul_severity + density + density:bg) + (1|plotf),
                      prior <- c(prior(normal(3, 1), class = "Intercept"),
                                 prior(normal(0, 1), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                      control = list(max_treedepth = 15)) 
mod_check_fun(growthJulSevD2Mod)

# early Aug
growthEAugSevD2Mod <- brm(data = growthD2Dat2, family = gaussian,
                         plant_growth ~ foc * (early_aug_severity + density + density:bg) + (1|plotf),
                         prior <- c(prior(normal(3, 1), class = "Intercept"),
                                    prior(normal(0, 1), class = "b")), # use default for sigma
                         iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                         control = list(max_treedepth = 15)) 
mod_check_fun(growthEAugSevD2Mod)

# late Aug
growthLAugSevD2Mod <- brm(data = growthD2Dat2, family = gaussian,
                          plant_growth ~ foc * (late_aug_severity + density + density:bg) + (1|plotf),
                          prior <- c(prior(normal(3, 1), class = "Intercept"),
                                     prior(normal(0, 1), class = "b")), # use default for sigma
                          iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                          control = list(max_treedepth = 15)) 
mod_check_fun(growthLAugSevD2Mod)

# save models
save(growthJulSevD1Mod, file = "output/focal_growth_jul_severity_model_2018_density_exp.rda")
save(growthJulSevD2Mod, file = "output/focal_growth_jul_severity_model_2019_density_exp.rda")
save(growthEAugSevD2Mod, file = "output/focal_growth_early_aug_severity_model_2019_density_exp.rda")
save(growthLAugSevD2Mod, file = "output/focal_growth_late_aug_severity_model_2019_density_exp.rda")


#### severity coefficients ####

# July
mv_jul_sev = "jul_severity = 0"
evS_jul_sev = "jul_severity + focs:jul_severity = 0"
evA_jul_sev = "jul_severity + foca:jul_severity = 0"

# early Aug
mv_early_aug_sev = "early_aug_severity = 0"
evS_early_aug_sev = "early_aug_severity + focs:early_aug_severity = 0"
evA_early_aug_sev = "early_aug_severity + foca:early_aug_severity = 0"

# late Aug
mv_late_aug_sev = "late_aug_severity = 0"
evS_late_aug_sev = "late_aug_severity + focs:late_aug_severity = 0"
evA_late_aug_sev = "late_aug_severity + foca:late_aug_severity = 0"

growthJulSevD1coefs <- hypothesis(growthJulSevD1Mod, c(mv_jul_sev, evS_jul_sev, evA_jul_sev))
growthJulSevD2coefs <- hypothesis(growthJulSevD2Mod, c(mv_jul_sev, evS_jul_sev, evA_jul_sev))
# reduced seedling growth
growthEAugSevD2coefs <- hypothesis(growthEAugSevD2Mod, c(mv_early_aug_sev, evS_early_aug_sev, evA_early_aug_sev))
# increased Mv and adult growth
growthLAugSevD2coefs <- hypothesis(growthLAugSevD2Mod, c(mv_late_aug_sev, evS_late_aug_sev, evA_late_aug_sev))
# increased Mv and seedling growth

# save
sevCoefs <- growthJulSevD1coefs[[1]] %>%
  mutate(year = 2018,
         month = "July") %>%
  full_join(growthJulSevD2coefs[[1]] %>%
              mutate(year = 2019,
                     month = "July")) %>%
  full_join(growthEAugSevD2coefs[[1]] %>%
              mutate(year = 2019,
                     month = "early August")) %>%
  full_join(growthLAugSevD2coefs[[1]] %>%
              mutate(year = 2019,
                     month = "late August")) %>%
  mutate(focal = rep(c("Mv", "Ev first-year", "Ev adult"), 4)) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob)) %>%
  relocate(year, month, focal)

# save
write_csv(sevCoefs, "output/focal_growth_severity_effect_2018_2019_density_exp.csv")


#### fungicide-only models ####

# average fungicide effect across background plot types
# keep plot 1

# initial visualizations
ggplot(growthD1Dat, aes(x = foc, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean")

ggplot(growthD2Dat, aes(x = foc, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean")

# fit models
growthFungD1Mod <- brm(data = growthD1Dat, family = gaussian,
                       plant_growth ~ foc * fungicide + (1|plotf),
                       prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                                  prior(normal(0, 1), class = "b")), # use default for sigma
                       iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                       control = list(adapt_delta = 0.999)) 
mod_check_fun(growthFungD1Mod)

growthFungD2Mod <- brm(data = growthD2Dat, family = gaussian,
                       plant_growth ~ foc * fungicide + (1|plotf),
                       prior <- c(prior(normal(2, 1), class = "Intercept"),
                                  prior(normal(0, 1), class = "b")), # use default for sigma
                       iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(growthFungD2Mod)

# save models
save(growthFungD1Mod, file = "output/focal_growth_fungicide_model_2018_density_exp.rda")
save(growthFungD2Mod, file = "output/focal_growth_fungicide_model_2019_density_exp.rda")


#### fungicide effects ####

# fungicide - control
mv_fung_eff = "fungicide = 0"
evS_fung_eff = "fungicide + focs:fungicide = 0"
evA_fung_eff = "fungicide + foca:fungicide = 0"

growthFungEff <- hypothesis(growthFungD1Mod,
                            c(mv_fung_eff,
                              evS_fung_eff,
                              evA_fung_eff))[[1]] %>%
  mutate(Year = 2018, Response = "tillers") %>%
  full_join(hypothesis(growthFungD2Mod,
                       c(mv_fung_eff,
                         evS_fung_eff,
                         evA_fung_eff))[[1]] %>%
              mutate(Year = 2019, Response = "biomass")) %>%
  mutate(Focal = rep(c("Mv", "Ev first-year", "Ev adult"), 2)) %>%
  select(Year, Response, Focal, Estimate:CI.Upper) %>%
  arrange(Year, Focal)

write_csv(growthFungEff, "output/focal_growth_fungicide_effect_2018_2019_density_exp.csv")
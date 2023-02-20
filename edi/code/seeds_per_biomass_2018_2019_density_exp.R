##### outputs ####

# mv_seeds_per_biomass_untransformed_model_2019_density_exp.rda
# evS_seeds_per_biomass_untransformed_model_2019_density_exp.rda
# evA_seeds_per_biomass_untransformed_model_2019_density_exp.rda


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import biomass data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# import seed data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}

mod_fit_fun <- function(dat, mod, treatCol, ranEff = site, xCol, minX, maxX, yCol, f2t = F){
  
  outDat <- dat %>%
    select({{treatCol}}) %>%
    unique() %>%
    mutate("{{ranEff}}" := "D5") %>%
    expand_grid(tibble("{{xCol}}" := seq(minX, maxX, length.out = 100))) %>%
    mutate(pred = fitted(mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(mod, newdata = ., allow_new_levels = T)[, "Q97.5"])
  
  if(f2t == T){
    outDat <- outDat %>%
      mutate(fungicide = dplyr::recode(fungicide, "0" = "water", "1" = "fungicide"))
    
    dat <- dat %>%
      mutate(fungicide = dplyr::recode(fungicide, "0" = "water", "1" = "fungicide"))
  }
  
  outPlot <- ggplot(dat, aes(x = {{xCol}}, y = {{yCol}}, color = {{treatCol}}, fill = {{treatCol}})) +
    geom_ribbon(data = outDat, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(data = outDat, aes(y = pred)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
    stat_summary(geom = "point", fun = "mean", size = 2) +
    theme_bw()
  
  return(list(outDat, outPlot))
  
}


#### edit data ####

# Ev dat
evD2Dat <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  rename(biomass_weight.g = weight)

# all data
seedsBioD2Dat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g)) %>%
  mutate(ID = as.character(plant),
         age = "seedling")  %>%
  select(-plant) %>%
  full_join(evD2Dat) %>%
  filter(!is.na(biomass_weight.g) & !is.na(seeds)) %>%
  mutate(log_seeds = log(seeds + 1),
         log_bio = log(biomass_weight.g),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""))


#### initial visualizations ####

ggplot(seedsBioD2Dat, aes(biomass_weight.g, seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  facet_wrap(~ plant_group, scales = "free")

# divide data
mvSeedsBioD2Dat2 <- seedsBioD2Dat %>%
  filter(plant_group == "Mv_seedling")

evSSeedsBioD2Dat2 <- seedsBioD2Dat %>%
  filter(plant_group == "Ev_seedling")

evASeedsBioD2Dat2 <- seedsBioD2Dat %>%
  filter(plant_group == "Ev_adult")


### fit models ###
mvSeedsBioD2Mod2 <- brm(data = mvSeedsBioD2Dat2, family = gaussian,
                           seeds ~ 0 + biomass_weight.g + fungicide:biomass_weight.g + (1|plotf),
                           prior <- c(prior(normal(0, 100), class = "b")), # use default for sigma
                           iter = 6000, warmup = 1000, chains = 3)
mod_check_fun(mvSeedsBioD2Mod2)

evSSeedsBioD2Mod2 <- brm(data = evSSeedsBioD2Dat2, family = gaussian,
                         seeds ~ 0 + biomass_weight.g + fungicide:biomass_weight.g + (1|plotf),
                         prior <- c(prior(normal(0, 100), class = "b")), # use default for sigma
                         iter = 6000, warmup = 1000, chains = 3)
mod_check_fun(evSSeedsBioD2Mod2)

evASeedsBioD2Mod2 <- brm(data = evASeedsBioD2Dat2, family = gaussian,
                         seeds ~ 0 + biomass_weight.g + fungicide:biomass_weight.g + (1|site),
                         prior <- c(prior(normal(0, 100), class = "b")), # use default for sigma
                         iter = 6000, warmup = 1000, chains = 3,
                         control = list(adapt_delta = 0.999))
mod_check_fun(evASeedsBioD2Mod2)

# simulated data
mvSeedsBioD2Sim2 <- mod_fit_fun(dat = mvSeedsBioD2Dat2, mod = mvSeedsBioD2Mod2, treatCol = fungicide,
                                xCol = biomass_weight.g, 
                                minX = min(mvSeedsBioD2Dat2$biomass_weight.g), 
                                maxX = max(mvSeedsBioD2Dat2$biomass_weight.g), 
                                yCol = seeds, f2t = T)
mvSeedsBioD2Sim2[[2]]

evSSeedsBioD2Sim2 <- mod_fit_fun(dat = evSSeedsBioD2Dat2, mod = evSSeedsBioD2Mod2, treatCol = fungicide,
                                 xCol = biomass_weight.g, 
                                 minX = min(evSSeedsBioD2Dat2$biomass_weight.g), 
                                 maxX = max(evSSeedsBioD2Dat2$biomass_weight.g), 
                                 yCol = seeds, f2t = T)
evSSeedsBioD2Sim2[[2]]

evASeedsBioD2Sim2 <- mod_fit_fun(dat = evASeedsBioD2Dat2, mod = evASeedsBioD2Mod2, treatCol = fungicide,
                                xCol = biomass_weight.g, 
                                minX = min(evASeedsBioD2Dat2$biomass_weight.g), 
                                maxX = max(evASeedsBioD2Dat2$biomass_weight.g), 
                                yCol = seeds, f2t = T)
evASeedsBioD2Sim2[[2]]

# save
save(mvSeedsBioD2Mod2, file = "output/mv_seeds_per_biomass_untransformed_model_2019_density_exp.rda")
save(evSSeedsBioD2Mod2, file = "output/evS_seeds_per_biomass_untransformed_model_2019_density_exp.rda")
save(evASeedsBioD2Mod2, file = "output/evA_seeds_per_biomass_untransformed_model_2019_density_exp.rda")
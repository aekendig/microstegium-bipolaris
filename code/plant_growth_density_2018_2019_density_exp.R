##### info ####

# file: plant_growth_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/4/21
# goal: analyses of plant growth as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


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

# combine and separate data
mvD1Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(growthD1Dat) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev())

mvD2Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(mvBioD2Dat %>%
              rename(ID = plant) %>%
              mutate(ID = as.character(ID)) %>%
              full_join(evBioD2Dat %>%
                          rename(biomass_weight.g = weight))) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev(),
         pg_trt = paste(substr(sp, 1, 1), substr(age, 1, 1), substr(treatment, 1, 1), sep = "_") %>%
           fct_rev())


#### visualizations ####

# 2018 Mv density height and tiller
ggplot(mvD1Dat, aes(Mv_seedling_density, height_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(mvD1Dat, aes(Mv_seedling_density, tiller_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")
# not a strong effect on either

# 2019 Mv density biomass
ggplot(mvD2Dat, aes(Mv_seedling_density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# 2019 site variation
mvD2Dat %>%
  mutate(plot1 = ifelse(plot == 1, "1", "other")) %>%
  ggplot(aes(site, biomass_weight.g, color = treatment, shape = plot1)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")


#### Mv density model ####

# remove plot 1 (large environmental effect)
mvD2Dat2 <- mvD2Dat %>%
  filter(plot != 1)

# split by species (couldn't specify appropriate priors with them all together)
mvMvD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Mv seedling")
mvEvSD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Ev seedling")
mvEvAD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Ev adult")

# intercept prior
mvD2Dat2 %>%
  filter(Mv_seedling_density == 7) %>%
  group_by(pg_trt) %>%
  summarise(bio = mean(biomass_weight.g))

x <- seq(-1, 25, length.out = 100)
y <- dgamma(x, shape = 7, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# alpha prior
x <- seq(0, 10, length.out = 100)
y <- dexp(x, 0.5)
plot(x, y, type = "l")

# fit models
mvMvD2Mod <- brm(data = mvMvD2Dat, family = gaussian,
               bf(biomass_weight.g ~ b0/(1 + alpha * Mv_seedling_density), 
                  b0 ~ 0 + treatment, 
                  alpha ~ 0 + treatment, 
                  nl = T),
               prior <- c(prior(gamma(20, 1), nlpar = "b0", lb = 0),
                          prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma
               iter = 6000, warmup = 1000, chains = 3) 

summary(mvMvD2Mod)
pp_check(mvMvD2Mod, nsamples = 100)

mvEvSD2Mod <- update(mvMvD2Mod, newdata = mvEvSD2Dat,
                     prior = set_prior("gamma(2.5, 1)", nlpar = "b0", lb = 0))

prior_summary(mvEvSD2Mod)
summary(mvEvSD2Mod)
pp_check(mvEvSD2Mod, nsamples = 100)

mvEvAD2Mod <- update(mvMvD2Mod, newdata = mvEvAD2Dat,
                     prior = set_prior("gamma(12, 1)", nlpar = "b0", lb = 0))

prior_summary(mvEvAD2Mod)
summary(mvEvAD2Mod)
pp_check(mvEvAD2Mod, nsamples = 100)

# simulated data
mvMvD2Sim <- mvMvD2Dat %>%
  select(treatment) %>%
  unique() %>%
  expand_grid(tibble(Mv_seedling_density = seq(0, 67, length.out = 100))) %>%
  mutate(pred = fitted(mvMvD2Mod, newdata = ., re_formula = NA)[, "Estimate"],
         lower = fitted(mvMvD2Mod, newdata = ., re_formula = NA)[, "Q2.5"],
         upper = fitted(mvMvD2Mod, newdata = ., re_formula = NA)[, "Q97.5"])

mvEvSD2Sim <- mvEvSD2Dat %>%
  select(treatment) %>%
  unique() %>%
  expand_grid(tibble(Mv_seedling_density = seq(0, 67, length.out = 100))) %>%
  mutate(pred = fitted(mvEvSD2Mod, newdata = ., re_formula = NA)[, "Estimate"],
         lower = fitted(mvEvSD2Mod, newdata = ., re_formula = NA)[, "Q2.5"],
         upper = fitted(mvEvSD2Mod, newdata = ., re_formula = NA)[, "Q97.5"])

mvEvAD2Sim <- mvEvAD2Dat %>%
  select(treatment) %>%
  unique() %>%
  expand_grid(tibble(Mv_seedling_density = seq(0, 67, length.out = 100))) %>%
  mutate(pred = fitted(mvEvAD2Mod, newdata = ., re_formula = NA)[, "Estimate"],
         lower = fitted(mvEvAD2Mod, newdata = ., re_formula = NA)[, "Q2.5"],
         upper = fitted(mvEvAD2Mod, newdata = ., re_formula = NA)[, "Q97.5"])

# model fit
ggplot(mvMvD2Dat, aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvMvD2Sim, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line(data = mvMvD2Sim, aes(y = pred)) +
  theme_bw()

ggplot(mvEvSD2Dat, aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvEvSD2Sim, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line(data = mvEvSD2Sim, aes(y = pred)) +
  theme_bw()

ggplot(mvEvAD2Dat, aes(x = Mv_seedling_density, y = biomass_weight.g, color = treatment, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvEvAD2Sim, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line(data = mvEvAD2Sim, aes(y = pred)) +
  theme_bw()


#### figure ####

# combine predicted datasets
mvD2Sim <- mvMvD2Sim %>%
  mutate(plant_group = "Mv seedling") %>%
  full_join(mvEvSD2Sim %>%
              mutate(plant_group = "Ev seedling")) %>%
  full_join(mvEvAD2Sim %>%
              mutate(plant_group = "Ev adult")) %>%
  mutate(plant_group = fct_rev(plant_group))

# use fungicide treatment
mvFunD2Sim <- mvD2Sim %>%
  filter(treatment == "fungicide")

mvFunD2Dat <- mvD2Dat %>%
  filter(treatment == "fungicide")

ggplot(mvFunD2Dat, aes(x = Mv_seedling_density, y = biomass_weight.g, color = plant_group, fill = plant_group)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_ribbon(data = mvFunD2Sim, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line(data = mvFunD2Sim, aes(y = pred)) +
  theme_bw()


#### output ####
save(mvMvD2Mod, file = "output/mv_biomass_mv_density_model_2019_density_exp.rda")
save(mvEvSD2Mod, file = "output/evS_biomass_mv_density_model_2019_density_exp.rda")
save(mvEvAD2Mod, file = "output/evA_biomass_mv_density_model_2019_density_exp.rda")
##### info ####

# file: mv_plot_biomass_seeds_density_2019_density_exp
# author: Amy Kendig
# date last edited: 11/24/21
# goal: mv plot-scale biomass and seeds (i.e., abundance) vs. density
# model code from plot_biomass_density_2019_dens_exp.R


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import plot information
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import plot biomass and seed data
plotD2Dat <- read_csv("intermediate-data/mv_plot_biomass_seeds_2019_density_exp.csv")
# plot_data_processing_2019_density_exp

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# plant group densities
plotDens <- plots %>%
  select(plot, treatment, background, background_density)

# combine data
# add focal Mv to background density
plotD2Dat2 <- plotD2Dat %>%
  left_join(plotDens) %>%
  filter(background %in% c("none", "Mv seedling")) %>%
  mutate(density = background_density + 3) %>%
  rename(biomass = biomass_mv)


#### fit model ####

# initial visualization
ggplot(plotD2Dat2, aes(x = density, y = biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2))

ggplot(plotD2Dat2, aes(x = density, y = seeds, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2))

# priors
plotD2Dat2 %>%
  filter(plot == 1) %>%
  group_by(treatment) %>%
  summarise(b0 = mean(biomass/3),
            s0 = mean(seeds/3))

x <- seq(-1, 20, length.out = 100)
y <- dgamma(x, shape = 14, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

x <- seq(0, 10, length.out = 100)
y <- dexp(x, 0.5)
plot(x, y, type = "l")

x <- seq(800, 1300, length.out = 100)
y <- dnorm(x, mean = 1050, sd = 100) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# model
mvBioDensMod <- brm(data = plotD2Dat2, family = gaussian,
                    bf(biomass ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(14, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(mvBioDensMod)

mvSeedDensMod <- brm(data = plotD2Dat2, family = gaussian,
                    bf(seeds ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(normal(1050, 100), nlpar = "b0"),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(mvSeedDensMod)


#### treatment effects

b0_eff = "b0_treatmentfungicide - b0_treatmentwater = 0"
alpha_eff = "alpha_treatmentfungicide - alpha_treatmentwater = 0"

mvBioEff <- hypothesis(mvBioDensMod, c(b0_eff, alpha_eff))
mvSeedEff <- hypothesis(mvSeedDensMod, c(b0_eff, alpha_eff))

mvEff <- mvBioEff[[1]] %>%
  mutate(response = "biomass") %>%
  full_join(mvSeedEff[[1]] %>%
              mutate(response = "seeds")) %>%
  mutate(parameter = rep(c("b0", "alpha"), 2)) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob)) %>%
  relocate(response, parameter)

# save outputs
save(mvBioDensMod, file = "output/mv_plot_biomass_density_model_2019_dens_exp.rda")
save(mvSeedDensMod, file = "output/mv_plot_seed_density_model_2019_dens_exp.rda")
write_csv(mvEff, "output/mv_plot_biomass_seeds_density_treatment_effect_2019_dens_exp.csv")

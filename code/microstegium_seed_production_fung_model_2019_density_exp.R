##### info ####

# file: microstegium_seed_production_fung_model_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: estimate Microstegium seed production based on fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
seedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
# adj values substitute averages of other two plants in the plot when a plant is missing data
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# add columns
# remove missing data
mvSeedD2Dat <- seedD2Dat %>%
  left_join(plotDens) %>%
  mutate(log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot)) %>%
  filter(!is.na(seeds))


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- mvSeedD2Dat %>%
  select(-c(background, background_sp)) %>%
  left_join(plotsD2)

# non-transformed seeds
ggplot(vizDat, aes(background_density, seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(~background, scales = "free") +
  theme_bw()

# log-transformed seeds
ggplot(vizDat, aes(background_density, log_seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(~background, scales = "free") +
  theme_bw()


#### fit regression ####

# select plots with no background (negative effect on growth)
mvSeedD2Dat0 <- mvSeedD2Dat %>%
  filter(background == "none")

# initial fit
mvSeedFuD2Mod1 <- brm(log_seeds ~ fungicide + (1|site),
                  data = mvSeedD2Dat0, family = gaussian,
                  prior <- c(prior(normal(6.5, 3), class = "Intercept"),
                             prior(normal(0, 10), class = "b"),
                             prior(cauchy(0, 1), class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 8 divergent transitions
summary(mvSeedFuD2Mod1)

# increase chains
mvSeedFuD2Mod2 <- update(mvSeedFuD2Mod1, chains = 3,
                       control = list(adapt_delta = 0.999))
summary(mvSeedFuD2Mod2)
plot(mvSeedFuD2Mod2)
pp_check(mvSeedFuD2Mod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(fungicide = c(1, 0)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))

# simulate fit
fitDat <- simDat %>%
  mutate(log_seeds =fitted(mvSeedFuD2Mod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_seeds_lower = fitted(mvSeedFuD2Mod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_seeds_upper = fitted(mvSeedFuD2Mod2, newdata = ., re_formula = NA)[, "Q97.5"],
         seeds = exp(log_seeds) - 1,
         seeds_lower = exp(log_seeds_lower) - 1,
         seeds_upper = exp(log_seeds_upper) - 1)

# summarise by site
vizDat2 <- mvSeedD2Dat0  %>%
  group_by(site, treatment) %>%
  summarise(seeds = mean(seeds),
            log_seeds = mean(log_seeds))

# fit figure
ggplot(fitDat, aes(treatment, seeds, color = treatment)) +
  geom_errorbar(aes(ymin = seeds_lower, ymax = seeds_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = vizDat2, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
  theme_bw()

(mvSeedFuD2Plot <- ggplot(fitDat, aes(treatment, log_seeds, color = treatment)) +
    geom_errorbar(aes(ymin = log_seeds_lower, ymax = log_seeds_upper), width = 0) +
    geom_point(size = 4, shape = 16) +
    geom_point(data = vizDat2, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
    theme_bw())



#### output ####
save(mvSeedFuD2Mod2, file = "output/microstegium_seed_fung_model_2019_density_exp.rda")
save(mvSeedFuD2Plot, file = "output/microstegium_seed_fung_figure_2019_density_exp.rda")

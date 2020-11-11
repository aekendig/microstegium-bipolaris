##### info ####

# file: elymus_adult_seed_production_fung_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/10/20
# goal: estimate Elymus adult seed production based on fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
seedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
seedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# summer survival 2018
# make survival 1 if the plant produced seeds in summer
# remove NA's 
evAGsSurvD1Dat <- survD1Dat %>%
  filter(month == "September" & sp == "Ev" & focal == 1 & age == "adult") %>%
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival))

# add columns
# remove missing data
evASeedD1Dat <- seedD1Dat %>%
  filter(ID == "A" & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evAGsSurvD1Dat %>%
               filter(survival == 1) %>%
               select(site, plot, treatment, sp, age, ID, survival)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 1")

# check for notes about seeds
unique(seedD2Dat$spikelet_notes)
filter(seedD2Dat, is.na(spikelet_weight.g))
# no missing data

# add columns
# remove missing data
evASeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         age = ifelse(ID == "A", "adult", "seedling"),
         yearf = "year 2") %>%
  filter(age == "adult")

# combine years
evASeedDat <- full_join(evASeedD1Dat, evASeedD2Dat)


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- evASeedDat %>%
  select(-c(background, background_sp)) %>%
  left_join(plotsD2)

# non-transformed seeds
ggplot(vizDat, aes(background_density, seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()

# log-transformed seeds
ggplot(vizDat, aes(background_density, log_seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()


#### fit regression ####

# select plots with no background (negative effect on growth)
evASeedDat0 <- evASeedDat %>%
  filter(background == "none")

# initial fit
evASeedFuMod1 <- brm(log_seeds ~ fungicide * yearf + (1|site),
                  data = evASeedDat0, family = gaussian,
                  prior <- c(prior(normal(4, 3), class = "Intercept"),
                             prior(normal(0, 10), class = "b"),
                             prior(cauchy(0, 1), class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 2 divergent transitions
summary(evASeedFuMod1)

# increase chains and adapt delta
evASeedFuMod2 <- update(evASeedFuMod1, chains = 3,
                     control = list(adapt_delta = 0.99999))
summary(evASeedFuMod2)
plot(evASeedFuMod2)
pp_check(evASeedFuMod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(fungicide = rep(c(1, 0), 2),
                 yearf = rep(c("year 1", "year 2"), each = 2)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))

# simulate fit
fitDat <- simDat %>%
  mutate(log_seeds =fitted(evASeedFuMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_seeds_lower = fitted(evASeedFuMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_seeds_upper = fitted(evASeedFuMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         seeds = exp(log_seeds) - 1,
         seeds_lower = exp(log_seeds_lower) - 1,
         seeds_upper = exp(log_seeds_upper) - 1)

# summarise by site
vizDat2 <- evASeedDat0  %>%
  group_by(yearf, site, treatment) %>%
  summarise(seeds = mean(seeds),
            log_seeds = mean(log_seeds))

# fit figure
ggplot(fitDat, aes(treatment, seeds, color = treatment)) +
  geom_errorbar(aes(ymin = seeds_lower, ymax = seeds_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = vizDat2, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
  facet_wrap(~ yearf) +
  theme_bw()

(evASeedFuPlot <- ggplot(fitDat, aes(treatment, log_seeds, color = treatment)) +
    geom_errorbar(aes(ymin = log_seeds_lower, ymax = log_seeds_upper), width = 0) +
    geom_point(size = 4, shape = 16) +
    geom_point(data = vizDat2, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
    facet_wrap(~ yearf) +
    theme_bw())



#### output ####
save(evASeedFuMod2, file = "output/elymus_adult_seed_fung_model_2018_2019_density_exp.rda")
save(evASeedFuPlot, file = "output/elymus_adult_seed_fung_figure_2018_2019_density_exp.rda")

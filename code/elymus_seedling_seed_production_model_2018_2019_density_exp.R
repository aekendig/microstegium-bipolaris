##### info ####

# file: elymus_seedling_seed_production_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: estimate Elymus seedling seed production based on density and fungicide treatments


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
plotsD1 <- read_csv("intermediate-data/plot_densities_2018_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens1 <- plotsD1 %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

plotDens2 <- plotsD2 %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# summer survival 2018
# make survival 1 if the plant produced seeds in summer
# remove NA's 
evSGsSurvD1Dat <- survD1Dat %>%
  filter(month == "September" & sp == "Ev" & focal == 1 & age == "seedling") %>%
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival))

# add columns
# remove missing data
evSSeedD1Dat <- seedD1Dat %>%
  filter(ID %in% c("1", "2", "3") & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evSGsSurvD1Dat %>%
               filter(survival == 1) %>%
               select(site, plot, treatment, sp, age, ID, survival)) %>%
  left_join(plotDens1) %>%
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
evSSeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  left_join(plotDens2) %>%
  mutate(seeds = replace_na(seeds, 0),
         log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         age = ifelse(ID == "A", "adult", "seedling"),
         yearf = "year 2") %>%
  filter(age == "seedling")

# combine years
evSSeedDat <- full_join(evSSeedD1Dat, evSSeedD2Dat)


#### initial visualizations ####

# non-transformed seeds
ggplot(evSSeedDat, aes(background_density, seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()

# log-transformed seeds
ggplot(evSSeedDat, aes(background_density, log_seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()


#### fit regression ####

# remove plots with no background (negative effect on growth)
evSSeedDat2 <- evSSeedDat %>%
  mutate(treatment = ifelse(treatment == "water", "control", treatment),
         treatment_year = paste(treatment, yearf, sep = "_")) %>%
  filter(background != "none")

# initial fit
evSSeedMod1 <- brm(bf(log_seeds ~ logS - log(1 + alphaA * mv_seedling_density + alphaS * ev_seedling_density + alphaP * ev_adult_density),
                        logS ~ yearf * treatment + (1|site),
                        alphaA ~ 0 + treatment_year,
                        alphaS ~ 0 + treatment_year,
                        alphaP ~ 0 + treatment_year,
                        nl = T),
                     data = evSSeedDat2, family = gaussian,
                     prior <- c(prior(normal(2, 10), nlpar = "logS", class = "b", coef = "Intercept"),
                                prior(normal(0, 10), nlpar = "logS", class = "b"),
                                prior(exponential(0.5), nlpar = "alphaA", lb = 0),
                                prior(exponential(0.5), nlpar = "alphaS", lb = 0),
                                prior(exponential(0.5), nlpar = "alphaP", lb = 0),
                                prior(cauchy(0, 1), nlpar = "logS", class = "sd"),
                                prior(cauchy(0, 1), class = "sigma")),
                     iter = 6000, warmup = 1000, chains = 1)
# 2 divergent transition
summary(evSSeedMod1)

# increase chains and adapt delta
evSSeedMod2 <- update(evSSeedMod1, chains = 3,
                      control = list(adapt_delta = 0.99999,
                                     max_treedepth = 15))
summary(evSSeedMod2)
plot(evSSeedMod2)
pp_check(evSSeedMod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(mv_seedling_density = c(seq(0, 64, length.out = 100), rep(0, 200)),
                 ev_seedling_density = c(rep(0, 100), seq(0, 16, length.out = 100), rep(0, 100)),
                 ev_adult_density = c(rep(0, 200), seq(0, 8, length.out = 100)),
                 background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 100)) %>%
  expand_grid(tibble(fungicide = rep(c(1, 0), 2),
                 yearf = rep(c("year 1", "year 2"), each = 2))) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         background_density = case_when(background == "Mv seedling" ~ mv_seedling_density,
                                        background == "Ev seedling" ~ ev_seedling_density,
                                        TRUE ~ ev_adult_density),
         treatment_year = paste(treatment, yearf, sep = "_"))

# simulate fit
fitDat <- simDat %>%
  mutate(log_seeds = fitted(evSSeedMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_seeds_lower = fitted(evSSeedMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_seeds_upper = fitted(evSSeedMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
evSSeedDat3 <- evSSeedDat2 %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment))

# fit figure
(evSSeedPlot <- ggplot(fitDat, aes(background_density, log_seeds, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = log_seeds_lower, ymax = log_seeds_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = evSSeedDat3, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = evSSeedDat3, geom = "point", size = 2, fun = "mean") +
  facet_grid(yearf ~ background, scales = "free_x") +
  theme_bw())



#### output ####
save(evSSeedMod2, file = "output/elymus_seedling_seed_model_2018_2019_density_exp.rda")
save(evSSeedPlot, file = "output/elymus_seedling_seed_figure_2018_2019_density_exp.rda")

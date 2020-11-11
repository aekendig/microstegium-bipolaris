##### info ####

# file: microstegium_seed_production_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/10/20
# goal: estimate Microstegium seed production based on density and fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
seedD1Dat <- read_csv("./intermediate-data/mv_processed_seeds_2018_density_exp.csv")
seedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
# adj values substitute averages of other two plants in the plot when a plant is missing data
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# check notes
unique(seedD2Dat$process_notes)
unique(seedD2Dat$analysis_notes)

# missing data?
filter(seedD2Dat, is.na(seeds)) %>%
  data.frame()
# 3 plants

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# add columns
# remove missing data
mvSeedD1Dat <- seedD1Dat %>%
  left_join(plotDens) %>%
  mutate(seeds = bag_seeds,
         log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 1")

# add columns
# remove missing data
mvSeedD2Dat <- seedD2Dat %>%
  left_join(plotDens) %>%
  mutate(log_seeds = log(seeds),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 2") %>%
  filter(!is.na(seeds))

# combine years
mvSeedDat <- full_join(mvSeedD1Dat, mvSeedD2Dat) %>%
  mutate(seeds_100 = seeds / 100)


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- mvSeedDat %>%
  select(-c(background, background_sp)) %>%
  left_join(plotsD2)

# non-transformed seeds
ggplot(vizDat, aes(background_density, seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()

ggplot(vizDat, aes(background_density, seeds_100, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_grid(yearf~background, scales = "free") +
  theme_bw()


#### fit regression ####

# remove plots with no background (negative effect on growth)
mvSeedDat2 <- mvSeedDat %>%
  filter(background != "none") %>%
  mutate(treatment = ifelse(treatment == "water", "control", treatment),
         treatment_year = paste(treatment, yearf, sep = "_"))

# initial fit
# using seeds_100 because the raw seed values are very large and cause many divergent transitions
mvSeedMod1 <- brm(bf(seeds_100 ~ maxS/(1 + gammaA * mv_seedling_density + gammaS * ev_seedling_density + gammaP * ev_adult_density),
                     maxS ~ yearf * treatment + (1|site/plotr),
                     gammaA ~ 0 + treatment_year,
                     gammaS ~ 0 + treatment_year,
                     gammaP ~ 0 + treatment_year,
                     nl = T),
                  data = mvSeedDat2, family = gaussian,
                  prior <- c(prior(normal(8, 10), nlpar = "maxS", class = "b", coef = "Intercept"),
                             prior(normal(10, 10), nlpar = "maxS", class = "b", coef = "yearfyear2"),
                             prior(normal(0, 10), nlpar = "maxS", class = "b"),
                             prior(exponential(0.5), nlpar = "gammaA", lb = 0),
                             prior(exponential(0.5), nlpar = "gammaS", lb = 0),
                             prior(exponential(0.5), nlpar = "gammaP", lb = 0),
                             prior(cauchy(0, 1), nlpar = "maxS", class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 19 divergent transitions
summary(mvSeedMod1)

# increase chains and adapt delta
mvSeedMod2 <- update(mvSeedMod1, chains = 3,
                     control = list(adapt_delta = 0.999))
summary(mvSeedMod2)
plot(mvSeedMod2)
pp_check(mvSeedMod2, nsamples = 50)


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
  mutate(seeds = fitted(mvSeedMod2, newdata = ., re_formula = NA)[, "Estimate"] * 100,
         seeds_lower = fitted(mvSeedMod2, newdata = ., re_formula = NA)[, "Q2.5"] * 100,
         seeds_upper = fitted(mvSeedMod2, newdata = ., re_formula = NA)[, "Q97.5"] * 100,
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
mvSeedDat3 <- mvSeedDat2 %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment))

# fit figure
(mvSeedPlot <- ggplot(fitDat, aes(background_density, seeds, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = seeds_lower, ymax = seeds_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = mvSeedDat3, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = mvSeedDat3, geom = "point", size = 2, fun = "mean") +
  facet_grid(yearf ~ background, scales = "free_x") +
  theme_bw())


#### output ####
save(mvSeedMod2, file = "output/microstegium_seed_model_2018_2019_density_exp.rda")
save(mvSeedPlot, file = "output/microstegium_seed_figure_2018_2019_density_exp.rda")

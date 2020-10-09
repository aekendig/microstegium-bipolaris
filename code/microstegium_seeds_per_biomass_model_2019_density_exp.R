##### info ####

# file: microstegium_seeds_per_biomass_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/8/20
# goal: analyze Microstegium seeds per unit biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
seedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
# adj values substitute averages of other two plants in the plot when a plant is missing data


#### edit data ####

# add columns
# remove missing data
mvSeedD2Dat <- seedD2Dat %>%
  rename(bio.g = biomass_weight.g) %>%
  mutate(seeds = replace_na(seeds, 0),
         seeds_per_bio = seeds / bio.g,
         log_seeds_per_bio = log(seeds_per_bio),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot)) %>%
  filter(!is.na(seeds_per_bio))


#### initial visualizations ####

# figure
ggplot(mvSeedD2Dat, aes(log_bio.g, log_seeds_per_bio, color = treatment)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw()
# not enough zeros to fit a separate regression

filter(mvSeedD2Dat, seeds_per_bio == 0)
# 3


#### normal regression ####

# subset data for no zeros
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  filter(seeds_per_bio > 0)

# initial fit
mvSeedBioD2Mod1 <- brm(log_seeds_per_bio ~ fungicide * log_bio.g + (1|site/plotr),
                       data = mvSeedD2Dat2, family = gaussian,
                       prior <- c(prior(normal(0, 10), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sigma)),
                       iter = 6000, warmup = 1000, chains = 1)
# 57 divergent transitions
summary(mvSeedBioD2Mod1)

# increase chains and adapt delta
mvSeedBioD2Mod2 <- update(mvSeedBioD2Mod1, chains = 3,
                          control = list(adapt_delta = 0.99999, max_treedepth = 15))
summary(mvSeedBioD2Mod2)
plot(mvSeedBioD2Mod2)
pp_check(mvSeedBioD2Mod2, nsamples = 50)


#### visualize normal regression ####

# simulate data
simDat <- tibble(log_bio.g = seq(min(mvSeedD2Dat2$log_bio.g), max(mvSeedD2Dat2$log_bio.g), length.out = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         site = NA,
         plotr = NA)

# simulate fit
fitDat <- simDat %>%
  mutate(log_seeds_per_bio = fitted(mvSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         log_seeds_per_bio_lower = fitted(mvSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         log_seeds_per_bio_upper = fitted(mvSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# fit figure
ggplot(mvSeedD2Dat2, aes(log_bio.g, log_seeds_per_bio, color = treatment, fill = treatment)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_ribbon(data = fitDat, aes(ymin = log_seeds_per_bio_lower, ymax = log_seeds_per_bio_upper), alpha = 0.5, color = NA) +
  geom_line(data = fitDat, size = 1.5) +
  theme_bw()
# seeds slightly decrease with biomass, no fungicide effect


#### output ####

save(mvSeedBioD2Mod2, file = "output/microstegium_seeds_per_biomass_model_2019_density_exp.rda")

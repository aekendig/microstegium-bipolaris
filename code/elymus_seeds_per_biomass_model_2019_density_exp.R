##### info ####

# file: elymus_seeds_per_biomass_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/8/20
# goal: analyze Elymus seeds per unit biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
seedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")


#### edit data ####

# add columns
# remove missing data
evSeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              rename(bio.g = weight)) %>%
  mutate(seeds = replace_na(seeds, 0),
         seeds_per_bio = seeds / bio.g,
         seeds_prod = ifelse(seeds > 0, 1, 0),
         log_seeds_per_bio = log(seeds_per_bio),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  filter(!is.na(seeds_per_bio))


#### initial visualizations ####

# figure
ggplot(evSeedD2Dat, aes(log_bio.g, log_seeds_per_bio, color = treatment, shape = age)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw()


# histogram
ggplot(evSeedD2Dat, aes(seeds_per_bio)) +
  geom_histogram() + 
  theme_bw()
# zero inflated

# yes/no seeds
ggplot(evSeedD2Dat, aes(log_bio.g, seeds_prod, color = treatment, shape = age)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw()
# a minimum biomass needed for seeds to be produced


#### logistic regression ####

# initial fit
evSeedD2Mod1 <- brm(seeds_prod ~ fungicide * log_bio.g + (1|site/plotr),
                    data = evSeedD2Dat, family = bernoulli,
                    prior = c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b),
                              prior(cauchy(0, 1), class = sd)),
                    iter = 6000, warmup = 1000, chains = 1)
# 7 divergent transitions
summary(evSeedD2Mod1)

# increase chains and adapt delta
evSeedD2Mod2 <- update(evSeedD2Mod1, chains = 3,
                        control = list(adapt_delta = 0.999))
summary(evSeedD2Mod2)
plot(evSeedD2Mod2)
pp_check(evSeedD2Mod2, nsamples = 50)


#### visualize logistic regression ####

# simulation data
simDat <- tibble(log_bio.g = seq(min(evSeedD2Dat$log_bio.g), max(evSeedD2Dat$log_bio.g), length.out = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         site = NA,
         plotr = NA)

# simulate fit
fitDat <- simDat %>%
  mutate(seeds_prod = fitted(evSeedD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         seeds_prod_lower = fitted(evSeedD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         seeds_prod_upper = fitted(evSeedD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# fit figure
ggplot(evSeedD2Dat, aes(log_bio.g, seeds_prod, color = treatment, fill = treatment)) +
  geom_point(size = 2, alpha = 0.5, aes(shape = age)) +
  geom_ribbon(data = fitDat, aes(ymin = seeds_prod_lower, ymax = seeds_prod_upper), alpha = 0.5, color = NA) +
  geom_line(data = fitDat, size = 1.5) +
  theme_bw()


#### normal regression ####

# subset data for no zeros
evSeedD2Dat2 <- evSeedD2Dat %>%
  filter(seeds_per_bio > 0)

# initial fit
evSeedBioD2Mod1 <- brm(log_seeds_per_bio ~ fungicide * log_bio.g + (1|site/plotr),
                        data = evSeedD2Dat2, family = gaussian,
                        prior <- c(prior(normal(0, 10), class = Intercept),
                                   prior(normal(0, 10), class = b),
                                   prior(cauchy(0, 1), class = sigma)),
                        iter = 6000, warmup = 1000, chains = 1)
# 38 divergent transitions
summary(evSeedBioD2Mod1)

# increase chains and adapt delta
evSeedBioD2Mod2 <- update(evSeedBioD2Mod1, chains = 3,
                          control = list(adapt_delta = 0.999))
summary(evSeedBioD2Mod2)
plot(evSeedBioD2Mod2)
pp_check(evSeedBioD2Mod2, nsamples = 50)


#### visualize normal regression ####

# simulate data
simDat2 <- tibble(log_bio.g = seq(min(evSeedD2Dat2$log_bio.g), max(evSeedD2Dat2$log_bio.g), length.out = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         site = NA,
         plotr = NA)

# simulate fit
fitDat2 <- simDat2 %>%
  mutate(log_seeds_per_bio = fitted(evSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         log_seeds_per_bio_lower = fitted(evSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         log_seeds_per_bio_upper = fitted(evSeedBioD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# fit figure
ggplot(evSeedD2Dat2, aes(log_bio.g, log_seeds_per_bio, color = treatment, fill = treatment)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_ribbon(data = fitDat2, aes(ymin = log_seeds_per_bio_lower, ymax = log_seeds_per_bio_upper), alpha = 0.5, color = NA) +
  geom_line(data = fitDat2, size = 1.5) +
  theme_bw()
# seeds slightly increase with biomass, slightly stronger with fungicide, but barely


#### output ####

save(evSeedD2Mod2, file = "output/elymus_seeds_produced_model_2019_density_exp.rda")
save(evSeedBioD2Mod2, file = "output/elymus_seeds_per_biomass_model_2019_density_exp.rda")

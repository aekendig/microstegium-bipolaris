##### info ####

# file: elymus_seedling_ngs_survival_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: estimate Elymus seedling non-growing season survival based on density and fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# winter survival 2018-2019
evSNgsSurvD1Dat <- survD1Dat %>%
  filter((month == "April" | month == "September") & focal == 1 & sp == "Ev" & age == "seedling") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot)) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(plotDens)

# average survival
mean(evSNgsSurvD1Dat$survival)


#### fit regression ####

# Initially fit regression with plant density information included. Removed this because (1) this is survival during the time when plants are not taking up resources and (2) the densities were from 2018 when they were not as well maintained. See Github archive for code.

# initial fit
evSNgsSurvD1Mod1 <- brm(survival ~ fungicide + (1|site/plotr),
                        data = evSNgsSurvD1Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
# 22 divergent transitions
summary(evSNgsSurvD1Mod1)

# increase chains
evSNgsSurvD1Mod2 <- update(evSNgsSurvD1Mod1, chains = 3,
                           control = list(adapt_delta = 0.999))
summary(evSNgsSurvD1Mod2)
plot(evSNgsSurvD1Mod2)
pp_check(evSNgsSurvD1Mod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(1, 0)) %>%
  mutate(survival = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))

# modify dataset so that zero plots are repeated
vizDat <- evSNgsSurvD1Dat %>%
  group_by(site, treatment) %>%
  summarise(survival = mean(survival))

# fit figure
evSNgsSurvD1Plot <- ggplot(fitDat, aes(treatment, survival, color = treatment)) +
  geom_errorbar(aes(ymin = survival_lower, ymax = survival_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  geom_point(data = vizDat, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  theme_bw()



#### output ####

save(evSNgsSurvD1Mod2, file = "output/elymus_seedling_ngs_survival_model_2019_density_exp.rda")
save(evSNgsSurvD1Plot, file = "output/elymus_seedling_ngs_survival_figure_2019_density_exp.rda")

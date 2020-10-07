##### info ####

# file: elymus_adult_ngs_survival_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/6/20
# goal: estimate Elymus adult non-growing season survival based on density and fungicide treatments


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
evANgsSurvD1Dat <- survD1Dat %>%
  filter((month == "April" | month == "September") & focal == 1 & sp == "Ev" & age == "adult") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(plotDens)

# average survival
mean(evANgsSurvD1Dat$survival)
# very high


#### fit regression ####

# initial fit
evANgsSurvD1Mod1 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site),
                      data = evANgsSurvD1Dat,
                      family = bernoulli,
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(normal(0, 10), class = b),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)
# 5 divergent transitions and several other warnings
# likely because there are too few zeros
summary(evANgsSurvD1Mod1)

# intercept-only model
evANgsSurvD1Mod2 <- brm(survival ~ 1 + (1|site),
                        data = evANgsSurvD1Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 3,
                        control = list(adapt_delta = 0.999))
summary(evANgsSurvD1Mod2)
plot(evANgsSurvD1Mod2)
pp_check(evANgsSurvD1Mod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(site = c("D1", "D2", "D3", "D4")) %>%
  mutate(survival = fitted(evANgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(evANgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(evANgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# average b
vizDat <- evANgsSurvD1Dat %>%
  group_by(site) %>%
  summarise(survival_obs = mean(survival),
            survival_lower_obs = mean_cl_boot(survival)$ymin,
            survival_upper_obs = mean_cl_boot(survival)$ymax) %>%
  full_join(fitDat)

# fit figure
ggplot(vizDat, aes(survival_obs, survival, color = site)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_errorbar(aes(ymin = survival_lower, ymax = survival_upper), width = 0) +
  geom_errorbarh(aes(xmin = survival_lower_obs, xmax = survival_upper_obs), height = 0) +
  geom_point(size = 2) +
  theme_bw()
# density tends to increase survival without disease and decrease survival with disease
# fits are messy


#### output ####

save(evANgsSurvD1Mod2, file = "output/elymus_adult_ngs_survival_model_2019_density_exp.rda")

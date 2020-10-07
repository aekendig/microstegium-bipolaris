##### info ####

# file: elymus_seedling_ngs_survival_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/6/20
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
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April) %>%
  left_join(plotDens)

# average survival
mean(evSNgsSurvD1Dat$survival)


#### fit regression ####

# initial fit
evSNgsSurvD1Mod1 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site/plot),
                      data = evSNgsSurvD1Dat,
                      family = bernoulli,
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(normal(0, 10), class = b),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)
# 13 divergent transitions
summary(evSNgsSurvD1Mod1)

# increase chains
evSNgsSurvD1Mod2 <- update(evSNgsSurvD1Mod1, chains = 3,
                           control = list(adapt_delta = 0.999))
summary(evSNgsSurvD1Mod2)
plot(evSNgsSurvD1Mod2)
pp_check(evSNgsSurvD1Mod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(mv_seedling_density = c(seq(0, 64, length.out = 100), rep(0, 200)),
                 ev_seedling_density = c(rep(0, 100), seq(0, 16, length.out = 100), rep(0, 100)),
                 ev_adult_density = c(rep(0, 200), seq(0, 8, length.out = 100)),
                 background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"),
         site = NA,
         background_density = case_when(background == "Mv seedling" ~ mv_seedling_density,
                                        background == "Ev seedling" ~ ev_seedling_density,
                                        TRUE ~ ev_adult_density))

# simulate fit
fitDat <- simDat %>%
  mutate(survival = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(evSNgsSurvD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# modify dataset so that zero plots are repeated
vizDat <- evSNgsSurvD1Dat %>%
  select(-c(background, background_sp)) %>%
  full_join(plotsD2)

# fit figure
ggplot(fitDat, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = vizDat, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = vizDat, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()
# density tends to increase survival without disease and decrease survival with disease
# fits are messy


#### output ####

save(evSNgsSurvD1Mod2, file = "output/elymus_seedling_ngs_survival_model_2019_density_exp.rda")

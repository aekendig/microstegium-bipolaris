##### info ####

# file: elymus_seedling_gs_survival_bh_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/4/20
# goal: estimate Elymus seedling growing season survival average for intercept of Beverton-Holt model


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# summer survival 2018
# make survival 1 if the plant produced seeds in summer
# remove NA's 
evSGsSurvD1Dat <- survD1Dat %>%
  filter(month == "September" & sp == "Ev" & focal == 1 & age == "seedling") %>%
  select(-month) %>%
  left_join(plotDens) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 1",
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(!is.na(survival))

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3)))

# summer survival 2019
# use the latest replacement date
evSGsSurvD2Dat <- survD2Dat %>%
  filter(focal == 1 & sp == "Ev" & age == "seedling") %>%
  group_by(site, plot, treatment, ID) %>%
  summarise(replace_date = max(replace_date),
            survival = 0) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plotDens) %>%
  mutate(survival = replace_na(survival, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 2")

# combine years
evSGsSurvDat <- full_join(evSGsSurvD1Dat, evSGsSurvD2Dat)


#### fit regression ####

# initial fit
evSGsSurvBhMod1 <- brm(survival ~ fungicide * yearf + (1|site/plotr),
                      data = evSGsSurvDat,
                      family = bernoulli,
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(normal(0, 10), class = b),
                                prior(cauchy(0, 1), class = sd)),
                      iter = 6000, warmup = 1000, chains = 1)
summary(evSGsSurvBhMod1)

# increase chains
evSGsSurvBhMod2 <- update(evSGsSurvBhMod1, chains = 3,
                          control = list(adapt_delta = 0.999))
summary(evSGsSurvBhMod2)
plot(evSGsSurvBhMod2)
pp_check(evSGsSurvBhMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = rep(c(1, 0), 2),
                 yearf = rep(c("year 1", "year 2"), each = 2)) %>%
  mutate(survival = fitted(evSGsSurvBhMod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(evSGsSurvBhMod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(evSGsSurvBhMod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))

# summarise by site
vizDat <- evSGsSurvDat  %>%
  group_by(yearf, site, treatment) %>%
  summarise(survival = mean(survival))

# fit figure
(evSGsSurvBhPlot <- ggplot(fitDat, aes(treatment, survival, color = treatment)) +
    geom_errorbar(aes(ymin = survival_lower, ymax = survival_upper), width = 0) +
    geom_point(size = 4, shape = 16) +
    geom_point(data = vizDat, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
    facet_wrap(~ yearf) +
    theme_bw())


#### output ####

save(evSGsSurvBhMod2, file = "output/elymus_seedling_gs_survival_bh_model_2018_2019_density_exp.rda")
save(evSGsSurvBhPlot, file = "output/elymus_seedling_gs_survival_bh_figure_2018_2019_density_exp.rda")

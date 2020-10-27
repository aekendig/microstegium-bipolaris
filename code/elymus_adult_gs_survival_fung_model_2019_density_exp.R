##### info ####

# file: elymus_adult_gs_survival_fung_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate Elymus adult growing season survival based on fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  mutate(ID = "A")

# summer survival 2019
# use the latest replacement date
evAGsSurvD2Dat <- survD2Dat %>%
  filter(focal == 1 & sp == "Ev" & age == "adult") %>%
  group_by(site, plot, treatment, ID) %>%
  summarise(replace_date = max(replace_date),
            survival = 0) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plotDens) %>%
  mutate(survival = replace_na(survival, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0))


#### fit regression ####

# initial fit
evAGsSurvD2FuMod1 <- brm(survival ~ fungicide + (1|site),
                       data = evAGsSurvD2Dat,
                       family = bernoulli,
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b),
                                 prior(cauchy(0, 1), class = sd)),
                       iter = 6000, warmup = 1000, chains = 1)
summary(evAGsSurvD2FuMod1)

# increase chains
evAGsSurvD2FuMod2 <- update(evAGsSurvD2FuMod1, chains = 3,
                            control = list(adapt_delta = 0.99))
summary(evAGsSurvD2FuMod2)
plot(evAGsSurvD2FuMod2)
pp_check(evAGsSurvD2FuMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(0, 1)) %>%
  mutate(survival = fitted(evAGsSurvD2FuMod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(evAGsSurvD2FuMod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(evAGsSurvD2FuMod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]) %>%
    mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))
  
# summarise by site
vizDat <- evAGsSurvD2Dat %>%
  group_by(site, treatment) %>%
  summarise(survival = mean(survival))  

# fit figure
(evAGsSurvD2FuPlot <- ggplot(fitDat, aes(treatment, survival, color = treatment)) +
    geom_errorbar(aes(ymin = survival_lower, ymax = survival_upper), width = 0) +
    geom_point(size = 4, shape = 16) +
    geom_point(data = vizDat, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
    theme_bw())


#### output ####

save(evAGsSurvD2FuMod2, file = "output/elymus_adult_gs_survival_fung_model_2019_density_exp.rda")
save(evAGsSurvD2FuPlot, file = "output/elymus_adult_gs_survival_fung_figure_2019_density_exp.rda")

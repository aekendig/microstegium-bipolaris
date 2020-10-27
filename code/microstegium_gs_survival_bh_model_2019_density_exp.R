##### info ####

# file: microstegium_gs_survival_bh_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate Microstegium growing season survival average for intercept of Beverton-Holt model


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
  expand_grid(ID = as.character(c(1, 2, 3)))

# summer survival 2019
# use the latest replacement date
mvGsSurvD2Dat <- survD2Dat %>%
  filter(focal == 1 & sp == "Mv") %>%
  group_by(site, plot, treatment, ID) %>%
  summarise(replace_date = max(replace_date),
            survival = 0) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plotDens) %>%
  mutate(survival = replace_na(survival, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot))


#### fit regression ####

# initial fit
mvGsSurvD2BhMod1 <- brm(survival ~ fungicide + (1|site/plotr),
                        data = mvGsSurvD2Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
summary(mvGsSurvD2BhMod1)

# increase chains
mvGsSurvD2BhMod2 <- update(mvGsSurvD2BhMod1, chains = 3,
                           control = list(adapt_delta = 0.99))
summary(mvGsSurvD2BhMod2)
plot(mvGsSurvD2BhMod2)
pp_check(mvGsSurvD2BhMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(1, 0)) %>%
  mutate(survival = fitted(mvGsSurvD2BhMod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvGsSurvD2BhMod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvGsSurvD2BhMod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"]) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "water"))

# summarise by site
vizDat <- mvGsSurvD2Dat %>%
  group_by(site, treatment) %>%
  summarise(survival = mean(survival))

# fit figure
(mvGsSurvD2BhPlot <- ggplot(fitDat, aes(treatment, survival, color = treatment)) +
    geom_errorbar(aes(ymin = survival_lower, ymax = survival_upper), width = 0) +
    geom_point(size = 4, shape = 16) +
    geom_point(data = vizDat, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
    theme_bw())


#### output ####

save(mvGsSurvD2BhMod2, file = "output/microstegium_gs_survival_bh_model_2019_density_exp.rda")
save(mvGsSurvD2BhPlot, file = "output/microstegium_gs_survival_bh_figure_2019_density_exp.rda")

##### info ####

# file: microstegium_gs_survival_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/23/20
# goal: estimate Microstegium growing season survival based on density and fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


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
mvGsSurvD2Mod1 <- brm(survival ~ fungicide * (mv_seedling_density + ev_seedling_density + ev_adult_density) + (1|site/plotr),
                        data = mvGsSurvD2Dat,
                        family = bernoulli,
                        prior = c(prior(normal(0, 10), class = Intercept),
                                  prior(normal(0, 10), class = b),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 6000, warmup = 1000, chains = 1)
# 1 divergent transition
summary(mvGsSurvD2Mod1)

# increase chains and max tree depth (got a divergent transitions warning)
mvGsSurvD2Mod2 <- update(mvGsSurvD2Mod1, chains = 3,
                         control = list(adapt_delta = 0.999, max_treedepth = 15))
summary(mvGsSurvD2Mod2)
plot(mvGsSurvD2Mod2)
pp_check(mvGsSurvD2Mod2, nsamples = 50)


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
  mutate(survival = fitted(mvGsSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         survival_lower = fitted(mvGsSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         survival_upper = fitted(mvGsSurvD2Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"])

# modify dataset so that zero plots are repeated
vizDat <- mvGsSurvD2Dat %>%
  select(-c(background, background_sp)) %>%
  full_join(plotsD2)

# fit figure
(mvGsSurvD2Plot <- ggplot(fitDat, aes(background_density, survival, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = survival_lower, ymax = survival_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = vizDat, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = vizDat, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw())
# survival is on average lower with water
# fits are reasonable


#### output ####

save(mvGsSurvD2Mod2, file = "output/microstegium_gs_survival_model_2019_density_exp.rda")
save(mvGsSurvD2Plot, file = "output/microstegium_gs_survival_figure_2019_density_exp.rda")

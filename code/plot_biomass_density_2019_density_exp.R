##### info ####

# file: plot_biomass_density_2019_density_exp
# author: Amy Kendig
# date last edited: 3/11/21
# goal: plot-level biomass and density relationships


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# biomass
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R

# load individual growth models
load("output/mv_biomass_mv_density_model_2019_density_exp.rda")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1))

# plot biomass
# use average of other species in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         sp_age = paste(sp, age, sep = "_")) %>%
  pivot_wider(-c(sp, age),
              names_from = sp_age,
              names_glue = "{sp_age}_biomass",
              values_from = biomass.g) %>%
  mutate(Mv_seedling_biomass = replace_na(Mv_seedling_biomass, 0),
         Ev_seedling_biomass = replace_na(Ev_seedling_biomass, 0),
         Ev_adult_biomass = replace_na(Ev_adult_biomass, 0)) %>%
  select(-none_seedling_biomass) %>%
  left_join(mvBioD2Dat %>%
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment, ) %>%
              summarise(Mv_seedling_biomass_foc = sum(biomass_weight.g)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID %in% c("1", "2", "3")) %>%
              group_by(site, plot, treatment) %>%
              mutate(weight_adj = mean(weight, na.rm = T)) %>%
              ungroup() %>%
              mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                        TRUE ~ weight)) %>%
              group_by(site, plot, treatment) %>%
              summarise(Ev_seedling_biomass_foc = sum(weight)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID == "A") %>%
              select(site, plot, treatment, weight) %>%
              rename(Ev_adult_biomass_foc = weight)) %>%
  mutate(Mv_seedling_biomass = Mv_seedling_biomass + Mv_seedling_biomass_foc,
         Ev_seedling_biomass = Ev_seedling_biomass + Ev_seedling_biomass_foc,
         Ev_adult_biomass = Ev_adult_biomass + Ev_adult_biomass_foc) %>%
  select(-c(Mv_seedling_biomass_foc, Ev_seedling_biomass_foc, Ev_adult_biomass_foc))


# combine
plotD2Dat <- plotDens %>%
  full_join(plotBioD2Dat) %>%
  select(site, plot, treatment, background, density_level, background_density, Mv_seedling_density:Ev_adult_density, Mv_seedling_biomass:Ev_adult_biomass) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0))

# divide
mvPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% 1:4)

evSPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% c(1, 5:7))

evAPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% c(1, 8:10))


#### initial visualizations ####

# figures
ggplot(mvPlotD2Dat, aes(Mv_seedling_density, Mv_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evSPlotD2Dat, aes(Ev_seedling_density, Ev_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evAPlotD2Dat, aes(Ev_adult_density, Ev_adult_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Mv model ####

# priors
mvPlotD2Dat %>%
  filter(plot == 1) %>%
  summarise(b0 = mean(Mv_seedling_biomass/3))

x <- seq(-1, 10, length.out = 100)
y <- dgamma(x, shape = 1, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

x <- seq(0, 10, length.out = 100)
y <- dexp(x, 0.5)
plot(x, y, type = "l")

# model
mvBioDensMod <- brm(data = mvPlotD2Dat, family = gaussian,
                    bf(Mv_seedling_biomass ~ (Mv_seedling_density * b0)/(1 + alpha * Mv_seedling_density), 
                       b0 ~ 0 + treatment, 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(14, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3) 

summary(mvBioDensMod)
pp_check(mvBioDensMod, nsamples = 100)
plot(mvBioDensMod)

# simulate values
mvBioDensSim <- mvPlotD2Dat %>%
  select(treatment) %>%
  unique() %>%
  expand_grid(tibble(Mv_seedling_density = seq(0, 67, length.out = 100))) %>%
  mutate(pred = fitted(mvBioDensMod, newdata = ., re_formula = NA)[, "Estimate"],
         lower = fitted(mvBioDensMod, newdata = ., re_formula = NA)[, "Q2.5"],
         upper = fitted(mvBioDensMod, newdata = ., re_formula = NA)[, "Q97.5"])

# figure
ggplot(mvPlotD2Dat, aes(Mv_seedling_density, Mv_seedling_biomass, color = treatment)) +
  geom_ribbon(data = mvBioDensSim, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = mvBioDensSim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Ev seedling model ####

# priors
evSPlotD2Dat %>%
  filter(plot == 1) %>%
  summarise(b0 = mean(Ev_seedling_biomass/3))

x <- seq(-1, 10, length.out = 100)
y <- dgamma(x, shape = 1, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

x <- seq(0, 10, length.out = 100)
y <- dexp(x, 0.5)
plot(x, y, type = "l")

# model
# evSBioDensMod <- brm(data = evSPlotD2Dat, family = gaussian,
#                     bf(Ev_seedling_biomass ~ (Ev_seedling_density * b0)/(1 + alpha * Ev_seedling_density)^(-sc), 
#                        b0 ~ 0 + treatment, 
#                        alpha ~ 0 + treatment, 
#                        sc ~ 0 + treatment, 
#                        nl = T),
#                     prior <- c(prior(gamma(1, 1), nlpar = "b0", lb = 0),
#                                prior(exponential(0.5), nlpar = "alpha", lb = 0),
#                                prior(exponential(0.5), nlpar = "sc", lb = 0)), # use default for sigma
#                     iter = 6000, warmup = 1000, chains = 3) 
# 407 divergences, fit is basically linear

evSBioDensMod <- brm(data = evSPlotD2Dat, family = gaussian,
                     log(Ev_seedling_biomass) ~ Ev_seedling_density * fungicide, 
                     prior <- c(prior(normal(0, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3) 

summary(evSBioDensMod)
pp_check(evSBioDensMod, nsamples = 100)
plot(evSBioDensMod)

# simulate values
evSBioDensSim <- evSPlotD2Dat %>%
  select(treatment, fungicide) %>%
  unique() %>%
  expand_grid(tibble(Ev_seedling_density = seq(0, 19, length.out = 100))) %>%
  mutate(pred = fitted(evSBioDensMod, newdata = ., re_formula = NA)[, "Estimate"] %>% exp(),
         lower = fitted(evSBioDensMod, newdata = ., re_formula = NA)[, "Q2.5"] %>% exp(),
         upper = fitted(evSBioDensMod, newdata = ., re_formula = NA)[, "Q97.5"] %>% exp())

# figure
ggplot(evSPlotD2Dat, aes(Ev_seedling_density, Ev_seedling_biomass, color = treatment)) +
  geom_ribbon(data = evSBioDensSim, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = evSBioDensSim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Ev adult model ####

# priors
evAPlotD2Dat %>%
  filter(plot == 1) %>%
  summarise(b0 = mean(Ev_adult_biomass/3))

x <- seq(-1, 10, length.out = 100)
y <- dgamma(x, shape = 2, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

x <- seq(0, 10, length.out = 100)
y <- dexp(x, 0.5)
plot(x, y, type = "l")

# model
# evABioDensMod <- brm(data = evAPlotD2Dat, family = gaussian,
#                     bf(Ev_adult_biomass ~ (Ev_adult_density * b0)/(1 + alpha * Ev_adult_density),
#                        b0 ~ 0 + treatment,
#                        alpha ~ 0 + treatment,
#                        nl = T),
#                     prior <- c(prior(gamma(2, 1), nlpar = "b0", lb = 0),
#                                prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma
#                     iter = 6000, warmup = 1000, chains = 3)
# huge error around estimates, almost linear

evABioDensMod <- brm(data = evAPlotD2Dat, family = gaussian,
                     log(Ev_adult_biomass) ~ Ev_adult_density * fungicide, 
                     prior <- c(prior(normal(2, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3) 

summary(evABioDensMod)
pp_check(evABioDensMod, nsamples = 100)
plot(evABioDensMod)

# simulate values
evABioDensSim <- evAPlotD2Dat %>%
  select(treatment, fungicide) %>%
  unique() %>%
  expand_grid(tibble(Ev_adult_density = seq(0, 9, length.out = 100))) %>%
  mutate(pred = fitted(evABioDensMod, newdata = ., re_formula = NA)[, "Estimate"] %>% exp(),
         lower = fitted(evABioDensMod, newdata = ., re_formula = NA)[, "Q2.5"] %>% exp(),
         upper = fitted(evABioDensMod, newdata = ., re_formula = NA)[, "Q97.5"] %>% exp())

# figure
ggplot(evAPlotD2Dat, aes(Ev_adult_density, Ev_adult_biomass, color = treatment)) +
  geom_ribbon(data = evABioDensSim, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = evABioDensSim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### figure ####

# combine datasets
figDat <- mvPlotD2Dat %>%
  rename(biomass = Mv_seedling_biomass,
         density = Mv_seedling_density) %>%
  mutate(plant_group = "Mv seedling") %>%
  full_join(evSPlotD2Dat %>%
              rename(biomass = Ev_seedling_biomass,
                     density = Ev_seedling_density) %>%
              mutate(plant_group = "Ev seedling") ) %>%
  full_join(evAPlotD2Dat %>%
              rename(biomass = Ev_adult_biomass,
                     density = Ev_adult_density) %>%
              mutate(plant_group = "Ev adult") ) %>%
  mutate(plant_group = fct_rev(plant_group),
         treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())

simDat <- mvBioDensSim %>%
  rename(density = Mv_seedling_density) %>%
  mutate(plant_group = "Mv seedling") %>%
  full_join(evSBioDensSim %>%
              rename(density = Ev_seedling_density) %>%
              mutate(plant_group = "Ev seedling") ) %>%
  full_join(evABioDensSim %>%
              rename(density = Ev_adult_density) %>%
              mutate(plant_group = "Ev adult") ) %>%
  mutate(plant_group = fct_rev(plant_group),
         treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())

# labels
labDat <- simDat %>%
  group_by(plant_group) %>%
  summarise(biomass = max(upper)) %>%
  ungroup() %>%
  mutate(treatment = "fungicide")

# figure
pdf("output/plot_biomass_density_2019_density_exp.pdf", width = 10.5, height = 5)
ggplot(figDat, aes(density, biomass, group = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3), aes(color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2,  position = position_dodge(0.3), aes(color = treatment)) +
  geom_ribbon(data = simDat, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = simDat, aes(y = pred, color = treatment)) +
  geom_text(data = labDat, x = 0, hjust = 0, aes(label = plant_group), size = 5) +
  facet_wrap(~plant_group, scales = "free") +
  xlab("Density") +
  ylab(expression(paste("Biomass (g ", m^-2, ")", sep = ""))) +
  scale_color_viridis_d(end = 0.6, name = "Treatment") +
  scale_fill_viridis_d(end = 0.6, name = "Treatment") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.22, 0.15),
        strip.background = element_blank(),
        strip.text = element_blank())
dev.off()


#### values for text ####

posterior_samples(mvBioDensMod) %>%
  transmute(fungicide = (67 * b_b0_treatmentfungicide) / (1 + b_alpha_treatmentfungicide * 67),
            water = (67 * b_b0_treatmentwater) / (1 + b_alpha_treatmentwater * 67)) %>%
  mutate(fung_eff = (fungicide - water) / water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "biomass") %>%
  group_by(treatment) %>%
  mean_hdi(biomass)

posterior_samples(evSBioDensMod) %>%
  rename(b_fungicide_density = "b_Ev_seedling_density:fungicide") %>%
  transmute(water = exp(b_Intercept + b_Ev_seedling_density * 19),
            fungicide = exp(b_Intercept + b_Ev_seedling_density * 19 + b_fungicide + b_fungicide_density * 19)) %>%
  mutate(fung_eff = (fungicide - water) / water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "biomass") %>%
  group_by(treatment) %>%
  mean_hdi(biomass)

posterior_samples(evABioDensMod) %>%
  rename(b_fungicide_density = "b_Ev_adult_density:fungicide") %>%
  transmute(water = exp(b_Intercept + b_Ev_adult_density * 9),
            fungicide = exp(b_Intercept + b_Ev_adult_density * 9 + b_fungicide + b_fungicide_density * 9)) %>%
  mutate(fung_eff = (fungicide - water) / water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "biomass") %>%
  group_by(treatment) %>%
  mean_hdi(biomass)


#### output ####
# plot-level datasheet
write_csv(plotD2Dat, "intermediate-data/plot_biomass_density_2019_density_exp.csv")

# models
save(mvBioDensMod, file = "output/mv_plot_biomass_density_model_2019_density_exp.rda")
save(evSBioDensMod, file = "output/ev_seedling_plot_biomass_density_model_2019_density_exp.rda")
save(evABioDensMod, file = "output/ev_adult_plot_biomass_density_model_2019_density_exp.rda")

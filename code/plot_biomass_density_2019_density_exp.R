##### info ####

# file: plot_biomass_density_2019_density_exp
# author: Amy Kendig
# date last edited: 3/16/21
# goal: plot-level biomass and density relationships


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)

# import data
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R
envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") 
# covariate_data_processing_2018_density_exp

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
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = as.factor(plot)) %>%
  left_join(envD1Dat %>%
              select(site, plot, treatment, soil_moisture_oct.prop, canopy_cover.prop))

# divide
mvPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% 1:4) %>%
  mutate(soil_moisture_oct.c = scale(soil_moisture_oct.prop, center = T, scale = F)[,1],
         canopy_cover.c = scale(canopy_cover.prop, center = T, scale = F)[,1],
         log_tot_bio = log(Mv_seedling_biomass))

evSPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% c(1, 5:7)) %>%
  mutate(soil_moisture_oct.c = scale(soil_moisture_oct.prop, center = T, scale = F)[,1],
         canopy_cover.c = scale(canopy_cover.prop, center = T, scale = F)[,1],
         log_tot_bio = log(Ev_seedling_biomass))

evAPlotD2Dat <- plotD2Dat %>%
  filter(plot %in% c(1, 8:10)) %>%
  mutate(soil_moisture_oct.c = scale(soil_moisture_oct.prop, center = T, scale = F)[,1],
         canopy_cover.c = scale(canopy_cover.prop, center = T, scale = F)[,1],
         log_tot_bio = log(Ev_adult_biomass))


#### initial visualizations ####

# biomass by density
ggplot(mvPlotD2Dat, aes(Mv_seedling_density, Mv_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evSPlotD2Dat, aes(Ev_seedling_density, Ev_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evAPlotD2Dat, aes(Ev_adult_density, Ev_adult_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# biomass by soil moisture
ggplot(mvPlotD2Dat, aes(soil_moisture_oct.prop, Mv_seedling_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

ggplot(evSPlotD2Dat, aes(soil_moisture_oct.prop, Ev_seedling_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

ggplot(evAPlotD2Dat, aes(soil_moisture_oct.prop, Ev_adult_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

# biomass by canopy cover
ggplot(mvPlotD2Dat, aes(canopy_cover.prop, Mv_seedling_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

ggplot(evSPlotD2Dat, aes(canopy_cover.prop, Ev_seedling_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

ggplot(evAPlotD2Dat, aes(canopy_cover.prop, Ev_adult_biomass, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = F)

# soil and canopy
ggplot(plotD2Dat, aes(canopy_cover.prop, soil_moisture_oct.prop)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x)

cor.test(~ canopy_cover.prop + soil_moisture_oct.prop, data = plotD2Dat)
# not significantly correlated

# soil moisture by density (because it was measured at the end of the Y1 growing season)
ggplot(mvPlotD2Dat, aes(Mv_seedling_density, soil_moisture_oct.prop, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evSPlotD2Dat, aes(Ev_seedling_density, soil_moisture_oct.prop, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(evAPlotD2Dat, aes(Ev_adult_density, soil_moisture_oct.prop, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Mv models ####

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
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(14, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 

summary(mvBioDensMod)
pp_check(mvBioDensMod, nsamples = 100)
plot(mvBioDensMod)

# simulate values
mvBioDensSim <- mvPlotD2Dat %>%
  select(treatment) %>%
  unique() %>%
  mutate(site = "D5") %>%
  expand_grid(tibble(Mv_seedling_density = seq(0, 67, length.out = 100))) %>%
  mutate(pred = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

# figure
ggplot(mvPlotD2Dat, aes(Mv_seedling_density, Mv_seedling_biomass, color = treatment)) +
  geom_ribbon(data = mvBioDensSim, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = mvBioDensSim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# evaluate environmental covariates
mvPlotD2Dat %>%
  filter(plot == 1 & fungicide == 0) %>%
  summarise(bio = mean(log_tot_bio))

mvBioDensMod2 <- brm(data = mvPlotD2Dat, family = gaussian,
                     log_tot_bio ~ plotf * fungicide + soil_moisture_oct.c + canopy_cover.c + (1|site), 
                     prior <- c(prior(normal(3, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 

summary(mvBioDensMod2)
#  no sig effect of soil moisture or canopy cover


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

evSPlotD2Dat %>%
  filter(plot == 1 & fungicide == 0) %>%
  summarise(bio = mean(log_tot_bio))

evSBioDensMod <- brm(data = evSPlotD2Dat, family = gaussian,
                     log_tot_bio ~ Ev_seedling_density * fungicide + soil_moisture_oct.c + canopy_cover.c + (1|site), 
                     prior <- c(prior(normal(0, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 

summary(evSBioDensMod)
# no significant effect of soil or canopy

evSBioDensMod2 <- update(evSBioDensMod, formula. = ~ . - soil_moisture_oct.c - canopy_cover.c)
summary(evSBioDensMod2)
pp_check(evSBioDensMod2, nsamples = 100)
plot(evSBioDensMod2)

# simulate values
evSBioDensSim <- evSPlotD2Dat %>%
  select(treatment, fungicide) %>%
  unique() %>%
  mutate(site = "D5") %>%
  expand_grid(tibble(Ev_seedling_density = seq(0, 19, length.out = 100))) %>%
  mutate(pred = fitted(evSBioDensMod2, newdata = ., allow_new_levels = T)[, "Estimate"] %>% exp(),
         lower = fitted(evSBioDensMod2, newdata = ., allow_new_levels = T)[, "Q2.5"] %>% exp(),
         upper = fitted(evSBioDensMod2, newdata = ., allow_new_levels = T)[, "Q97.5"] %>% exp())

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
                     log(Ev_adult_biomass) ~ Ev_adult_density * fungicide + soil_moisture_oct.c + canopy_cover.c + (1|site), 
                     prior <- c(prior(normal(2, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.9999)) 

summary(evABioDensMod)
# no significant effect of soil or canopy

evABioDensMod2 <- update(evABioDensMod, formula. = ~ . - soil_moisture_oct.c - canopy_cover.c)
summary(evABioDensMod2)
pp_check(evABioDensMod2, nsamples = 100)
plot(evABioDensMod2)

# simulate values
evABioDensSim <- evAPlotD2Dat %>%
  select(treatment, fungicide) %>%
  unique() %>%
  mutate(site = "D5") %>%
  expand_grid(tibble(Ev_adult_density = seq(0, 9, length.out = 100))) %>%
  mutate(pred = fitted(evABioDensMod2, newdata = ., allow_new_levels = T)[, "Estimate"] %>% exp(),
         lower = fitted(evABioDensMod2, newdata = ., allow_new_levels = T)[, "Q2.5"] %>% exp(),
         upper = fitted(evABioDensMod2, newdata = ., allow_new_levels = T)[, "Q97.5"] %>% exp())

# figure
ggplot(evAPlotD2Dat, aes(Ev_adult_density, Ev_adult_biomass, color = treatment)) +
  geom_ribbon(data = evABioDensSim, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
  geom_line(data = evABioDensSim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### combined figure ####

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
pdf("output/plot_biomass_density_presentation_2019_density_exp.pdf", width = 10.5, height = 5)
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


#### separate figure ####

# change treatment levels
mvPlotD2Dat2 <- mvPlotD2Dat %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())
mvBioDensSim2 <- mvBioDensSim %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())

evSPlotD2Dat2 <- evSPlotD2Dat %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())
evSBioDensSim2 <- evSBioDensSim %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())

evAPlotD2Dat2 <- evAPlotD2Dat %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())
evABioDensSim2 <- evABioDensSim %>%
  mutate(treatment = case_when(treatment == "water" ~ "water (control)",
                               TRUE ~ "fungicide") %>%
           fct_rev())

# figure settings #
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

mv_fig <- ggplot(mvPlotD2Dat2, aes(Mv_seedling_density, Mv_seedling_biomass, group = treatment)) +
  geom_ribbon(data = mvBioDensSim2, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, color = NA) +
  geom_line(data = mvBioDensSim2, aes(y = pred, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, position = position_dodge(1), aes(color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2,  position = position_dodge(1), aes(color = treatment)) +
  xlab(expression(paste(italic("M. vimineum"), " density (", m^-2, ")", sep = ""))) +
  ylab(expression(paste(italic("M. vimineum"), " biomass (g ", m^-2, ")", sep = ""))) +
  scale_color_viridis_d(end = 0.6, name = "Treatment") +
  scale_fill_viridis_d(end = 0.6, name = "Treatment") +
  fig_theme

evS_fig <- ggplot(evSPlotD2Dat2, aes(Ev_seedling_density, Ev_seedling_biomass, group = treatment)) +
  geom_ribbon(data = evSBioDensSim2, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, color = NA) +
  geom_line(data = evSBioDensSim2, aes(y = pred, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, position = position_dodge(0.3), aes(color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2,  position = position_dodge(0.3), aes(color = treatment)) +
  xlab(expression(paste(italic("E. virginicus"), " seedling density (", m^-2, ")", sep = ""))) +
  ylab(expression(paste(italic("E. virginicus"), " seedling biomass (g ", m^-2, ")", sep = ""))) +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  fig_theme +
  theme(legend.position = c(0.35, 0.88),
        axis.title.y = element_text(size = 10, hjust = -1),
        axis.title.x = element_text(size = 10, hjust = 0.85))

evA_fig <- ggplot(evAPlotD2Dat2, aes(Ev_adult_density, Ev_adult_biomass, group = treatment)) +
  geom_ribbon(data = evABioDensSim2, aes(y = pred, ymin = lower, ymax = upper, fill = treatment), alpha = 0.3, color = NA) +
  geom_line(data = evABioDensSim2, aes(y = pred, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, position = position_dodge(0.1), aes(color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2,  position = position_dodge(0.1), aes(color = treatment)) +
  xlab(expression(paste(italic("E. virginicus"), " adult density (", m^-2, ")", sep = ""))) +
  ylab(expression(paste(italic("E. virginicus"), " adult biomass (g ", m^-2, ")", sep = ""))) +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  fig_theme +
  theme(axis.title.x = element_text(size = 10, hjust = 0.75))

pdf("output/plot_biomass_density_2019_density_exp.pdf", width = 7, height = 3)
plot_grid(mv_fig, evS_fig, evA_fig,
          nrow = 1,
          labels = LETTERS[1:3])
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

posterior_samples(evSBioDensMod2) %>%
  rename(b_fungicide_density = "b_Ev_seedling_density:fungicide") %>%
  transmute(water = exp(b_Intercept + b_Ev_seedling_density * 19),
            fungicide = exp(b_Intercept + b_Ev_seedling_density * 19 + b_fungicide + b_fungicide_density * 19)) %>%
  mutate(fung_eff = (fungicide - water) / water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "biomass") %>%
  group_by(treatment) %>%
  mean_hdi(biomass)

posterior_samples(evABioDensMod2) %>%
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
save(evSBioDensMod2, file = "output/ev_seedling_plot_biomass_density_model_2019_density_exp.rda")
save(evABioDensMod, file = "output/ev_adult_plot_biomass_density_model_2019_density_exp.rda")

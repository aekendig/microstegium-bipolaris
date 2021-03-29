##### info ####

# file: ev_fungicide_figure_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/28/21
# goal: figure of statistically significant fungicide effects


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)

# load models
load("output/focal_growth_no_background_model_2019_density_exp.rda")
load("output/evS_tillers_mv_density_model_2018_density_exp.rda")
load("output/ev_germination_model_2018_2019_density_exp.rda")

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())


#### model outputs ####

# growth
summary(growthD2Mod)

growthDat <- posterior_samples(growthD2Mod) %>%
  rename(b_fungicide_seedling = "b_fungicide:plant_groupEv_seedling") %>%
  transmute(water = exp(b_Intercept + b_plant_groupEv_seedling),
            fung = exp(b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_seedling),
            perc = (fung - water)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "pred") %>%
  group_by(treatment) %>% mean_hdi(pred) %>%
  mutate(treatment = recode(treatment, 
                            "fung" = "Fungicide",
                            "water" = "Water (control)") %>%
           fct_rev())

growthDat2 <- growthDat %>%
  filter(treatment != "perc")

# Mv competition
summary(evSMvD1Mod)

compDat <- posterior_samples(evSMvD1Mod) %>%
  transmute(water = b_alpha_treatmentwater,
            fung = b_alpha_treatmentfungicide,
            perc = (fung - water)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "pred") %>%
  group_by(treatment) %>% mean_hdi(pred) %>%
  mutate(treatment = recode(treatment, 
                            "fung" = "Fungicide",
                            "water" = "Water (control)") %>%
           fct_rev())

compDat2 <- compDat %>%
  filter(treatment != "perc")

# germination
summary(evGermMod)

germDat <- posterior_samples(evGermMod) %>%
  transmute(water = exp(b_Intercept) / (1 + exp(b_Intercept)),
            fung = exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide)),
            perc = (fung - water)/water) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "pred") %>%
  group_by(treatment) %>% mean_hdi(pred) %>%
  mutate(treatment = recode(treatment, 
                            "fung" = "Fungicide",
                            "water" = "Water (control)") %>%
           fct_rev())

germDat2 <- germDat %>%
  filter(treatment != "perc")


#### treatment figure ####

growthFig <- ggplot(growthDat2, aes(treatment, pred, color = treatment)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0) +
  geom_point(size = 2) +
  ylab(expression(paste("Biomass (g ", plant^-1, ")", sep = ""))) +
  annotate(geom = "text", x = 0.5, y = max(growthDat2$.upper), label = "Fungicide: +26%", hjust = 0, vjust = 1, size = 2.5) +
  scale_color_viridis_d(end = 0.5) +
  fig_theme +
  theme(axis.text.x = element_blank())

compFig <- ggplot(compDat2, aes(treatment, pred, color = treatment)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0) +
  geom_point(size = 2) +
  ylab(expression(paste(italic("M. vimineum "), alpha, sep = ""))) +
  annotate(geom = "text", x = 1.4, y = max(compDat2$.upper), label = "Fungicide: -89%", hjust = 0, size = 2.5, vjust = 1) +
  scale_color_viridis_d(end = 0.5) +
  fig_theme +
  theme(axis.text.x = element_blank())

germFig <- ggplot(germDat2, aes(treatment, pred, color = treatment)) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0) +
  geom_point(size = 2) +
  ylab("Germination") +
  annotate(geom = "text", x = 1.4, y = max(germDat2$.upper), label = "Fungicide: -30%", hjust = 0, size = 2.5, vjust = 1) +
  scale_color_viridis_d(end = 0.5) +
  fig_theme

# output 

pdf("output/ev_fungicide_figure_2018_2019_density_exp.pdf", width = 2, height = 5.5)
plot_grid(growthFig, compFig, germFig,
          ncol = 1,
          labels = LETTERS[1:3],
          rel_heights = c(1, 1, 1.1))
dev.off()


#### fungicide effect figures ####

# combine datasets
combDat <- growthDat %>%
  mutate(trait = "Biomass") %>%
  full_join(compDat %>%
              mutate(trait = "Competition")) %>%
  full_join(germDat %>%
              mutate(trait = "Germination")) %>%
  filter(treatment == "perc") %>%
  mutate(trait = fct_relevel(trait, "Biomass", "Germination"))

pdf("output/ev_fungicide_effect_figure_2018_2019_density_exp.pdf", width = 2.75, height = 2.75)
ggplot(combDat, aes(trait, pred)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0) +
  geom_point(size = 2) +
  ylab(expression(paste("Fungicide effect on ", italic("E. virginicus")))) +
  fig_theme
dev.off()
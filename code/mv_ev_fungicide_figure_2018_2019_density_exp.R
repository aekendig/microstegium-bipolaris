##### info ####

# file: mv_ev_fungicide_figure_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/12/21
# goal: figure of all fungicide effects


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)

# load models
load("output/mv_germination_model_2018_density_exp.rda")
load("output/ev_germination_model_2018_2019_density_exp.rda")
load("output/focal_growth_no_background_model_2018_density_exp.rda")
load("output/focal_growth_no_background_model_2019_density_exp.rda")
load("output/mv_growing_season_survival_model_2018_density_exp.rda")
load("output/evS_growing_season_survival_model_2018_density_exp.rda")
load("output/evA_growing_season_survival_model_2018_density_exp.rda")
load("output/mv_growing_season_survival_model_2019_density_exp.rda")
load("output/evS_growing_season_survival_model_2019_density_exp.rda")
load("output/evA_growing_season_survival_model_2019_density_exp.rda")
load("output/evS_winter_survival_model_2018_density_exp.rda")
load("output/evA_winter_survival_model_2018_density_exp.rda")

# load data
mvSeedsBioD2Dat <- read_csv("intermediate-data/mv_seeds_per_biomass_2019_density_exp.csv")
evSSeedsBioD2Dat <- read_csv("intermediate-data/evS_seeds_per_biomass_2019_density_exp.csv")
evASeedsBioD2Dat <- read_csv("intermediate-data/evA_seeds_per_biomass_2019_density_exp.csv")
mvSeedsBioD2Sim <- read_csv("intermediate-data/mv_seeds_per_biomass_sim_2019_density_exp.csv")
evSSeedsBioD2Sim <- read_csv("intermediate-data/evS_seeds_per_biomass_sim_2019_density_exp.csv")
evASeedsBioD2Sim <- read_csv("intermediate-data/evA_seeds_per_biomass_sim_2019_density_exp.csv")

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

col_pal = c("#440154FF", "#472F7DFF", "#39568CFF", "#35B779FF")
shape_pal = c(21, 23, 22, 24)


#### sample processing function ####

# for binomial models
samp_bin_fun <- function(mod){
  
  out <- posterior_samples(mod) %>%
    mutate(water = exp(b_Intercept) / (1 + exp(b_Intercept)),
           fung = exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide))) %>%
    transmute(prop_change = (fung - water) / water) %>%
    median_hdci(prop_change)

  return(out)
  
}


#### posterior samples ####

# 2018 Mv germination 
mvGermD1Samps <- samp_bin_fun(mvGermD1Mod) %>%
  mutate(plant_group = "Mv",
         response = "germination",
         year = 2018)

# 2018-2019 Ev germination
evGermSamps <- posterior_samples(evGermMod) %>%
  rename(b_fungicide_19 = "b_fungicide:yearf2019") %>%
  mutate(water_18 = exp(b_Intercept) / (1 + exp(b_Intercept)),
         fung_18 = exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide)),
         water_19 = exp(b_Intercept + b_yearf2019) / (1 + exp(b_Intercept + b_yearf2019)),
         fung_19 = exp(b_Intercept + b_fungicide + b_yearf2019 + b_fungicide_19) / (1 + exp(b_Intercept + b_fungicide + b_yearf2019 + b_fungicide_19))) %>%
  transmute(prop_18 = (fung_18 - water_18) / water_18,
            prop_19 = (fung_19 - water_19) / water_19) %>%
  pivot_longer(cols = everything(),
               names_to = "year",
               values_to = "prop_change") %>%
  group_by(year) %>% median_hdci(prop_change) %>%
  mutate(year = recode(year, 
                       "prop_18" = "2018",
                       "prop_19" = "2019") %>% 
           as.numeric(),
         plant_group = "Ev",
         response = "germination")

# 2018 growth
growthD1Samps <- posterior_samples(growthD1Mod) %>%
  rename(b_fungicide_EvS = "b_fungicide:plant_groupEv_seedling",
         b_fungicide_EvA = "b_fungicide:plant_groupEv_adult") %>%
  mutate(water_Mv = exp(b_Intercept),
            fung_Mv = exp(b_Intercept + b_fungicide),
            water_EvS = exp(b_Intercept + b_plant_groupEv_seedling),
            fung_EvS = exp(b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_EvS),
            water_EvA = exp(b_Intercept + b_plant_groupEv_adult),
            fung_EvA = exp(b_Intercept + b_plant_groupEv_adult + b_fungicide + b_fungicide_EvA)) %>%
  transmute(prop_Mv = (fung_Mv - water_Mv) / water_Mv,
            prop_EvS = (fung_EvS - water_EvS) / water_EvS,
            prop_EvA = (fung_EvA - water_EvA) / water_EvA) %>%
  pivot_longer(cols = everything(),
               names_to = "plant_group",
               values_to = "prop_change") %>%
  group_by(plant_group) %>% median_hdci(prop_change) %>%
  mutate(plant_group = recode(plant_group, 
                              "prop_EvA" = "Ev adult",
                              "prop_EvS" = "Ev seedling",
                              "prop_Mv" = "Mv"),
         year = 2018,
         response = "growth")

# 2019 growth
growthD2Samps <- posterior_samples(growthD2Mod) %>%
  rename(b_fungicide_EvS = "b_fungicide:plant_groupEv_seedling",
         b_fungicide_EvA = "b_fungicide:plant_groupEv_adult") %>%
  mutate(water_Mv = exp(b_Intercept),
         fung_Mv = exp(b_Intercept + b_fungicide),
         water_EvS = exp(b_Intercept + b_plant_groupEv_seedling),
         fung_EvS = exp(b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_EvS),
         water_EvA = exp(b_Intercept + b_plant_groupEv_adult),
         fung_EvA = exp(b_Intercept + b_plant_groupEv_adult + b_fungicide + b_fungicide_EvA)) %>%
  transmute(prop_Mv = (fung_Mv - water_Mv) / water_Mv,
            prop_EvS = (fung_EvS - water_EvS) / water_EvS,
            prop_EvA = (fung_EvA - water_EvA) / water_EvA) %>%
  pivot_longer(cols = everything(),
               names_to = "plant_group",
               values_to = "prop_change") %>%
  group_by(plant_group) %>% median_hdci(prop_change) %>%
  mutate(plant_group = recode(plant_group, 
                              "prop_EvA" = "Ev adult",
                              "prop_EvS" = "Ev seedling",
                              "prop_Mv" = "Mv"),
         year = 2019,
         response = "growth")

# 2018 Mv survival
mvSurvD1Samps <- samp_bin_fun(mvSurvD1Mod) %>%
  mutate(plant_group = "Mv",
         response = "survival",
         year = 2018)

# 2018 EvS survival
evSSurvD1Samps <- samp_bin_fun(evSSurvD1Mod) %>%
  mutate(plant_group = "Ev seedling",
         response = "survival",
         year = 2018)

# 2018 EvA survival
evASurvD1Samps <- samp_bin_fun(evASurvD1Mod) %>%
  mutate(plant_group = "Ev adult",
         response = "survival",
         year = 2018)

# 2019 Mv survival
mvSurvD2Samps <- samp_bin_fun(mvSurvD2Mod) %>%
  mutate(plant_group = "Mv",
         response = "survival",
         year = 2019)

# 2019 EvS survival
evSSurvD2Samps <- samp_bin_fun(evSSurvD2Mod) %>%
  mutate(plant_group = "Ev seedling",
         response = "survival",
         year = 2019)

# 2019 EvA survival
evASurvD2Samps <- samp_bin_fun(evASurvD2Mod) %>%
  mutate(plant_group = "Ev adult",
         response = "survival",
         year = 2019)

# winter EvS survival
evSWinSurvD1Samps <- samp_bin_fun(evSWinSurvD1Mod) %>%
  mutate(plant_group = "Ev seedling",
         response = "winter survival",
         year = 2018)

# 2019 EvA survival
evAWinSurvD1Samps <- samp_bin_fun(evAWinSurvD1Mod) %>%
  mutate(plant_group = "Ev adult",
         response = "winter survival",
         year = 2018)


#### combine data ####

dat <- mvGermD1Samps %>%
  full_join(evGermSamps) %>%
  full_join(growthD1Samps) %>%
  full_join(growthD2Samps) %>%
  full_join(mvSurvD1Samps) %>%
  full_join(evSSurvD1Samps) %>%
  full_join(evASurvD1Samps) %>%
  full_join(mvSurvD2Samps) %>%
  full_join(evSSurvD2Samps) %>%
  full_join(evASurvD2Samps) %>%
  full_join(evSWinSurvD1Samps) %>%
  full_join(evAWinSurvD1Samps) %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))


#### fungicide effect figures ####

# 2018
feD1Fig <- dat %>%
  filter(year == 2018) %>%
  ggplot(aes(response, prop_change, color = plant_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = plant_group), width = 0, position = position_dodge(0.5), alpha = 0.7) +
  geom_point(position = position_dodge(0.5), aes(shape = plant_group, fill = sig)) +
  geom_text(x = 1, y = 0.92, label = "2018", color = "black", size = 2.5, check_overlap = T, hjust = 1) +
  scale_shape_manual(values = shape_pal, name = "Plant group") +
  scale_color_manual(values = col_pal, name = "Plant group") +
  scale_fill_manual(values = c("white", "black"), guide = F) +
  ylab("Proportional change\ndue to fungicide") +
  xlab("Plant response") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1))

# 2019
feD2Fig <- dat %>%
  filter(year == 2019) %>%
  ggplot(aes(response, prop_change, color = plant_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = plant_group), width = 0, position = position_dodge(0.5), alpha = 0.7) +
  geom_point(position = position_dodge(0.5), aes(shape = plant_group, fill = sig)) +
  geom_text(x = 1, y = 0.54, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 1) +
  scale_shape_manual(values = shape_pal, name = "Plant group") +
  scale_color_manual(values = col_pal, name = "Plant group") +
  scale_fill_manual(values = c("white", "black"), guide = F) +
  ylab("Proportional change\ndue to fungicide") +
  xlab("Plant response") +
  fig_theme +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1)) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))


#### biomass seeds figure ####

# combine data
bioSeedsDat <- mvSeedsBioD2Dat %>%
  mutate(ID = as.character(ID)) %>%
  full_join(evSSeedsBioD2Dat %>%
              mutate(ID = as.character(ID))) %>%
  full_join(evASeedsBioD2Dat) %>%
  filter(treatment == "water") %>% # weak/no fungicide effects
  mutate(plant_group = str_replace(plant_group, "_", " "),
         plant_group = recode(plant_group, "Mv seedling" = "Mv"))

bioSeedsSim <- mvSeedsBioD2Sim %>%
  mutate(plant_group = "Mv") %>%
  full_join(evSSeedsBioD2Sim %>%
              mutate(plant_group = "Ev seedling")) %>%
  full_join(evASeedsBioD2Sim %>%
              mutate(plant_group = "Ev adult")) %>%
  filter(fungicide == "water")

seedsBioFig <- ggplot(bioSeedsSim, aes(x = log_bio, y = pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = plant_group), alpha = 0.5) +
  geom_line(aes(color = plant_group)) +
  geom_point(data = bioSeedsDat, alpha = 0.7, size = 0.7, aes(y = log_seeds, shape = plant_group, color = plant_group)) +
  geom_text(x = -1.35, y = 9.25, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
  scale_shape_manual(values = shape_pal[2:4], guide = F) +
  scale_color_manual(values = col_pal[2:4], guide = F) +
  scale_fill_manual(values = col_pal[2:4], guide = F) +
  ylab("Seeds per plant") +
  xlab("Plant biomass (g)") +
  fig_theme +
  theme(legend.position = c(0.8, 0.2),
        axis.title.x = element_text(size = 10))

#### full figure ####

# extract legend
leg = get_legend(feD2Fig)

# legend and third figure
figC <- plot_grid(seedsBioFig, leg, 
                  nrow = 2,
                  rel_heights = c(1, 0.4))

pdf("output/mv_ev_fungicide_effect_figure_2018_2019_density_exp.pdf", width = 7, height = 3)
plot_grid(feD1Fig, feD2Fig + theme(legend.position = "none"), figC,
          labels = c("A", "B", "C"),
          rel_widths = c(1, 0.9, 1),
          nrow = 1)
dev.off()
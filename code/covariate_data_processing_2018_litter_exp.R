##### info ####

# file: covariate_data_processing_2018_litter_exp
# author: Amy Kendig
# date last edited: 3/6/20
# goal: create a dataset of covariates


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)
library(cowplot)

# import data
soil <- read_csv("./data/soil_moisture_jun_2018_litter_exp.csv")
trees <- read_csv("./data/canopy_cover_trees_jun_2018_density_litter_exp.csv")
canopy <- read_csv("./data/canopy_cover_oct_2018_litter_exp.csv")
plots <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


#### edit data ####

# use proportion
soil2 <- soil %>%
  mutate(soil_moisture.prop = soil_moisture.vwc / 100) %>%
  select(site, plot, soil_moisture.prop)

# look at canopy notes
unique(canopy$processing_notes)

#  proportion canopy cover
canopy2 <- canopy %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                       type == "c" ~ (count * 1.04) / 100)) %>%
  select(site, plot, canopy_cover.prop)
  
# site level proportion canopy cover, 
canopy_site <- trees %>%
  filter(experiment == "litter" & !is.na(count)) %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                       type == "c" ~ (count * 1.04) / 100))  %>%
  group_by(site) %>%
  summarise(canopy_cover.prop = mean(canopy_cover.prop)) %>%
  ungroup()

# site level basal tree stand area
stand_site <- trees %>%
  filter(experiment == "litter" & !is.na(trees)) %>%
  mutate(stand_area.m2ha = trees * 10)  %>%
  group_by(site) %>%
  summarise(stand_area.m2ha = mean(stand_area.m2ha)) %>%
  ungroup()

# merge plot data
plot_dat <- full_join(plots, soil2) %>%
  full_join(canopy2)

# merge site data
site_dat <- full_join(canopy_site, stand_site)


#### correlations ####

# plot data
# plot_dat %>%
#   select(litter, litter_weight.g, soil_moisture.prop, canopy_cover.prop) %>%
#   ggpairs()
# # weak relationships among the variables
# 
# # site data
# cor.test(site_dat$canopy_cover.prop, site_dat$stand_area.m2ha)
# uncorrelated


#### site figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# base theme
base_theme <- theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none")

# soil moisture
fig_sp <- ggplot(plot_dat, aes(x = site, y = soil_moisture.prop)) +
  stat_summary(geom = "point", fun.y = mean, size = 3) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, width = 0.1) +
  base_theme +
  xlab("Site") +
  ylab("Soil moisture")

# plot canopy cover
fig_cp <- ggplot(plot_dat, aes(x = site, y = canopy_cover.prop)) +
  stat_summary(geom = "point", fun.y = mean, size = 3) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, width = 0.1) +
  base_theme +
  xlab("Site") +
  ylab("Canopy cover")

# basal stand area
fig_bs <- ggplot(site_dat, aes(x = site, y = stand_area.m2ha)) +
  geom_point(size = 3) +
  base_theme +
  xlab("Site") +
  ylab(expression(paste("Basal stand area (", m^2, ha^{-1}, ")", sep = "")))

# site canopy cover
fig_cs <- ggplot(site_dat, aes(x = site, y = canopy_cover.prop)) +
  geom_point(size = 3) +
  base_theme +
  xlab("Site") +
  ylab("Canopy cover")

# combine and save
pdf("./output/site_covariates_2018_litter_exp.pdf", width = 4, height = 4)
plot_grid(fig_sp, fig_cp, fig_bs, fig_cs,
         nrow = 2,
         labels = letters[1:4],
         label_size = lg_txt)
dev.off()


#### save data ####
write_csv(plot_dat, "./intermediate-data/plot_covariates_2018_litter_exp.csv")
write_csv(site_dat, "./intermediate-data/site_covariates_2018_litter_exp.csv")

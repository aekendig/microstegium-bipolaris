#### leaf scan data needs to be checked - see leaf_scan_data_processing scripts ####


#### info ####

# file: plot_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 4/13/21
# goal: combine plot-scale data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
fSevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
evBgSevD1Dat <- read_csv("intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
mvBgSevD1Dat <- read_csv("intermediate-data/mv_background_leaf_scans_2018_density_exp.csv")
# all above: leaf_scans_data_processing_2018_density_exp.R
bgDisD1Dat <- read_csv("data/bg_disease_late_aug_2018_density_exp.csv")
mvBioD1Dat <- read_csv("intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")
# mv_biomass_data_processing_2018_density_exp.R
mvSeedsD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp


#### Mv severity ####

# tiller counts
mvDisD1Dat <- bgDisD1Dat %>%
  filter(plot %in% c(2, 3, 4)) %>%
  select(site, plot, treatment, tiller, leaves_tot, leaves_infec) %>%
  mutate(focal = 0) %>%
  full_join(fSevD1Dat %>%
              filter(month == "late_aug" & plot %in% c(2, 3, 4) & sp == "Mv") %>%
              select(site, plot, treatment, ID, leaves_tot, leaves_infec, focal)) %>%
  filter(!is.na(leaves_tot)) %>%
  group_by(site, treatment) %>%
  summarise(leaves_tot = mean(leaves_tot),
            leaves_infec = mean(leaves_infec)) %>%
  ungroup()

ls_mv_late_aug2 %>%
  filter(month == "late_aug" & plot == 2 & treatment == "water" & site == "D3")

# severity
mvSevD1Dat <- mvBgSevD1Dat %>%
  filter(month == "late_aug" & plot %in% c(2, 3, 4)) %>%
  full_join(fSevD1Dat %>%
              filter(month == "late_aug" & plot %in% c(2, 3, 4) & sp == "Mv")) %>%
  group_by(site, treatment) %>%
  summarise(leaf_area.pix = sum(leaf_area.pix) / sum(leaf_count), # manually checked leaf counts
            lesion_area.pix = sum(lesion_area.pix) / sum(leaf_count)) %>%
  ungroup() %>%
  full_join(mvDisD1Dat) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(site, treatment, severity) %>%
  pivot_wider(names_from = treatment,
              values_from = severity)

# t-test
t.test(mvSevD1Dat$fungicide, mvSevD1Dat$water, paired = T)
# not significantly different


#### Mv biomass ####

# g per m2
# mv plots
# make wide
mvBioD1Dat2 <- mvBioD1Dat %>%
  mutate(biomass.g_m2 = bio.g * 0.25 * 0.49) %>%
  filter(plot %in% c(2, 3, 4))

mvBioD1DatW <- mvBioD1Dat2 %>%
  select(site, plot, treatment, biomass.g_m2) %>%
  pivot_wider(names_from = treatment,
              values_from = biomass.g_m2)

# figure
ggplot(mvBioD1Dat2, aes(x = treatment, y = biomass.g_m2)) +
  geom_boxplot()

# t-test
t.test(mvBioD1DatW$fungicide, mvBioD1DatW$water, paired = T)


#### Mv seeds ####

# seeds per m2
# mv plots
# make wide
mvSeedsD1Dat2 <- mvSeedsD1Dat %>%
  mutate(seeds_soil = seeds_soil * (0.05 * 0.25),
         seeds_bio = seeds_bio * (0.49 * 0.25),
         seeds = seeds_soil + seeds_bio) %>%
  filter(plot %in% c(2, 3, 4))

mvSeedsD1DatW <- mvSeedsD1Dat2 %>%
  select(site, plot, treatment, seeds) %>%
  pivot_wider(names_from = treatment,
              values_from = seeds)

# figure
ggplot(mvSeedsD1Dat2, aes(x = treatment, y = seeds)) +
  geom_boxplot()

# t-test
t.test(mvSeedsD1DatW$fungicide, mvSeedsD1DatW$water, paired = T)
# samples not representative of density treatments, from high density areas

# high density plots
mvSeedsD1DatW2 <- mvSeedsD1DatW %>%
  filter(plot == 4)

t.test(mvSeedsD1DatW2$fungicide, mvSeedsD1DatW2$water, paired = T)
# similar difference, not sig
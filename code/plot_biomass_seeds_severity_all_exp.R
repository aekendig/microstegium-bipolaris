#### info ####

# file: plot_biomass_seeds_severity_all_exp
# author: Amy Kendig
# date last edited: 4/12/21
# goal: plot-level biomass/seeds vs. severity for multiple experiments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
# ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
flory <- read_csv("intermediate-data/flory_2011_extracted_figure.csv")
stricker <- read_csv("intermediate-data/stricker_2016_extracted_figure.csv")
# two above: extracting_data_from_figures.R


#### edit data ####

# Ev seeds
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0), # only have seed number, not spikelet number
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  filter((plot == 7 & age == "seedling") | (plot == 10 & age == "adult")) %>% # high density plots
  group_by(site, treatment, age) %>%
  summarise(seeds = mean(seeds)) %>% # plot average
  ungroup() %>%
  mutate(plot_seeds = case_when(age == "adult" ~ seeds * 8,
                                age == "seedling" ~ seeds * 16), # scale up with density
         sp = "Ev") %>%
  group_by(treatment, sp, age) %>%
  summarise(seed_heads = mean(plot_seeds),
            seed_heads_se = sd(plot_seeds)/sqrt(length(plot_seeds))) %>% # experiment average
  ungroup()

# Mv seeds
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, flowers) %>%
  mutate(age = "seedling")  %>%
  filter(plot == 4) %>% # high density plots
  group_by(site, treatment, sp, age) %>%
  summarise(flowers = mean(flowers, na.rm = T)) %>% # plot average
  ungroup() %>%
  mutate(plot_flowers = flowers * 64) %>% # scale up with density
  group_by(treatment, sp, age) %>%
  summarise(seed_heads = mean(plot_flowers),
            seed_heads_se = sd(plot_flowers)/sqrt(length(plot_flowers))) %>% # experiment average
  ungroup()

# plot biomass
# use average of others in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         sp_age = paste(sp, age, sep = "_")) %>%
  filter((plot == 4 & sp == "Mv") | (plot == 7 & sp_age == "Ev_seedling") | (plot == 10 & sp_age == "Ev_adult")) %>% # high density plots
  left_join(mvBioD2Dat %>% # add focal biomass
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g),
                     age = "seedling") %>%
              group_by(site, plot, treatment, sp, age) %>%
              summarise(biomass_foc = sum(biomass_weight.g)) %>%
              ungroup() %>%
              full_join(evBioD2Dat %>%
                          filter(ID %in% c("1", "2", "3")) %>%
                          group_by(site, plot, treatment) %>%
                          mutate(weight_adj = mean(weight, na.rm = T)) %>%
                          ungroup() %>%
                          mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                                    TRUE ~ weight),
                                 sp = "Ev",
                                 age = "seedling") %>%
                          group_by(site, plot, treatment, sp, age) %>%
                          summarise(biomass_foc = sum(weight)) %>%
                          ungroup()) %>%
              full_join(evBioD2Dat %>%
                          filter(ID == "A") %>%
                          select(site, plot, treatment, weight) %>%
                          rename(biomass_foc = weight) %>%
                          mutate(sp = "Ev",
                                 age = "adult"))) %>%
  mutate(biomass_tot = biomass.g + biomass_foc) %>% # plot value
  group_by(treatment, sp, age, sp_age) %>%
  summarise(biomass_g_m2 = mean(biomass_tot),
            biomass_g_m2_se = sd(biomass_tot)/sqrt(length(biomass_tot))) %>% # experiment average
  ungroup()

# severity
plotSevD2Dat <- sevD2Dat %>%
  mutate(sp_age = paste(sp, age, sep = "_"),
         severity = 100 * (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 100, 100, severity)) %>%
  filter(month == "early_aug") %>%
  filter((plot == 4 & sp == "Mv") | (plot == 7 & sp_age == "Ev_seedling") | (plot == 10 & sp_age == "Ev_adult")) %>%
  group_by(site, treatment, sp, age, sp_age) %>%
  summarise(severity = mean(severity, na.rm = T)) %>% # plot average
  ungroup()%>%
  group_by(treatment, sp, age) %>%
  summarise(severity_se = sd(severity)/sqrt(length(severity)),
            severity = mean(severity)) %>% # experiment average
  ungroup()

# combine data
dat <- evSeedD2Dat2 %>%
  full_join(mvSeedD2Dat2) %>%
  mutate(sp_age = paste(sp, age, sep = "_"),
         year = 2019,
         site = 1) %>%
  full_join(plotBioD2Dat) %>%
  full_join(plotSevD2Dat) %>%
  mutate(experiment = "current study") %>%
  full_join(flory %>%
              mutate(sp = "Mv",
                     age = "seedling",
                     sp_age = paste(sp, age, sep = "_"),
                     experiment = "Flory et al. 2011",
                     year = 2010)) %>%
  full_join(stricker %>%
              mutate(sp = "Mv",
                     age = "seedling",
                     sp_age = paste(sp, age, sep = "_"),
                     experiment = "Stricker et al. 2016")) %>%
  mutate(year_site = paste(year, site, sep = "_"),
         treatment = case_when(treatment == "control" ~ "water",
                               TRUE ~ treatment))


#### effect size ####

# pooled sd function
pooled_sd <- function(s1, s2){
  out <- sqrt((s1^2 + s2^2)/2)
  return(out)
}

# need sample sizes for cohen's d to be more accurate
dat2 <- dat %>%
  pivot_wider(names_from = treatment,
              values_from = c(severity, severity_se, biomass_g_m2, biomass_g_m2_se, seed_heads, seed_heads_se),
              names_glue = "{.value}_{treatment}") %>%
  mutate(biomass_g_m2_pooled_sd = pooled_sd(biomass_g_m2_se_fungicide, biomass_g_m2_se_water),
         seed_heads_pooled_sd = pooled_sd(seed_heads_fungicide, seed_heads_water),
         biomass_g_m2_d = (biomass_g_m2_fungicide - biomass_g_m2_water) / biomass_g_m2_pooled_sd,
         seed_heads_d = (seed_heads_fungicide - seed_heads_water) / seed_heads_pooled_sd,
         biomass_g_m2_eff = (biomass_g_m2_fungicide - biomass_g_m2_water) / biomass_g_m2_water,
         seed_heads_eff = (seed_heads_fungicide - seed_heads_water) / seed_heads_water,
         severity_eff = (severity_fungicide - severity_water) / severity_water)

# mv only
mvDat2 <- dat2 %>%
  filter(sp == "Mv")


#### figure ####

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

study_col = "#FDE725FF"

# distributions
ggplot(mvDat2, aes(x = severity_eff)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = filter(mvDat2, experiment == "current study")$severity_eff, color = study_col) +
  xlab("Proportional change in lesions (%)") +
  ylab("Replicates") +
  fig_theme

ggplot(mvDat2, aes(x = biomass_g_m2_eff)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = filter(mvDat2, experiment == "current study")$biomass_g_m2_eff, color = study_col) +
  xlab(expression(paste("Proportional change in biomass (g ", m^-2, ")", sep = ""))) +
  ylab("Replicates") +
  fig_theme

ggplot(mvDat2, aes(x = seed_heads_eff)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = filter(mvDat2, experiment == "current study")$seed_heads_eff, color = study_col) +
  xlab("Proportional change in seed heads") +
  ylab("Replicates") +
  fig_theme

ggplot(dat2, aes(severity_water, biomass_g_m2_eff, group = year_site, color = experiment, shape = sp_age)) +
  geom_point()

bio_fig <- dat2 %>%
  filter(sp == "Mv") %>%
  ggplot(aes(severity_water, biomass_g_m2_eff, group = year_site, color = experiment)) +
  geom_point() +
  theme_bw() +
  fig_theme + 
  theme(legend.position = c(0.65, 0.83)) +
  xlab("Lesion area without fungicide (%)") +
  ylab("Fungicide effect on biomass")

ggplot(dat2, aes(severity_water, seed_heads_eff, group = year_site, color = experiment, shape = sp_age)) +
  geom_point()

seed_fig <- dat2 %>%
  filter(sp == "Mv") %>%
  ggplot(aes(severity_water, seed_heads_eff, group = year_site, color = experiment)) +
  geom_point() +
  theme_bw() +
  fig_theme +
  theme(axis.title.y = element_text(size = 10, hjust = 0)) +
  xlab("Lesion area without fungicide (%)") +
  ylab("Fungicide effect on seed heads")

pdf("output/plot_biomass_seeds_severity_all_exp.pdf",
    width = 5, height = 2.5)
plot_grid(bio_fig, seed_fig,
          labels = c("A", "B"),
          nrow = 1)
dev.off()

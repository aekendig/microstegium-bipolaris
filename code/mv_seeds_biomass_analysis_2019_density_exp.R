##### info ####

# file: mv_seeds_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 2/10/20
# goal: evaluate the effects of density treatments and environmental covariates on the biomass and seed production of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
mv_bio <- read_csv("./data/focal_mv_seeds_biomass_2019_density_exp.csv")
mv_bag <- read_csv("./data/focal_mv_bag_seeds_2019_density_exp.csv")
plots <- read_csv("./data/plot_treatments_for_figures_2018_2019_density_exp.csv")


#### edit data ####

# look at notes
unique(mv_bag$process_notes)

# look at seed counts
unique(mv_bag$seeds)

# average seed count by plant
seed_count <- mv_bag %>%
  group_by(site, plot, treatment, sp, plant) %>%
  summarize(seeds_per_flower = mean(seeds, na.rm = T))
  
# combine with plot data
# add three to all density measurements (three focals)
micro <- mv_bio %>%
  left_join(plots) %>%
  left_join(seed_count) %>%
  mutate(background_density_tot = case_when(background == "Ev adult" ~ background_density + 1,
                   TRUE ~ background_density + 3),
         flower_seeds = flowers * seeds_per_flower)

#### raw data figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# color palettes
col_pal <- c("white", "black")

# base summary figure
base_fig <- ggplot(micro, aes(x = background_density_tot)) +
  stat_summary(aes(group = treatment), geom = "errorbar", fun.data = mean_cl_boot, width = 5, position = position_dodge(2)) +
  stat_summary(aes(fill = treatment), geom = "point", fun.y = mean, position = position_dodge(2), shape = 21, size = 3) +
  facet_wrap(~background, scales = "free") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  xlab("Plant density") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = lg_txt))

# base scatterplot figure
scat_fig <- ggplot(micro, aes(x = biomass_weight.g, y = flower_seeds)) +
  geom_point(shape = 21, size = 2, aes(fill = treatment)) +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none")

# biomass figure
base_fig %+% aes(y = biomass_weight.g) +
  theme(legend.position = c(0.9, 0.8)) +
  ylab("Biomass (g)")

# seed figure
base_fig %+% aes(y = flower_seeds) +
  theme(legend.position = c(0.9, 0.8)) +
  ylab("Chasmogamous seeds")

# biomass and seed relationship
scat_fig +
  xlab("Biomass (g)") +
  ylab("Chasmogamous seeds")


#### output intermediate data ####
write_csv(micro, "intermediate-data/mv_seeds_biomass_covariates_2019_density_exp.csv")

##### info ####

# file: mv_seeds_biomass_analysis_2018_density_exp
# author: Amy Kendig
# date last edited: 2/13/20
# goal: evaluate the effects of density treatments and environmental covariates on the biomass and seed production of Microstegium


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
mv_bio <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")
mv_lit <- read_csv("./data/mv_litter_biomass_apr_2019_density_exp.csv")
mv_bag <- read_csv("./data/mv_bag_seed_2018_density_exp.csv")
plots <- read_csv("./data/plot_treatments_for_figures_2018_2019_density_exp.csv")


#### edit data ####

# look at notes
unique(mv_bio$processing_notes)

# edit biomass values and join with other data
micro <- mv_bio %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment)) %>%
  left_join(mv_lit %>%
              select(site, plot, treatment, mv.g) %>%
              rename(litter.g = mv.g)) %>%
  left_join(plots) %>%
  mutate(background_density_tot = background_density + 3,
         litter_conversion = litter.g/bio.g)


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

# litter conversion
ggplot(micro, aes(x = treatment, y = litter_conversion)) +
  geom_boxplot()


#### statistics ####

# t-test of treatment on litter conversion
t.test(litter_conversion ~ treatment, data = filter(micro, !is.na(litter_conversion)))

#### output intermediate data ####
write_csv(micro, "intermediate-data/mv_seeds_biomass_covariates_2018_density_exp.csv")

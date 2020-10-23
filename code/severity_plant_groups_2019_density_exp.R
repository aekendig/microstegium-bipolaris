##### info ####

# file: severity_plant_groups_2019_density_exp
# author: Amy Kendig
# date last edited: 10/19/20
# goal: demonstrate effect of fungicide treatments on severity of each plant group and on Mv outside of plots


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
plot_dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")


#### edit data ####

# combine data
# assign plant group
# calculate severity 
dat <- plot_dat %>%
  full_join(edge_dat) %>%
  mutate(plant_group = case_when(ID == "Edge" ~ "Mv outside plot",
                                 sp == "Mv" & ID != "Edge" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "Ev seedling",
                                 sp == "Ev" & age == "adult" ~ "Ev adult"),
         severity = lesion_area.pix / leaf_area.pix)


#### figure ####
pdf("output/severity_plant_groups_2019_density_exp.pdf", width = 5, height = 3)
dat %>%
  filter(month == "early_aug") %>%
  ggplot(aes(plant_group, severity, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.2)) +
  xlab("Plant group") +
  ylab ("Proportion of leaf surface with lesions") +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, hjust = 0.9),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.8, 0.8))
dev.off()

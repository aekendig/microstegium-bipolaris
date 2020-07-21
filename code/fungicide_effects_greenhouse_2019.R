##### info ####

# file: fungicide_effects_greenhouse_2019
# author: Amy Kendig
# date last edited: 7/6/20
# goal: see how fungicide treatment affects Mv biomass and seed head production, Ev biomass


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
mv_dat <- read_csv("./data/mv_biomass_seeds_height_jun_2019_fungicide_exp.csv")
ev_dat <- read_csv("./data/ev_biomass_dec_2019_fungicide_exp.csv")


#### edit data ####

# look at notes
unique(mv_dat$notes) # dead plant
unique(ev_dat$notes)

# remove the dead plant
# fungicide column
mv_dat2 <- mv_dat %>%
  filter(is.na(notes)) %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0"))

# combine live and dead weight
# fungicide column
ev_dat2 <- ev_dat %>%
  group_by(treatment, pot, sp) %>%
  summarise(weight.g = sum(weight.g)) %>%
  ungroup() %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0"))


#### visualize ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        plot.title = element_text(size = 10, hjust = 0.5))

# save to pdf
pdf("output/fungicide_effects_greenhouse_2019.pdf")

# Mv biomass
mv_dat2 %>%
  ggplot(aes(x = treatment, y = weight.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " aboveground biomass (g)", sep = ""))) +
  xlab("Treatment") +
  temp_theme
# fungicide reduced weight

# Ev biomass
ev_dat2 %>%
  ggplot(aes(x = treatment, y = weight.g)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Elymus), " aboveground biomass (g)", sep = ""))) +
  xlab("Treatment") +
  temp_theme
# fungicide slightly increased weight

# seeds
mv_dat2 %>%
  ggplot(aes(x = treatment, y = seed_heads)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " seed heads", sep = ""))) +
  xlab("Treatment") +
  temp_theme
# fungicide reduced seed heads

# height
mv_dat2 %>%
  ggplot(aes(x = treatment, y = height.in)) +
  geom_point(color = "gray") +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", color = "black", width = 0.05) +
  ylab(expression(paste(italic(Microstegium), " height (in.)", sep = ""))) +
  xlab("Treatment") +
  temp_theme
# no effect of fungicide on height

# biomass/seeds

mv_dat2 %>%
  ggplot(aes(x = weight.g, y = seed_heads, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab(expression(paste(italic(Microstegium), " aboveground biomass (g)", sep = ""))) +
  ylab(expression(paste(italic(Microstegium), " seed heads", sep = ""))) +
  temp_theme

# close pdf
dev.off()

#### statistical models ####

# Mv biomass
mv_fung_bio_mod <- t.test(weight.g ~ fungicide, data = mv_dat2)
mv_fung_bio_mod
# fungicide mean: 21.05
# water mean: 23.80
# no significant effect of fungicide

# Ev biomass
ev_fung_bio_mod <- t.test(weight.g ~ fungicide, data = ev_dat2)
ev_fung_bio_mod
# fungicide mean: 13.84
# water mean: 13.39
# no significant effect of fungicide
  
# Mv seed heads
mv_fung_seed_mod <- t.test(seed_heads ~ fungicide, data = mv_dat2)
mv_fung_seed_mod
# fungicide mean: 30
# water mean: 37
# no significant effect of fungicide

# Mv height
mv_fung_height_mod <- t.test(height.in ~ fungicide, data = mv_dat2)
mv_fung_height_mod
# fungicide mean: 42.2
# water mean: 41.2
# no significant effect of fungicide

# Seeds by biomass
mv_seed_bio_mod <- lm(seed_heads ~ weight.g * fungicide, data = mv_dat2)
summary(mv_seed_bio_mod)
# water slope: 0.64
# fungicide slope: 0.82
# no significant effect of fungicide
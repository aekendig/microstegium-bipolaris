##### info ####

# file: mv_density_figure_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/24/21
# goal: figure of invader abundance vs. density and competitor fitness vs. invader density

# for 2018 figure, add annual survival of Ev as a function of Mv density

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)

# import models
load("output/focal_seed_density_model_2019_density_exp.rda")
load("output/focal_growth_density_model_2019_density_exp.rda")
load("output/mv_plot_biomass_density_model_2019_dens_exp.rda")
load("output/mv_plot_seed_density_model_2019_dens_exp.rda")

# import data
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

seedD2Dat <- read_csv("intermediate-data/focal_seed_density_data_2019_density_exp.csv")
seedD2Alpha <- read_csv("output/focal_seed_competition_coefficients_2018_2019_density_exp.csv")
# focal_seed_density_2018_2019_density_exp.R

growthD2Dat <- read_csv("intermediate-data/focal_growth_density_data_2019_density_exp.csv")
growthD2Alpha <- read_csv("output/focal_growth_competition_coefficients_2018_2019_density_exp.csv")
# focal_growth_2018_2019_density_exp.R

plotD2Dat <- read_csv("intermediate-data/mv_plot_biomass_seeds_2019_density_exp.csv")
plotD2Alpha <- read_csv("output/mv_plot_biomass_seeds_density_treatment_effect_2019_dens_exp.csv")
# mv_plot_biomass_seeds_density_2019_density_exp.R


#### edit data ####

# plant group densities
plotDens <- plots %>%
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1,
                             TRUE ~ background_density)) %>%
  select(plot, treatment, background, density)

# plot biomass/seed data
plotD2Dat2 <- plotD2Dat %>%
  left_join(plotDens) %>%
  filter(background %in% c("none", "Mv seedling")) %>%
  rename(biomass = biomass_mv) %>%
  mutate(focal = "Invader",
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)")) %>%
  select(-c(biomass_bg, biomass_foc_mv)) %>%
  pivot_longer(cols = c(seeds, biomass),
               names_to = "response",
               values_to = "value")

# raw data
rawD2Dat <- seedD2Dat %>%
  select(site, plot, treatment, fungicide, focal, foc, ID, background, bg,density, seeds, log_seeds) %>%
  mutate(response = "seeds") %>%
  rename(fitness = seeds,
         log_fitness = log_seeds) %>%
  full_join(growthD2Dat %>%
              select(site, plot, treatment, fungicide, focal, foc, ID, background, bg,density, biomass_weight.g, plant_growth) %>%
              mutate(response = "growth") %>%
              rename(fitness = biomass_weight.g,
                     log_fitness = plant_growth)) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         focal = fct_recode(focal, "First-year competitor" = "Ev seedling",
                            "Adult competitor" = "Ev adult",
                            "Invader" = "Mv") %>%
           fct_relevel("Invader"),
         background = fct_recode(background, "First-year competitor" = "Ev seedling",
                                 "Adult competitor" = "Ev adult",
                                 "Invader" = "Mv") %>%
           fct_relevel("Invader")) %>%
  filter(bg == "m" & foc != "m")

# template prediction data
predDatTemplate <- rawD2Dat %>%
  group_by(focal, foc, background, bg, treatment, fungicide) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

plotPredDatTemplate <- plotDens %>%
  filter(background %in% c("none", "Mv seedling")) %>%
  group_by(treatment) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

# prediction data
seedPredD2Dat <- predDatTemplate %>%
    mutate(log_fitness = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

growthPredD2Dat <- predDatTemplate %>%
  mutate(log_fitness = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

plotSeedPredD2Dat <- plotPredDatTemplate %>%
  mutate(value = fitted(mvSeedDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(mvSeedDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(mvSeedDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

plotBioPredD2Dat <- plotPredDatTemplate %>%
  mutate(value = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

# combine
predD2Dat <- seedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(growthPredD2Dat %>%
              mutate(response = "growth")) %>%
  mutate(fitness = exp(log_fitness))

plotPredD2Dat <- plotSeedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(plotBioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))


#### figure settings ####

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}

# dodge size
dodge_width <- 3


#### figure ####

# labels
plot_labels <- c(biomass = "Plot~biomass~(g~m^-2)",
                 seeds = "Plot~seeds~(m^-2)")


response_labels <- c(growth = "Plant growth (ln[biomass])",
                     seeds = "Seeds (ln[seeds + 1])")

# invader figure
inv_fig <- ggplot(plotPredD2Dat, aes(x = density, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  stat_summary(data = plotD2Dat2, geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(dodge_width), aes(color = treatment)) +
  stat_summary(data = plotD2Dat2, geom = "point", fun = "mean", size = 1.5, position = position_dodge(dodge_width), aes(color = treatment)) +
  facet_grid(rows = vars(response),
             cols = vars(focal),
             scales = "free",
             switch = "y",
             labeller = labeller(response = as_labeller(plot_labels, label_parsed))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = "none")

# competitor figure
comp_fig <- ggplot(predD2Dat, aes(x = density, y = log_fitness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  stat_summary(data = rawD2Dat, geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(dodge_width), aes(color = treatment)) +
  stat_summary(data = rawD2Dat, geom = "point", fun = "mean", size = 1.5, position = position_dodge(dodge_width), aes(color = treatment)) +
  facet_grid(rows = vars(response),
             cols = vars(focal),
             switch = "y",
             scales = "free",
             labeller = labeller(response = response_labels)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  fig_theme +
  theme(legend.position = c(0.82, 0.9))

comb_fig <- plot_grid(inv_fig, comp_fig,
                      nrow = 1,
                      rel_widths = c(0.65, 1))

x_title <- textGrob(expression(paste("Invader density (", m^-2, ")", sep = "")), 
                    gp = gpar(fontsize = 10))

pdf("output/mv_density_figure_2019_density_exp.pdf", width = 6, height =4.5)
grid.arrange(arrangeGrob(comb_fig, bottom = x_title))
dev.off()
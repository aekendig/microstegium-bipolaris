##### outputs ####

# Figure 2
# Figure S1
# Tables S1-S2 & S7-S9
# Results paragraphs 1-2


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(patchwork)
library(broom.mixed)

# import models
load("output/focal_seed_density_model_2018_density_exp.rda")
load("output/focal_seed_density_model_2019_density_exp.rda")
load("output/focal_growth_density_model_2018_density_exp.rda")
load("output/focal_growth_density_model_2019_density_exp.rda")
load("output/mv_plot_biomass_density_model_2019_dens_exp.rda")
load("output/mv_plot_seed_density_model_2019_dens_exp.rda")
load("output/mv_litter_establishment_model_2018_litter_exp.rda")
load("output/ev_litter_establishment_bh_model_2019_litter_exp.rda")

# import data
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

seedD1Dat <- read_csv("intermediate-data/focal_seed_density_data_2018_density_exp.csv")
seedD2Dat <- read_csv("intermediate-data/focal_seed_density_data_2019_density_exp.csv")

growthD1Dat <- read_csv("intermediate-data/focal_growth_density_data_2018_density_exp.csv")
growthD2Dat <- read_csv("intermediate-data/focal_growth_density_data_2019_density_exp.csv")

plotD2Dat <- read_csv("intermediate-data/mv_plot_biomass_seeds_2019_density_exp.csv")

mvLitDat <- read_csv("intermediate-data/mv_litter_establishment_data_2018_litter_exp.csv")
evLitDat <- read_csv("intermediate-data/ev_litter_establishment_data_2018_litter_exp.csv")

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}


#### edit raw data ####

# plant group densities
plotDens <- plots %>%
  select(plot, treatment, background, background_density)

# plot biomass/seed data
plotD2Dat2 <- plotD2Dat %>%
  left_join(plotDens) %>%
  filter(background %in% c("none", "Mv seedling")) %>%
  mutate(density = background_density + 3,
         focal = "Invader (Mv)",
         treatment = fct_relevel(treatment, "water")) %>%
  rename(biomass = biomass_mv) %>%
  select(-c(biomass_bg, biomass_foc_mv)) %>%
  pivot_longer(cols = c(seeds, biomass),
               names_to = "response",
               values_to = "value")

# raw data
rawD1Dat <- seedD1Dat %>%
  select(site, plot, treatment, fungicide, focal, foc, ID, background, bg, density, seeds, log_seeds) %>%
  mutate(response = "seeds") %>%
  rename(fitness = seeds,
         log_fitness = log_seeds) %>%
  full_join(growthD1Dat %>%
              select(site, plot, treatment, fungicide, focal, foc, ID, background, bg, density, tillers_jul, tillers_jun, plant_growth) %>%
              mutate(response = "growth") %>%
              rename(fitness = tillers_jul/tillers_jun,
                     log_fitness = plant_growth)) %>%
  mutate(treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control"),
         focal = fct_recode(focal, "First-year competitor (Ev)" = "Ev seedling",
                            "Adult competitor (Ev)" = "Ev adult",
                            "Invader (Mv)" = "Mv") %>%
           fct_relevel("Invader (Mv)", "First-year competitor (Ev)"),
         background = fct_recode(background, "First-year competitor (Ev)" = "Ev seedling",
                                 "Adult competitor (Ev)" = "Ev adult",
                                 "Invader (Mv)" = "Mv") %>%
           fct_relevel("Invader (Mv)", "First-year competitor (Ev)")) %>%
  filter(bg == "m")

rawD2Dat <- seedD2Dat %>%
  select(site, plot, treatment, fungicide, focal, foc, ID, background, bg, density, seeds, log_seeds) %>%
  mutate(response = "seeds") %>%
  rename(fitness = seeds,
         log_fitness = log_seeds) %>%
  full_join(growthD2Dat %>%
              select(site, plot, treatment, fungicide, focal, foc, ID, background, bg, density, biomass_weight.g, plant_growth) %>%
              mutate(response = "growth") %>%
              rename(fitness = biomass_weight.g,
                     log_fitness = plant_growth)) %>%
  mutate(treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control"),
         focal = fct_recode(focal, "First-year competitor (Ev)" = "Ev seedling",
                            "Adult competitor (Ev)" = "Ev adult",
                            "Invader (Mv)" = "Mv") %>%
           fct_relevel("Invader (Mv)", "First-year competitor (Ev)"),
         background = fct_recode(background, "First-year competitor (Ev)" = "Ev seedling",
                                 "Adult competitor (Ev)" = "Ev adult",
                                 "Invader (Mv)" = "Mv") %>%
           fct_relevel("Invader (Mv)", "First-year competitor (Ev)")) %>%
  filter(bg == "m")

mvLitDat2 <- mvLitDat %>%
  mutate(focal = "Invader (Mv)")

evLitDat2 <- evLitDat %>%
  mutate(focal = "Competitor (Ev)")


#### edit prediction data ####

# template prediction data
predDatTemplate <- rawD2Dat %>%
  group_by(focal, foc, background, bg, treatment, fungicide) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

plotPredDatTemplate <- plotD2Dat2 %>%
  group_by(treatment) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

# prediction data
seedPredD1Dat <- predDatTemplate %>%
  mutate(log_fitness = fitted(seedD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(seedD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(seedD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

seedPredD2Dat <- predDatTemplate %>%
    mutate(log_fitness = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(seedD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

growthPredD1Dat <- predDatTemplate %>%
  mutate(log_fitness = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

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

mvEstPredDat <- tibble(litter.g.cm2 = rep(seq(min(mvLitDat$litter.g.cm2), max(mvLitDat$litter.g.cm2), 
                                              length.out = 100), 2),
                       sterilized = rep(c(0, 1), each = 100)) %>%
  mutate(mv_planted_cor = round(mean(mvLitDat$mv_planted_cor))) %>%
  mutate(Est = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Estimate"]/mv_planted_cor,
         Est_lower = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Q2.5"]/mv_planted_cor,
         Est_upper = fitted(mvEstL1Mod2, newdata = ., type = "response", re_formula = NA)[, "Q97.5"]/mv_planted_cor,
         litter.g.m2 = litter.g.cm2 * 10000,
         sterilizedF = ifelse(sterilized == 0, "live", "sterilized"),
         focal = "Invader (Mv)")

evEstPredDat <- tibble(litter.g.cm2 = seq(min(evLitDat$litter.g.cm2), max(evLitDat$litter.g.cm2), 
                                              length.out = 100)) %>%
  mutate(Est = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Estimate"],
         Est_lower = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q2.5"],
         Est_upper = fitted(evEstL2BhMod2, newdata = ., type = "response", re_formula = NA)[, "Q97.5"],
         litter.g.m2 = litter.g.cm2 * 10000,
         focal = "Competitor (Ev)")

# combine
predD1Dat <- seedPredD1Dat %>%
  mutate(response = "seeds") %>%
  full_join(growthPredD1Dat %>%
              mutate(response = "growth")) %>%
  mutate(fitness = exp(log_fitness))

predD2Dat <- seedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(growthPredD2Dat %>%
              mutate(response = "growth")) %>%
  mutate(fitness = exp(log_fitness))

plotPredD2Dat <- plotSeedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(plotBioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  mutate(focal = "Invader (Mv)")


#### figure settings ####

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                  margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

col_pal = c("black", "#238A8D")
col_pal2 = c("black", "#8F4BED")

# dodge size
dodge_width <- 3


#### figure ####

# labels
plot_labels <- c(biomass = "Biomass~(g~m^-2)",
                 seeds = "Seed~production~(m^-2)")

response_labels <- c(growth = "Biomass (ln[g])",
                     seeds = "Seed production (ln[x+1])")

response_labels2 <- c(growth = "Tiller growth (ln[x])",
                     seeds = "Seed production (ln[x+1])")

focal_labels <- c("Invader (Mv)" = "Invader (Mv)",
                  "Adult competitor (Ev)" = "Adult competitor (Ev)",
                  "First-year competitor (Ev)" = "1st yr competitor (Ev)")

# abundance figure
inv_fig <- ggplot(plotPredD2Dat, aes(x = density, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = plotD2Dat2, aes(color = treatment), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  facet_grid(rows = vars(response),
             cols = vars(focal),
             scales = "free",
             switch = "y",
             labeller = labeller(response = as_labeller(plot_labels, label_parsed))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(x = expression(paste("Invader density (", m^-2, ")", sep = "")),
       title = "Invader abundance") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# per capita figure
percap_fig <- ggplot(predD2Dat, aes(x = density, y = log_fitness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = rawD2Dat, aes(color = treatment), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  facet_grid(rows = vars(response),
             cols = vars(focal),
             switch = "y",
             scales = "free",
             labeller = labeller(response = response_labels,
                                 focal = focal_labels)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  labs(x = expression(paste("Invader density (", m^-2, ")", sep = "")),
       title = "Invader per capita effects") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# litter figure
mvLitFig <- ggplot(mvLitDat2, aes(x = litter.g.m2, y = 100*prop_germ_den_cor, fill = sterilizedF, color = sterilizedF)) +
  geom_ribbon(data = mvEstPredDat, aes(y = 100*Est, ymin = 100*Est_lower, ymax = 100*Est_upper), alpha = 0.3, color = NA) +
  geom_line(data = mvEstPredDat, aes(y = 100*Est)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  scale_color_manual(values = col_pal2) +
  scale_fill_manual(values = col_pal2) +
  facet_grid(cols = vars(focal)) +
  labs(y = "Seedling emergence (%)",
       title = "Invader litter effects") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

evLitFig <- ggplot(evLitDat2, aes(x = litter.g.m2, y = 100*ev_prop_germ)) +
  geom_ribbon(data = evEstPredDat, aes(y = 100*Est, ymin = 100*Est_lower, ymax = 100*Est_upper), alpha = 0.3, color = NA) +
  geom_line(data = evEstPredDat, aes(y = 100*Est)) +
  geom_point(size = 0.5) +
  facet_grid(cols = vars(focal)) +
  labs(x = expression(paste("Invader litter (g ", m^-2, ")", sep = "")),
       y = "Seedling emergence (%)") +
  fig_theme

# legend figure
leg_dat <- tibble(x = rep(seq(0, 2), 3),
                  y = rep(seq(1, 3), 3),
                  lower = rep(seq(0, 2), 3),
                  upper = rep(seq(2, 4), 3),
                  trt = rep(c("control", "fungicide", "autoclaved"), each = 3)) %>%
  mutate(trt = fct_relevel(trt, "control", "fungicide"))

leg_fig <- ggplot(leg_dat, aes(x = x, y = y, color = trt, fill = trt)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_point(size = 2) +
  scale_color_manual(values = c(col_pal, col_pal2[2]), name = "Disease treatment") +
  scale_fill_manual(values = c(col_pal, col_pal2[2]), name = "Disease treatment") +
  fig_theme

leg <- get_legend(leg_fig)

comb_fig <- (inv_fig + percap_fig +
               (mvLitFig + (evLitFig + plot_layout(tag_level = 'new')) 
                + plot_layout(ncol = 1)) +
               plot_layout(width = c(0.36, 1, 0.35)) +
               plot_annotation(tag_levels = 'A') & 
               theme(plot.tag = element_text(size = 10, face = "bold")))
   
pdf("output/mv_density_figure_2019_density_exp.pdf", width = 7.08, height = 3.54)
plot_grid(comb_fig, leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### supplementary figure ####

# competitor figure
pdf("output/mv_per_capita_effects_figure_2018_density_exp.pdf", width = 4, height = 3.15)
ggplot(predD1Dat, aes(x = density, y = log_fitness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = rawD1Dat, aes(color = treatment), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  facet_grid(rows = vars(response),
             cols = vars(focal),
             switch = "y",
             scales = "free",
             labeller = labeller(response = response_labels2,
                                 focal = focal_labels)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  labs(x = expression(paste("Invader density (", m^-2, ")", sep = "")),
       title = "Invader per capita effects") +
  fig_theme +
  theme(axis.title.y = element_blank())
dev.off()


#### values for text ####

# Mv biomass
set.seed(184)
posterior_predict(mvBioDensMod,
                  newdata = filter(plotPredDatTemplate, density == 67),
                  allow_new_levels = T) %>%
  as_tibble(.name_repair = ~ c("water", "fungicide")) %>%
  mutate(fung_eff = 100 * (fungicide - water) / water) %>%
  median_hdi(fung_eff)

# evA seeds and biomass
as_draws_df(growthD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(alpha_fungicide = (density + foca_density + fungicide_density + foca_fungicide_density),
            alpha_water = (density + foca_density)) %>%
  mutate(fung_eff = 100 * (alpha_fungicide - alpha_water) / alpha_water,
         draws = 1:15000) %>%
  pivot_longer(cols = c(alpha_fungicide, alpha_water, fung_eff),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)

as_draws_df(seedD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(alpha_fungicide = (density + foca_density + fungicide_density + foca_fungicide_density),
            alpha_water = (density + foca_density)) %>%
  mutate(fung_eff = 100 * (alpha_fungicide - alpha_water) / alpha_water,
         draws = 1:15000) %>%
  pivot_longer(cols = c(alpha_fungicide, alpha_water, fung_eff),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)

growthD1Table <- as_draws_df(growthD1Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  mutate(alpha_fungicide_a = (density + foca_density + fungicide_density + foca_fungicide_density),
         alpha_water_a = (density + foca_density),
         alpha_fungicide_s = (density + focs_density + fungicide_density + focs_fungicide_density),
         alpha_water_s = (density + focs_density),
         alpha_fungicide_m = (density + fungicide_density),
         alpha_water_m = density) %>%
  transmute(EvA_fung_eff = 100 * (alpha_fungicide_a - alpha_water_a) / alpha_water_a,
            EvS_fung_eff = 100 * (alpha_fungicide_s - alpha_water_s) / alpha_water_s,
            Mv_fung_eff = 100 * (alpha_fungicide_m - alpha_water_m) / alpha_water_m,
            draws = 1:15000) %>%
  pivot_longer(cols = c(EvA_fung_eff, EvS_fung_eff, Mv_fung_eff),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)

seedD1Table <- as_draws_df(seedD1Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  mutate(alpha_fungicide_a = (density + foca_density + fungicide_density + foca_fungicide_density),
         alpha_water_a = (density + foca_density),
         alpha_fungicide_s = (density + focs_density + fungicide_density + focs_fungicide_density),
         alpha_water_s = (density + focs_density),
         alpha_fungicide_m = (density + fungicide_density),
         alpha_water_m = density) %>%
  transmute(EvA_fung_eff = 100 * (alpha_fungicide_a - alpha_water_a) / alpha_water_a,
            EvS_fung_eff = 100 * (alpha_fungicide_s - alpha_water_s) / alpha_water_s,
            Mv_fung_eff = 100 * (alpha_fungicide_m - alpha_water_m) / alpha_water_m,
            draws = 1:15000) %>%
  pivot_longer(cols = c(EvA_fung_eff, EvS_fung_eff, Mv_fung_eff),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)

perCapD1Table <- growthD1Table %>%
  mutate(response = "tiller growth") %>%
  full_join(seedD1Table %>%
              mutate(response = "seeds")) %>%
  select(response, variable, value, .lower, .upper)


#### tables ####

write_csv(tidy(mvBioDensMod), "output/mv_plot_biomass_density_model_2019_dens_exp.csv")
write_csv(tidy(mvSeedDensMod), "output/mv_plot_seed_density_model_2019_dens_exp.csv")
write_csv(tidy(mvEstL1Mod2), "output/microstegium_litter_establishment_model_2018_litter_exp.csv")
write_csv(tidy(evEstL2BhMod2), "output/elymus_litter_establishment_bh_model_2019_litter_exp.csv")
write_csv(perCapD1Table, "output/mv_per_capita_effects_2018_dens_exp.csv")
##### info ####

# file: focal_growth_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/5/21
# goal: analyses of plant growth as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)

# import plot information
plotsD <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")

# import growth data
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1)) %>%
  select(plot, treatment, background, density)

# missing data
filter(tillerD1Dat, is.na(tillers_jul)) # 0
filter(tillerD1Dat, is.na(tillers_jun)) # 0
filter(tillerD1Dat, tillers_jul == 0 | tillers_jun == 0) # these plants are dead
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling

# combine data
growthD1Dat <- tillerD1Dat %>%
  left_join(plotDens) %>%
  mutate(plant_growth = log(tillers_jul/tillers_jun),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(plot, str_sub(treatment, 1, 1), sep = "")) %>%
  filter(tillers_jun > 0 & tillers_jul > 0)

growthD2Dat <- mvBioD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  left_join(plotDens) %>%
  mutate(plant_growth = log(biomass_weight.g),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(plot, str_sub(treatment, 1, 1), sep = "")) %>%
  filter(!is.na(biomass_weight.g))


#### models ####

# initial visualization
ggplot(growthD1Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

ggplot(growthD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

# remove plot 1
growthD1Dat2 <- growthD1Dat %>%
  filter(plot > 1)

growthD2Dat2 <- growthD2Dat %>%
  filter(plot > 1)

# fit models
growthD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                   plant_growth ~ density * foc * bg * fungicide + (1|site),
                   prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3) 
# 120 divergent transitions and 9 transitions exceed max treedepth
growthD1Mod <- update(growthD1Mod, control = list(adapt_delta = 0.999999, max_treedepth = 15))
mod_check_fun(growthD1Mod)

growthD2Mod <- update(growthD1Mod, 
                       control = list(adapt_delta = 0.9999999, max_treedepth = 15), 
                       newdata = growthD2Dat2)
mod_check_fun(growthD2Mod)

# save models
save(growthD1Mod, file = "output/focal_growth_density_model_2018_density_exp.rda")
save(growthD2Mod, file = "output/focal_growth_density_model_2019_density_exp.rda")


#### intraspecific vs. interspecific ####

# common terms on both sides of = were deleted
# inter listed first in name
# intra on left side of =

# Mv 
evS_mv_ctrl_hyp = "-density:focs = 0"
evS_mv_fung_hyp = "-density:focs - density:focs:fungicide = 0"
evA_mv_ctrl_hyp = "-density:foca = 0"
evA_mv_fung_hyp = "-density:foca - density:foca:fungicide = 0"

# EvS 
mv_evS_ctrl_hyp = "0.5 * (density:focs + density:focs:bgs + density:foca + density:foca:bgs) = 0"
mv_evS_fung_hyp = "0.5 * (density:focs + density:focs:bgs + density:focs:fungicide + density:focs:bgs:fungicide + density:foca + density:foca:bgs + density:foca:fungicide + density:foca:bgs:fungicide) = 0"

# EvA 
mv_evA_ctrl_hyp = "0.5 * (density:focs + density:focs:bga + density:foca + density:foca:bga) = 0"
mv_evA_fung_hyp = "0.5 * (density:focs + density:focs:bga + density:focs:fungicide + density:focs:bga:fungicide + density:foca + density:foca:bga + density:foca:fungicide + density:foca:bga:fungicide) = 0"

growthD1hyps <- hypothesis(growthD1Mod, 
                           c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                             mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))
# none are significantly different from zero

growthD2hyps <- hypothesis(growthD2Mod, 
                           c(evS_mv_ctrl_hyp, evS_mv_fung_hyp, evA_mv_ctrl_hyp, evA_mv_fung_hyp,
                             mv_evS_ctrl_hyp, mv_evS_fung_hyp, mv_evA_ctrl_hyp, mv_evA_fung_hyp))
# none are significantly different from zero

write_csv(growthD1hyps[[1]], "output/focal_growth_intra_vs_intra_comp_2018_density_exp.csv")
write_csv(growthD2hyps[[1]], "output/focal_growth_intra_vs_intra_comp_2019_density_exp.csv")


#### intercepts ####

# are Ev adult intercepts > Mv in 2019 (they have larger alphas)
# common terms on both sides are deleted

mv_ctrl_int <- "bga = 0"
mv_fung_int <- "bga + bga:fungicide = 0"
evS_ctrl_int <- "bga + focs:bga = 0"
evS_fung_int <- "bga + focs:bga + bga:fungicide + focs:bga:fungicide = 0"
evA_ctrl_int <- "bga + foca:bga = 0"
evA_fung_int <- "bga + foca:bga + bga:fungicide + foca:bga:fungicide = 0"

intD2hyps <- hypothesis(growthD2Mod, 
                           c(mv_ctrl_int, mv_fung_int,
                             evS_ctrl_int, evS_fung_int,
                             evA_ctrl_int, evA_fung_int))

write_csv(intD2hyps[[1]], "output/evA_mv_competition_intercept_comparison_2019_density_exp.csv")


#### figure ####

# density function
dens_fun <- function(f, b, trt, yr){
  
  if(yr == "2018"){
    dat <- growthD1Dat2 %>% filter(foc == f & bg == b & treatment == trt)
  }else{
    dat <- growthD2Dat2 %>% filter(foc == f & bg == b & treatment == trt)
  }
  
  density = seq(min(dat$density), max(dat$density), length.out = 100)
  
  return(density)
}

# predicted data
predD1Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2018",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  mutate(plant_growth = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predD2Dat <- tibble(foc = c("a", "s", "m")) %>%
  expand_grid(tibble(bg = c("a", "s", "m"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  mutate(year = "2019",
         plotf = "A",
         fungicide = case_when(treatment == "water" ~ 0,
                               treatment == "fungicide" ~ 1)) %>%
  mutate(density = pmap(list(foc, bg, treatment, year), dens_fun)) %>%
  unnest(density) %>%
  mutate(plant_growth = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         background = fct_recode(bg, Mv = "m", "Ev seedling" = "s", "Ev adult" = "a"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

# Mv background
mv_mv_ctrl_alpha = "density = 0"
mv_mv_fung_alpha = "density + density:fungicide = 0"
evS_mv_ctrl_alpha = "density + density:focs = 0"
evS_mv_fung_alpha = "density + density:fungicide + density:focs + density:focs:fungicide = 0"
evA_mv_ctrl_alpha = "density + density:foca = 0"
evA_mv_fung_alpha = "density + density:fungicide + density:foca + density:foca:fungicide = 0"

# EvS background
evS_evS_ctrl_alpha = "density + density:focs + density:bgs + density:focs:bgs = 0"
evS_evS_fung_alpha = "density + density:focs + density:bgs + density:focs:bgs + density:fungicide + density:focs:fungicide + density:bgs:fungicide + density:focs:bgs:fungicide = 0"
mv_evS_ctrl_alpha = "density +  density:bgs = 0"
mv_evS_fung_alpha = "density +  density:bgs + density:fungicide + density:bgs:fungicide = 0"
evA_evS_ctrl_alpha = "density +  density:bgs + density:foca + density:foca:bgs = 0"
evA_evS_fung_alpha = "density +  density:bgs + density:foca + density:foca:bgs + density:fungicide +  density:bgs:fungicide + density:foca:fungicide + density:foca:bgs:fungicide = 0"

# EvA intra vs. inter
evA_evA_ctrl_alpha = "density + density:foca + density:bga + density:foca:bga = 0"
evA_evA_fung_alpha = "density + density:foca + density:bga + density:foca:bga + density:fungicide + density:foca:fungicide + density:bga:fungicide + density:foca:bga:fungicide = 0"
mv_evA_ctrl_alpha = "density +  density:bga = 0"
mv_evA_fung_alpha = "density +  density:bga + density:fungicide + density:bga:fungicide = 0"
evS_evA_ctrl_alpha = "density + density:bga + density:focs + density:focs:bga = 0"
evS_evA_fung_alpha = "density + density:bga + density:focs + density:focs:bga + density:fungicide + density:bga:fungicide + density:focs:fungicide + density:focs:bga:fungicide = 0"

growthD1alphas <- hypothesis(growthD1Mod, 
                           c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                             evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                             evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                             evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                             mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                             evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                             evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                             mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                             evS_evA_ctrl_alpha, evS_evA_fung_alpha))

growthD2alphas <- hypothesis(growthD2Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

# combine alphas
alphaDat <- growthD1alphas[[1]] %>%
  mutate(year = "2018",
         foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                        "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                        "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                        "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                        "s_a_ctrl", "s_a_fung")) %>%
  full_join(growthD2alphas[[1]] %>%
              mutate(year = "2019",
                     foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                                    "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                                    "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                                    "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                                    "s_a_ctrl", "s_a_fung"))) %>%
  select(-Hypothesis) %>%
  rowwise() %>%
  mutate(foc = str_split(foc_bg_trt, "_")[[1]][1],
         bg = str_split(foc_bg_trt, "_")[[1]][2],
         trt = str_split(foc_bg_trt, "_")[[1]][3]) %>%
  ungroup() %>%
  mutate(treatment = fct_recode(trt, "control (water)" = "ctrl",
                          "fungicide" = "fung"),
         sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))

# edit to save
alphaDatSave <- alphaDat %>%
  left_join(predDat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)")) %>%
  select(year, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper)

# save
write_csv(alphaDatSave, "output/focal_growth_competition_coefficients_2018_2019_density_exp.csv")

# combine with preddat
predDat2 <- predDat %>%
  left_join(alphaDat)

# raw data
figDat <- growthD1Dat2 %>%
  select(site, plot, treatment, sp, age, ID, focal, background, density, plant_growth) %>%
  mutate(year = "2018") %>%
  full_join(growthD2Dat2 %>%
              select(site, plot, treatment, sp, age, ID, focal, background, density, plant_growth) %>%
              mutate(year = "2019")) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# alphas for figure
# sig beta values
alphaDat2 <- alphaDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, foc, focal, bg, background) %>%
              summarise(density = max(density))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(upper = max(CI.Upper)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(growth = max(plant_growth))) %>%
              rowwise() %>%
              mutate(plant_growth = max(c(growth, upper))) %>%
              ungroup() %>%
              select(-c(growth, upper))) %>%
  mutate(comp = round(Estimate, 2))


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
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.1, 0, 0.2, 2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

yearText <- tibble(year = c("2018", "2019"),
                   density = c(3.1, 3.1),
                   plant_growth = c(0.87, 3.09),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# 2018
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# 2019
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  geom_text(data = filter(alphaDat2, year == "2019"), 
            aes(label = paste("alpha", " == ", comp, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg <- get_legend(pairD1Fig)

# combine plots
combFig <- plot_grid(pairD1Fig + theme(legend.position = "none"), pairD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9),
                     label_x = c(0, -0.01))

# combine
pdf("output/focal_growth_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()


#### supplementary figure ####

yearText2 <- yearText %>%
  mutate(density = c(-0.1, -0.1),
         plant_growth = c(0.45, 2.3))

# 2018
rawD1Fig <- ggplot(growthD1Dat, aes(x = density, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  geom_text(data = filter(yearText2, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_color_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# 2019
rawD2Fig <- ggplot(growthD2Dat, aes(x = density, y = plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  geom_text(data = filter(yearText2, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg2 <- get_legend(rawD1Fig)

# combine plots
combFig2 <- plot_grid(rawD1Fig + theme(legend.position = "none"), rawD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9),
                     label_x = c(0, -0.01))

# combine
pdf("output/focal_raw_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig2, leg2,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()


##### info ####

# file: direct_disease_figure_density_exp
# author: Amy Kendig
# date last edited: 7/22/20
# goal: figure of direct disease model results


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(glmmTMB)

# import models
load("output/mv_survival_high_mv_model_2018_2019_density_exp.rda")
load("output/ev_seedling_survival_high_mv_model_2018_2019_density_exp.rda")
load("output/ev_adult_survival_high_mv_model_2018_2019_density_exp.rda")
load("output/mv_biomass_high_mv_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_seedling_biomass_high_mv_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_adult_biomass_high_mv_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/mv_seeds_high_mv_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_seedling_seeds_high_mv_model_2018_2019_density_exp.rda")
load("output/ev_adult_seeds_high_mv_model_2018_2019_density_exp.rda")
load("output/mv_germination_high_mv_model_2018_density_exp.rda")

# examine models
summary(mv_hi_mv_surv_mod)
summary(evs_hi_mv_surv_mod)
summary(eva_hi_mv_surv_mod)
summary(mv_hi_mv_bio_fung_mod)
summary(evs_hi_mv_bio_fung_mod)
summary(eva_hi_mv_bio_fung_mod)
summary(mv_hi_mv_seed_fung_mod)
summary(evs_hi_mv_seed_mod)
summary(eva_hi_mv_seed_mod)
summary(mv_hi_mv_germ_mod)


#### figure settings ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.title = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.box.margin = margin(-10, -10, -10, -10))

# colors
col_pal = c("#C0A76D", "#55A48B")


#### mean and intervals ####

mean_fun <- function(mod, logit, plant_type, response){
  
  out <- posterior_samples(mod) %>%
    as_tibble() %>%
    mutate(ctrl = case_when(logit == T ~ exp(b_Intercept) / (1 + exp(b_Intercept)),
                            logit == F ~ exp(b_Intercept)),
           fung = case_when(logit == T ~ exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide)),
                            logit == F ~ exp(b_Intercept + b_fungicide)), 
           fung_eff = log(fung / ctrl)) %>%
    select(fung_eff) %>%
    mean_hdci() %>%
    mutate(plant_type = plant_type,
           response = response)
  
  return(out)
}


#### extract mean and intervals ####

# mean function
mod_eff <- mean_fun(mv_hi_mv_surv_mod, T, "Mv", "survival") %>%
  rbind(mean_fun(evs_hi_mv_surv_mod, T, "Ev seedling", "survival")) %>%
  rbind(mean_fun(eva_hi_mv_surv_mod, T, "Ev adult", "survival")) %>%
  rbind(mean_fun(mv_hi_mv_bio_fung_mod, F, "Mv", "biomass")) %>%
  rbind(mean_fun(evs_hi_mv_bio_fung_mod, F, "Ev seedling", "biomass")) %>%
  rbind(mean_fun(eva_hi_mv_bio_fung_mod, F, "Ev adult", "biomass")) %>%
  rbind(mean_fun(mv_hi_mv_seed_fung_mod, F, "Mv", "seeds")) %>%
  rbind(mean_fun(evs_hi_mv_seed_mod, F, "Ev seedling", "seeds")) %>%
  rbind(mean_fun(eva_hi_mv_seed_mod, F, "Ev adult", "seeds")) %>%
  rbind(mean_fun(mv_hi_mv_germ_mod, T, "Mv", "germination"))

# modify dataset
mod_eff2 <- mod_eff %>%
  mutate(plant_type = fct_relevel(plant_type, "Mv", "Ev seedling"),
         Response = fct_relevel(response, "survival", "biomass", "seeds", "germination"))


#### fungicide effects figure ####

fung_fig <- ggplot(mod_eff2, aes(x = plant_type, y = fung_eff)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = Response), width = 0, size = 0.1, position = position_dodge(0.5)) +
  geom_point(aes(shape = Response), size = 2, position = position_dodge(0.5)) +
  ylab("Fungicide effect:log(fungicide / control)") +
  scale_x_discrete("Focal plant type", labels = expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult"))) +
  temp_theme

pdf("output/direct_disease_hi_mv_figure_density_exp.pdf",
    width = 5, height = 3)
fung_fig
dev.off()

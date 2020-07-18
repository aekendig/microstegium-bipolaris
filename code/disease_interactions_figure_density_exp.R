##### info ####

# file: disease_interactions_figure_density_exp
# author: Amy Kendig
# date last edited: 7/16/20
# goal: figure of disease-mediated interactions


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)


#### figure settings ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text.y = element_text(size = 6, color="black"),
        axis.text.x = element_text(size = 6, color="black", angle = 20, hjust = 1),
        axis.title = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.text.align = 0,
        legend.position = "bottom",
        legend.direction = "horizontal")

# colors
col_pal = c("#a6611a", "#018571")


#### survival ####

### import models ###
load("output/mv_survival_mv_background_model_2019_density_exp.rda")
load("output/mv_survival_ev_seedling_background_model_2019_density_exp.rda")
load("output/mv_survival_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_seedling_survival_mv_background_model_2019_density_exp.rda")
load("output/ev_seedling_survival_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_seedling_survival_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_adult_survival_mv_background_model_2019_density_exp.rda")
load("output/ev_adult_survival_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_adult_survival_ev_adult_background_model_2019_density_exp.rda")

# indicate variables with sig CI
summary(mv_mv_surv_mod) # background
summary(mv_evs_surv_mod) # background
summary(mv_eva_surv_mod)
summary(evs_mv_surv_mod) # fung:background
summary(evs_evs_surv_mod) # background
summary(evs_eva_surv_mod) # fung:background
summary(eva_mv_surv_mod)
summary(eva_evs_surv_mod) # background
summary(eva_eva_surv_mod)

### function for density effects ###
surv_dens_fun <- function(mod, focal_type, plant_type){
  
  out <- posterior_samples(mod) %>%
    as_tibble() %>%
    rename("b_fungicide_background_density" = "b_fungicide:background_density") %>%
    mutate(water = exp(b_Intercept + b_background_density) / (1 + exp(b_Intercept + b_background_density)) - exp(b_Intercept) / (1 + exp(b_Intercept)),
           fungicide = exp(b_Intercept + b_fungicide + b_background_density + b_fungicide_background_density) / (1 + exp(b_Intercept + b_fungicide + b_background_density + b_fungicide_background_density)) - exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide))) %>%
    select(water, fungicide) %>%
    pivot_longer(cols = c("water", "fungicide"), names_to = "treatment", values_to = "density_effect") %>%
    group_by(treatment) %>%
    mean_hdci() %>%
    mutate(plant_type = plant_type,
           focal_type = focal_type)
  
  return(out)
}
# Used hdci because mv_mv_surv, mv_evs_surv, and mv_eva_surv had two hdi's, the one that is not well represented with hdci didn't make much sense anyway (lower and upper don't include effect). All other models are the same with hdci or hdi.

### density effects table ###

surv_dens_eff <- surv_dens_fun(mv_mv_surv_mod, "Mv", "Mv") %>%
  rbind(surv_dens_fun(mv_evs_surv_mod, "Mv", "Evs")) %>%
  rbind(surv_dens_fun(mv_eva_surv_mod, "Mv", "Eva")) %>%
  rbind(surv_dens_fun(evs_mv_surv_mod, "Evs", "Mv")) %>%
  rbind(surv_dens_fun(evs_evs_surv_mod, "Evs", "Evs")) %>%
  rbind(surv_dens_fun(evs_eva_surv_mod, "Evs", "Eva")) %>%
  rbind(surv_dens_fun(eva_mv_surv_mod, "Eva", "Mv")) %>%
  rbind(surv_dens_fun(eva_evs_surv_mod, "Eva", "Evs")) %>%
  rbind(surv_dens_fun(eva_eva_surv_mod, "Eva", "Eva"))

# modify table
surv_dens_eff2 <- surv_dens_eff %>%
  mutate(trt_plant_type = paste(plant_type, treatment, sep = "_") %>%
           fct_relevel("Mv_water", "Mv_fungicide", "Evs_water", "Evs_fungicide", "Eva_water", "Eva_fungicide"),
         focal_type = fct_relevel(focal_type, "Mv", "Evs", "Eva"),
         plant_type = fct_relevel(plant_type, "Mv", "Evs", "Eva"),
         Treatment = recode(treatment, "water" = "control (water)"),
         fung_eff = case_when(focal_type == "Evs" & plant_type == "Mv" & treatment == "fungicide" ~ "*    ",
                              focal_type == "Evs" & plant_type == "Eva" & treatment == "fungicide"  ~ "*    ",
                              TRUE ~ ""),
         line_min = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ -Inf,
                              TRUE ~ 0),
         line_max = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ Inf,
                              TRUE ~ 0))


### figure ###

surv_fig <- ggplot(surv_dens_eff2, aes(focal_type, density_effect)) +
  geom_errorbar(aes(ymin = line_min, ymax = line_max, group = trt_plant_type), size = 0.2, width = 0, position = position_dodge(1.2), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = trt_plant_type), width = 0, size = 0.1, position = position_dodge(0.8)) +
  geom_point(aes(shape = plant_type, fill = Treatment, group = trt_plant_type), size = 1.5, position = position_dodge(0.8)) +
  geom_text(aes(label = fung_eff, y = (.upper + 0.02), group = trt_plant_type), position = position_dodge(0.8)) +
  scale_x_discrete("Focal plant type", labels = expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult"))) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 24, 22), 
                     labels =expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult")),
                     name = "Background plant type") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Focal plant type") +
  ylab("Change in proportion surviving with one plant") +
  temp_theme
  

#### biomass ####

### import models ###
load("output/mv_biomass_mv_background_model_2019_density_exp.rda")
load("output/mv_biomass_ev_seedling_background_model_2019_density_exp.rda")
load("output/mv_biomass_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_seedling_biomass_mv_background_model_2019_density_exp.rda")
load("output/ev_seedling_biomass_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_seedling_biomass_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_adult_biomass_mv_background_model_2019_density_exp.rda")
load("output/ev_adult_biomass_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_adult_biomass_ev_adult_background_model_2019_density_exp.rda")

summary(mv_mv_bio_mod)

### test for treatment effect ###
coef_fun <- function(mod){
  
  out <- posterior_samples(mod) %>%
    as_tibble() %>%
    mutate(diff = b_alpha_Treatmentfungicide - b_alpha_Treatmentcontrolwater) %>%
    select(diff) %>%
    mean_hdi()
  
  return(out)
}

coef_fun(mv_mv_bio_mod)
coef_fun(mv_evs_bio_mod)
coef_fun(mv_eva_bio_mod)
coef_fun(evs_mv_bio_mod)
coef_fun(evs_evs_bio_mod)
coef_fun(evs_eva_bio_mod)
coef_fun(eva_mv_bio_mod)
coef_fun(eva_evs_bio_mod)
coef_fun(eva_eva_bio_mod)
# none

### extract coefficients ###
bio_dens_eff <- fixef(mv_mv_bio_mod)[c(3, 4), ] %>%
  as_tibble() %>%
  mutate(treatment = c("water", "fungicide"),
         focal_type = "Mv",
         plant_type = "Mv") %>%
  rbind(fixef(mv_evs_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Mv",
                 plant_type = "Evs")) %>%
  rbind(fixef(mv_eva_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Mv",
                 plant_type = "Eva")) %>%
  rbind(fixef(evs_mv_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Mv")) %>%
  rbind(fixef(evs_evs_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Evs")) %>%
  rbind(fixef(evs_eva_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Eva")) %>%
  rbind(fixef(eva_mv_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Mv")) %>%
  rbind(fixef(eva_evs_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Evs")) %>%
  rbind(fixef(eva_eva_bio_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Eva")) 

# modify table
bio_dens_eff2 <- bio_dens_eff %>%
  mutate(trt_plant_type = paste(plant_type, treatment, sep = "_") %>%
           fct_relevel("Mv_water", "Mv_fungicide", "Evs_water", "Evs_fungicide", "Eva_water", "Eva_fungicide"),
         focal_type = fct_relevel(focal_type, "Mv", "Evs", "Eva"),
         plant_type = fct_relevel(plant_type, "Mv", "Evs", "Eva"),
         Treatment = recode(treatment, "water" = "control (water)"),
         line_min = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ -Inf,
                              TRUE ~ 0),
         line_max = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ Inf,
                              TRUE ~ 0))

### figure ###

bio_fig <- ggplot(bio_dens_eff2, aes(focal_type, Estimate)) +
  geom_errorbar(aes(ymin = line_min, ymax = line_max, group = trt_plant_type), size = 0.2, width = 0, position = position_dodge(1.2), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5, group = trt_plant_type), width = 0, size = 0.1, position = position_dodge(0.8)) +
  geom_point(aes(shape = plant_type, fill = Treatment, group = trt_plant_type), size = 1.5, position = position_dodge(0.8)) +
  scale_x_discrete("Focal plant type", labels = expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult"))) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 24, 22), 
                     labels =expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult")),
                     name = "Background plant type") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Focal plant type") +
  ylab("Competitive effect on biomass (g)") +
  temp_theme +
  theme(legend.position = "none")


#### seeds ####

### import models ###
load("output/mv_seeds_mv_background_model_2019_density_exp.rda")
load("output/mv_seeds_ev_seedling_background_model_2019_density_exp.rda")
load("output/mv_seeds_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_seedling_seeds_mv_background_model_2019_density_exp.rda")
load("output/ev_seedling_seeds_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_seedling_seeds_ev_adult_background_model_2019_density_exp.rda")
load("output/ev_adult_seeds_mv_background_model_2019_density_exp.rda")
load("output/ev_adult_seeds_ev_seedling_background_model_2019_density_exp.rda")
load("output/ev_adult_seeds_ev_adult_background_model_2019_density_exp.rda")

summary(mv_mv_seed_mod)

### test for tretment effect ###
coef_fun(mv_mv_bio_mod)
coef_fun(mv_evs_bio_mod)
coef_fun(mv_eva_bio_mod)
coef_fun(evs_mv_bio_mod)
coef_fun(evs_evs_bio_mod)
coef_fun(evs_eva_bio_mod)
coef_fun(eva_mv_bio_mod)
coef_fun(eva_evs_bio_mod)
coef_fun(eva_eva_bio_mod)
# none

### extract coefficients ###
seed_dens_eff <- fixef(mv_mv_seed_mod)[c(3, 4), ] %>%
  as_tibble() %>%
  mutate(treatment = c("water", "fungicide"),
         focal_type = "Mv",
         plant_type = "Mv") %>%
  rbind(fixef(mv_evs_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Mv",
                 plant_type = "Evs")) %>%
  rbind(fixef(mv_eva_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Mv",
                 plant_type = "Eva")) %>%
  rbind(fixef(evs_mv_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Mv")) %>%
  rbind(fixef(evs_evs_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Evs")) %>%
  rbind(fixef(evs_eva_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Evs",
                 plant_type = "Eva")) %>%
  rbind(fixef(eva_mv_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Mv")) %>%
  rbind(fixef(eva_evs_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Evs")) %>%
  rbind(fixef(eva_eva_seed_mod)[c(3, 4), ] %>%
          as_tibble() %>%
          mutate(treatment = c("water", "fungicide"),
                 focal_type = "Eva",
                 plant_type = "Eva")) 

# modify table
seed_dens_eff2 <- seed_dens_eff %>%
  mutate(trt_plant_type = paste(plant_type, treatment, sep = "_") %>%
           fct_relevel("Mv_water", "Mv_fungicide", "Evs_water", "Evs_fungicide", "Eva_water", "Eva_fungicide"),
         focal_type = fct_relevel(focal_type, "Mv", "Evs", "Eva"),
         plant_type = fct_relevel(plant_type, "Mv", "Evs", "Eva"),
         Treatment = recode(treatment, "water" = "control (water)"),
         line_min = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ -Inf,
                              TRUE ~ 0),
         line_max = case_when(plant_type == "Eva" & treatment == "fungicide" & focal_type != "Eva" ~ Inf,
                              TRUE ~ 0))

### figure ###

seed_fig <- ggplot(seed_dens_eff2, aes(focal_type, Estimate)) +
  geom_errorbar(aes(ymin = line_min, ymax = line_max, group = trt_plant_type), size = 0.2, width = 0, position = position_dodge(1.2), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5, group = trt_plant_type), width = 0, size = 0.1, position = position_dodge(0.8)) +
  geom_point(aes(shape = plant_type, fill = Treatment, group = trt_plant_type), size = 1.5, position = position_dodge(0.8)) +
  scale_x_discrete("Focal plant type", labels = expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult"))) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 24, 22), 
                     labels =expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult")),
                     name = "Background plant type") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Focal plant type") +
  ylab("Competitive effect on seeds") +
  temp_theme +
  theme(legend.position = "none")


#### combine figures ###

# extract legend
leg = get_legend(surv_fig)

# Main figures
top_fig <- plot_grid(surv_fig + theme(legend.position = "none"),
                     bio_fig,
                     seed_fig,
                     nrow = 1,
                     labels = c("A", "B", "C"),
                     label_size = 10)

# with legend

pdf("output/disease_interactions_figure_density_exp.pdf",
    width = 7, height = 3.5)
plot_grid(top_fig, leg,
          nrow = 2,
          rel_heights = c(1, 0.07))
dev.off()
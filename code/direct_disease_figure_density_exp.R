##### info ####

# file: direct_disease_figure_density_exp
# author: Amy Kendig
# date last edited: 7/15/20
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
load("output/mv_survival_no_background_model_2018_2019_density_exp.rda")
load("output/ev_seedling_survival_no_background_model_2018_2019_density_exp.rda")
load("output/ev_adult_survival_no_background_model_2018_2019_density_exp.rda")
load("output/mv_biomass_no_background_model_2019_density_exp.rda")
load("output/mv_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_seedling_biomass_no_background_model_2019_density_exp.rda")
load("output/ev_seedling_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_adult_biomass_no_background_model_2019_density_exp.rda")
load("output/ev_adult_biomass_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/mv_seeds_no_background_model_2019_density_exp.rda")
load("output/mv_seeds_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
load("output/ev_seedling_seeds_no_background_model_2018_2019_density_exp.rda")
load("output/ev_adult_seeds_no_background_model_2018_2019_density_exp.rda")
load("output/mv_germination_no_background_model_2018_density_exp.rda")
# load("output/ev_seedling_survival_severity_model_jul_2019_density_exp.rda")
# load("output/ev_seedling_survival_severity_model_early_aug_2019_density_exp.rda")
load("output/ev_seedling_biomass_severity_model_jul_2019_density_exp.rda")

# examine models
summary(mv_no_back_surv_mod)
summary(evs_no_back_surv_mod)
summary(eva_no_back_surv_mod)
# pp_check(eva_no_back_surv_mod, nsamples = 100)
summary(mv_no_back_bio_mod)
summary(mv_no_back_bio_fung_mod)
summary(evs_no_back_bio_mod)
summary(evs_no_back_bio_fung_mod)
summary(eva_no_back_bio_mod)
summary(eva_no_back_bio_fung_mod)
summary(mv_no_back_seed_mod)
summary(mv_no_back_seed_fung_mod)
summary(evs_no_back_seed_mod)
summary(eva_no_back_seed_mod)
summary(mv_no_back_germ_mod)
# summary(evs_surv_sev_jul_19_mod)
# summary(eva_surv_sev_eau_19_mod)
summary(evs_bio_sev_jul_19_mod)

# import data
# evs_surv_sev_dat <- read_csv("intermediate-data/ev_seedling_survival_severity_data_2019_dens_exp.csv")
# eva_surv_sev_dat <- read_csv("intermediate-data/ev_adult_survival_severity_data_2019_dens_exp.csv")
evs_bio_sev_dat <- read_csv("intermediate-data/ev_seedling_biomass_severity_data_2019_dens_exp.csv")


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

# Used hdci because Ev adult survival had two hdi's, the one that is not well represented with hdci didn't make much sense anyway (lower and upper don't include effect). All other models are the same with hdci or hdi.


#### extract mean and intervals ####

# mean function
mod_eff <- mean_fun(mv_no_back_surv_mod, T, "Mv", "survival") %>%
  rbind(mean_fun(evs_no_back_surv_mod, T, "Ev seedling", "survival")) %>%
  rbind(mean_fun(eva_no_back_surv_mod, T, "Ev adult", "survival")) %>%
  rbind(mean_fun(mv_no_back_bio_fung_mod, F, "Mv", "biomass")) %>%
  rbind(mean_fun(evs_no_back_bio_fung_mod, F, "Ev seedling", "biomass")) %>%
  rbind(mean_fun(eva_no_back_bio_fung_mod, F, "Ev adult", "biomass")) %>%
  rbind(mean_fun(mv_no_back_seed_fung_mod, F, "Mv", "seeds")) %>%
  rbind(mean_fun(evs_no_back_seed_mod, F, "Ev seedling", "seeds")) %>%
  rbind(mean_fun(eva_no_back_seed_mod, F, "Ev adult", "seeds")) %>%
  rbind(mean_fun(mv_no_back_germ_mod, T, "Mv", "germination"))

# negative numbers: all three biomass, Ev seedling seeds

# modify dataset
mod_eff2 <- mod_eff %>%
  mutate(plant_type = fct_relevel(plant_type, "Mv", "Ev seedling"),
         Response = fct_relevel(response, "survival", "biomass", "seeds", "germination"))


#### fungicide effects figure ####

fung_fig <- ggplot(mod_eff2, aes(x = plant_type, y = fung_eff)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = Response), width = 0, size = 0.1, position = position_dodge(0.5)) +
  geom_point(aes(shape = Response), size = 2, position = position_dodge(0.5)) +
  ylab("Fungicide effect: log(fungicide / control)") +
  scale_x_discrete("Focal plant type", labels = expression(italic(Microstegium), paste(italic(Elymus), " seedling", sep = ""), paste(italic(Elymus), " adult"))) +
  temp_theme


#### fungicide effects figure for presentation ####

pdf("output/direct_disease_presentation_figure_density_exp.pdf",
    width = 4.5, height = 4)
ggplot(mod_eff2, aes(x = plant_type, y = fung_eff)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper, group = Response), width = 0, size = 0.1, position = position_dodge(0.5)) +
  geom_point(aes(shape = Response, fill = Response), size = 3, position = position_dodge(0.5)) +
  ylab("Fungicide effect: log(fungicide / control)") +
  scale_x_discrete("Focal plant type", labels = expression(atop(atop("Invasive", italic(Microstegium)), NA), atop(atop("Native", paste(italic(Elymus), " seedling", sep = "")), NA), atop(atop("Native", paste(italic(Elymus), " adult", sep = "")), NA))) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21:24)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, color="black"),
        axis.text.x = element_text(size = 12, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.position = "bottom",
        legend.direction = "horizontal")
dev.off()


# #### Ev seedling survival ####
# 
# # min and max values
# evs_surv_range <- evs_surv_sev_dat %>%
#   summarise(min = min(severity_jul, na.rm = T),
#             max = 1)
# 
# # prediction data
# evs_surv_sim_dat <- tibble(severity_jul = seq(evs_surv_range$min, evs_surv_range$max, length.out = 200)) %>%
#   mutate(fungicide = 0,
#          site = NA,
#          plot = NA,
#          Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
#   mutate(survival = predict(evs_surv_sev_jul_19_mod, newdata = ., type = "response", re.form = NA),
#          survival_se = predict(evs_surv_sev_jul_19_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)
# 
# # plot
# ggplot(evs_surv_sim_dat, aes(severity_jul, survival)) +
#   geom_ribbon(aes(ymin = survival - survival_se, ymax = survival + survival_se, fill = Treatment), alpha = 0.4) +
#   geom_line(aes(color = Treatment)) +
#   annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
#   scale_color_manual(values = col_pal) +
#   scale_fill_manual(values = col_pal) +
#   ylab(expression(paste(italic(Elymus), " seedling survival", sep = ""))) +
#   xlab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
#   temp_theme
# # some of these values don't make sense, the survival is 1 across most of the informed range
# 
# 
# #### Ev adult survival ####
# 
# # min and max values
# eva_surv_range <- eva_surv_sev_dat %>%
#   summarise(min = min(severity_early_aug, na.rm = T),
#             max = max(severity_early_aug, na.rm = T))
# 
# # prediction data
# eva_surv_sim_dat <- tibble(severity_early_aug = seq(eva_surv_range$min, eva_surv_range$max, length.out = 200)) %>%
#   mutate(fungicide = 0,
#          site_plot = NA,
#          Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
#   mutate(survival = predict(eva_surv_sev_eau_19_mod, newdata = ., type = "response", re.form = NA),
#          survival_se = predict(eva_surv_sev_eau_19_mod, newdata = ., type = "response", re.form = NA, se.fit = T)$se.fit)
# 
# # plot
# ggplot(eva_surv_sim_dat, aes(severity_early_aug, survival)) +
#   geom_ribbon(aes(ymin = survival - survival_se, ymax = survival + survival_se, fill = Treatment), alpha = 0.4) +
#   geom_line(aes(color = Treatment)) +
#   annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
#   scale_color_manual(values = col_pal) +
#   scale_fill_manual(values = col_pal) +
#   ylab(expression(paste(italic(Elymus), " adult survival", sep = ""))) +
#   xlab(expression(paste(italic(Elymus), " lesions", sep = ""))) +
#   temp_theme
# # not useful - most of the action is at the very high end of the spectrum where values are uncertain and the effect size is tiny


#### Ev seedling biomass ####

# min and max values
evs_bio_range <- evs_bio_sev_dat %>%
  summarise(min = min(severity_jul, na.rm = T),
            max = 1)

# prediction data
evs_bio_sim_dat <- tibble(severity_jul = seq(evs_bio_range$min, evs_bio_range$max, length.out = 200)) %>%
  mutate(fungicide = 0,
         site = NA,
         plot = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(veg_weight.g = predict(evs_bio_sev_jul_19_mod, newdata = ., re.form = NA),
         veg_weight_se = predict(evs_bio_sev_jul_19_mod, newdata = ., re.form = NA, se.fit = T)$se.fit,
         lesions = severity_jul * 100)

# plot
ev_bio_fig <- ggplot(evs_bio_sim_dat, aes(lesions, veg_weight.g)) +
  geom_ribbon(aes(ymin = veg_weight.g - veg_weight_se, ymax = veg_weight.g + veg_weight_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste(italic(Elymus), " seedling biomass (g)  ", sep = ""))) +
  xlab("Lesions (% leaf surface)") +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  temp_theme +
  theme(legend.position = "none")


#### Ev seedling biomass figure for presentation ####

pdf("output/ev_biomass_severity_figure_2019_density_exp.pdf", width = 3.5, height = 3.5)
ggplot(evs_bio_sim_dat, aes(lesions, veg_weight.g)) +
  geom_point(data = evs_bio_sev_dat, alpha = 0.3, aes(x = severity_jul*100, color = Treatment)) +
  geom_ribbon(aes(ymin = veg_weight.g - veg_weight_se, ymax = veg_weight.g + veg_weight_se, fill = Treatment), alpha = 0.4) +
  geom_line(aes(color = Treatment)) +
  #annotate(geom = "text", label = "July", x = Inf, y = Inf, hjust = 1.2, vjust = 1.2, size = 2) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ylab(expression(paste("Native ", italic(Elymus), " seedling biomass (g)", sep = ""))) +
  xlab("% leaf area diseased") +
  scale_x_continuous(breaks = c(0, 50, 100))  +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.position = "none")
dev.off()

#### combine figures ####

# legend
leg = get_legend(fung_fig + theme(legend.justification = c(0.1, 0.5)))

# right plots
right_plots = plot_grid(ev_bio_fig, leg, 
                        nrow = 2,
                        rel_heights = c(1, 0.6))

pdf("output/direct_disease_figure_density_exp.pdf",
    width = 4.3, height = 3)
plot_grid(fung_fig + theme(legend.position = "none"), right_plots,
          nrow = 1,
          rel_widths = c(1, 0.6),
          labels = c("A", "B"),
          label_size = 10)
dev.off()

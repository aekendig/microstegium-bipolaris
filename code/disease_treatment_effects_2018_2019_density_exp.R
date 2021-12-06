##### info ####

# file: disease_treatment_effects_2018_2019_density_exp.R
# author: Amy Kendig
# date last edited: 12/6/21
# goal: treatment effects for Mv and Ev


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(-0.1, 0, -0.1, 0, unit = "cm"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 7, hjust = 0.5))

col_pal = c("black", "#238A8DFF")


#### severity ####

load("output/severity_fungicide_model_2019_density_exp.rda")

sev <- as_draws_df(trtD2Mod) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(a_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            a_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide)),
            s_W = exp(b_Intercept + b_plant_typeEv_seedling)/(1 +  exp(b_Intercept + b_plant_typeEv_seedling)),
            s_F = exp(b_Intercept + b_plant_typeEv_seedling + b_fungicide + b_plant_typeEv_seedling_fungicide)/(1 +  exp(b_Intercept + b_plant_typeEv_seedling + b_fungicide + b_plant_typeEv_seedling_fungicide)),
            m_W = exp(b_Intercept + b_plant_typeMv_seedling)/(1 +  exp(b_Intercept + b_plant_typeMv_seedling)),
            m_F = exp(b_Intercept + b_plant_typeMv_seedling + b_fungicide + b_plant_typeMv_seedling_fungicide)/(1 +  exp(b_Intercept + b_plant_typeMv_seedling + b_fungicide + b_plant_typeMv_seedling_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = parameter, 
            Species = c(rep("Competitor", 2), rep("Invader", 2), rep("Competitor", 2)) %>%
              fct_relevel("Invader"),
            Age = c(rep("adult", 2), rep("first-year", 4)) %>%
                      fct_relevel("first-year"),
            Treatment = if_else(str_ends(parameter, "_W") == T, "control", "fungicide"),
            Estimate = estimate * 100,
            Lower = .lower * 100,
            Upper = .upper * 100)

sev_fig <- ggplot(sev, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.3)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.3)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "Disease severity (%)") +
  fig_theme +
  theme(axis.text.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


#### germination ####

load("output/mv_germination_fungicide_model_2018_density_exp.rda")

mv_germ <- as_draws_df(mvGermD1Mod3)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Species = "Invader",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper)


load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

ev_germ <- as_draws_df(evGermMod2)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Species = "Competitor",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper)

germ <- mv_germ %>%
  full_join(ev_germ) %>%
  mutate(Species = fct_relevel(Species, "Invader"))

germ_fig <- ggplot(germ, aes(x = Species, y = Estimate)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.3)) +
  geom_point(aes(fill = Treatment), size = 2, stroke = 0.3, shape = 21, 
             position = position_dodge(0.3)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(y = "Germination fraction") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())


#### establishment ####

# import model
load("output/survival_fungicide_model_2019_density_exp.rda")

est <- as_draws_df(survFungD2Mod)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(e_A_W = exp(b_Intercept)/(1 + exp(b_Intercept)),
            e_A_F = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide)),
            e_P_W = exp(b_Intercept + b_focs)/(1 + exp(b_Intercept + b_focs)),
            e_P_F = exp(b_Intercept + b_focs + b_fungicide)/(1 + exp(b_Intercept + b_focs + b_fungicide))) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Parameter = parameter,
            Species = rep(c("Invader", "Competitor"), each = 2),
            Treatment = if_else(str_detect(parameter, "W") == T, "control", "fungicide"),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader"))

est_fig <- ggplot(est, aes(x = Species, y = Estimate)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.3)) +
  geom_point(aes(fill = Treatment), size = 2, stroke = 0.3, shape = 21, 
             position = position_dodge(0.3)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(y = "Establishment fraction") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())


#### biomass ####

# load model
load("output/focal_growth_density_model_2019_density_exp.rda")

# sample model
bio <- as_draws_df(growthD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(r_A = Intercept,
            r_F = (Intercept + focs),
            r_P = (Intercept + foca),
            r_A_fung = (Intercept + fungicide),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide)) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Parameter = parameter,
            Species = rep(c("Invader",
                            "Competitor",
                            "Competitor"), each = 2),
            Age = c(rep("first-year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader"),
         Age = fct_relevel(Age, "first-year"))

bio_fig <- ggplot(bio, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.3)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.3)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "ln(Biomass) (g)") +
  fig_theme +
  theme(legend.position = "none")


#### seeds ####

# load model
load("output/focal_seed_density_model_2019_density_exp.rda")

# sample model
seed <- as_draws_df(seedD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(r_A = Intercept,
            r_F = (Intercept + focs),
            r_P = (Intercept + foca),
            r_A_fung = (Intercept + fungicide),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide)) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Parameter = parameter,
            Species = rep(c("Invader",
                            "Competitor",
                            "Competitor"), each = 2),
            Age = c(rep("first-year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader"),
         Age = fct_relevel(Age, "first-year"))

seed_fig <- ggplot(seed, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.3)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.3)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "ln(Seeds + 1)") +
  fig_theme +
  theme(legend.position = "none")



#### combine ####

# legend
leg <- get_legend(sev_fig)
top_row <- plot_grid(sev_fig + theme(legend.position = "none"), germ_fig, 
                     align = "h", axis = "l", ncol = 2,
                     labels = LETTERS[1:2],
                     label_size = 10)
mid_row <- plot_grid(est_fig, bio_fig,
                     align = "h", axis = "l", ncol = 2,
                     labels = LETTERS[3:4],
                     label_size = 10)
bottom_row <- plot_grid(seed_fig, leg, 
                        align = "h", axis = "l", ncol = 2,
                        labels = c(LETTERS[5], ""),
                        label_size = 10)

pdf("output/disease_treatment_effects_2018_2019_density_exp.pdf", width = 3.54, height = 4.33)
# plot_grid(top_row, mid_row, bottom_row, ncol = 1) # lines up middle plots
plot_grid(sev_fig + theme(legend.position = "none"), germ_fig,
          est_fig, bio_fig,
          ncol = 2,
          seed_fig, leg,
          labels = c(LETTERS[1:5], ""),
          label_size = 10)
dev.off()
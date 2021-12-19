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
library(lazerhawk)

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
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.box = "vertical",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 7, hjust = 0.5))

col_pal = c("black", "#238A8DFF")


#### severity ####

# load model
load("output/severity_fungicide_model_2019_density_exp.rda")

# draws
sev <- as_draws_df(trtD2Mod) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(a_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            a_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide)),
            s_W = exp(b_Intercept + b_plant_typeEv_seedling)/(1 +  exp(b_Intercept + b_plant_typeEv_seedling)),
            s_F = exp(b_Intercept + b_plant_typeEv_seedling + b_fungicide + b_plant_typeEv_seedling_fungicide)/(1 +  exp(b_Intercept + b_plant_typeEv_seedling + b_fungicide + b_plant_typeEv_seedling_fungicide)),
            m_W = exp(b_Intercept + b_plant_typeMv_seedling)/(1 +  exp(b_Intercept + b_plant_typeMv_seedling)),
            m_F = exp(b_Intercept + b_plant_typeMv_seedling + b_fungicide + b_plant_typeMv_seedling_fungicide)/(1 +  exp(b_Intercept + b_plant_typeMv_seedling + b_fungicide + b_plant_typeMv_seedling_fungicide))) 

# summarize for fig
sev2 <- sev %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = parameter, 
            Species = c(rep("Competitor", 2), rep("Invader", 2), rep("Competitor", 2)) %>%
              fct_relevel("Invader"),
            Age = c(rep("adult", 2), rep("1st year", 4)) %>%
                      fct_relevel("1st year"),
            Treatment = if_else(str_ends(parameter, "_W") == T, "control", "fungicide"),
            Estimate = estimate * 100,
            Lower = .lower * 100,
            Upper = .upper * 100)

# fig
sev_fig <- ggplot(sev2, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "Disease severity (%)") +
  fig_theme +
  theme(axis.text.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

# values for text
sev %>%
  transmute(mv = 100 * (m_F - m_W) / m_W,
            evS = 100 * (s_F - s_W) / s_W,
            evA = 100 * (a_F - a_W) / a_W) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate)

# table
write_csv(brms_SummaryTable(trtD2Mod), "output/severity_fungicide_model_2019_density_exp.csv")


#### germination ####

# load model
load("output/mv_germination_fungicide_model_2018_density_exp.rda")

# draws and summarize
mv_germ <- as_draws_df(mvGermD1Mod3)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) 

mv_germ2 <- mv_germ %>%
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

# load model
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# draws and summarize
ev_germ <- as_draws_df(evGermMod2)  %>%
  transmute(g_W = exp(b_Intercept)/(1 +  exp(b_Intercept)),
            g_F = exp(b_Intercept + b_fungicide)/(1 +  exp(b_Intercept + b_fungicide))) 

ev_germ2 <- ev_germ %>%
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

# combine
germ <- mv_germ2 %>%
  full_join(ev_germ2) %>%
  mutate(Species = fct_relevel(Species, "Invader"))

# fig
germ_fig <- ggplot(germ, aes(x = Species, y = Estimate)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment), size = 2, stroke = 0.3, shape = 21, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(y = "Germination fraction") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# values for text
mv_germ %>%
  transmute(eff = (g_F - g_W) / g_W) %>%
  mean_hdci(eff)

ev_germ %>%
  transmute(eff = (g_F - g_W) / g_W) %>%
  mean_hdci(eff)

# tables
write_csv(brms_SummaryTable(mvGermD1Mod3), "output/mv_germination_fungicide_model_2018_density_exp.csv")
write_csv(brms_SummaryTable(evGermMod2), "output/ev_germination_fungicide_model_2018_2019_density_exp.csv")


#### establishment ####

# load model
load("output/survival_fungicide_model_2019_density_exp.rda")

# draws
est <- as_draws_df(survFungD2Mod)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(e_A_W = exp(b_Intercept)/(1 + exp(b_Intercept)),
            e_A_F = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide)),
            e_P_W = exp(b_Intercept + b_focs)/(1 + exp(b_Intercept + b_focs)),
            e_P_F = exp(b_Intercept + b_focs + b_fungicide)/(1 + exp(b_Intercept + b_focs + b_fungicide)))

# summarize for fig
est2 <- est %>%
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

# fig
est_fig <- ggplot(est2, aes(x = Species, y = Estimate)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment), size = 2, stroke = 0.3, shape = 21, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(y = "Establishment fraction") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# values for text
est %>%
  transmute(mv = 100 * (e_A_F - e_A_W) / e_A_W,
            ev = 100 * (e_P_F - e_P_W) / e_P_W) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate)

# table
write_csv(brms_SummaryTable(survFungD2Mod), 
          "output/survival_fungicide_model_2019_density_exp.csv")


#### biomass ####

# load model
load("output/focal_growth_density_model_2019_density_exp.rda")

# draws
bio <- as_draws_df(growthD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(r_A = Intercept,
            r_F = (Intercept + focs),
            r_P = (Intercept + foca),
            r_A_fung = (Intercept + fungicide),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide))

# summarize for figure
bio2 <- bio %>%
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
            Age = c(rep("1st year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader"),
         Age = fct_relevel(Age, "1st year"))

# fig
bio_fig <- ggplot(bio2, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "ln(Biomass) (g)") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# values for text
bio %>%
  transmute(mv = 100 * (r_A_fung - r_A) / r_A,
            evS = 100 * (r_F_fung - r_F) / r_F,
            evA = 100 * (r_P_fung - r_P) / r_P) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate)

# table
write_csv(brms_SummaryTable(growthD2Mod), 
          "output/focal_growth_density_model_2019_density_exp.csv")


#### seeds ####

# load model
load("output/focal_seed_density_model_2019_density_exp.rda")

# draws
seed <- as_draws_df(seedD2Mod)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(r_A = Intercept,
            r_F = (Intercept + focs),
            r_P = (Intercept + foca),
            r_A_fung = (Intercept + fungicide),
            r_F_fung = (Intercept + focs + fungicide + focs_fungicide),
            r_P_fung = (Intercept + foca + fungicide + foca_fungicide)) 

# summarize for fig
seed2 <- seed %>%
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
            Age = c(rep("1st year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader"),
         Age = fct_relevel(Age, "1st year"))

seed_fig <- ggplot(seed2, aes(x = Species, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "ln(Seeds + 1)") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# values for text
seed %>%
  transmute(mv = 100 * (r_A_fung - r_A) / r_A,
            evS = 100 * (r_F_fung - r_F) / r_F,
            evA = 100 * (r_P_fung - r_P) / r_P) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate)

# table
write_csv(brms_SummaryTable(seedD2Mod), 
          "output/focal_seed_density_model_2019_density_exp.csv")


#### biomass comp. coeff. ####

# load model
load("output/focal_growth_biomass_model_2019_density_exp.rda")

# draws
alpha <- as_draws_df(growthD2Mod2)  %>%
  rename_with(str_replace_all, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = "b_", replacement = "") %>%
  transmute(alpha_AA = plot_biomass,
            alpha_FA = plot_biomass + focs_plot_biomass,
            alpha_PA = plot_biomass + foca_plot_biomass,
            alpha_AF = plot_biomass + plot_biomass_bgs,
            alpha_FF = plot_biomass + focs_plot_biomass + plot_biomass_bgs + focs_plot_biomass_bgs,
            alpha_PF = plot_biomass + foca_plot_biomass + plot_biomass_bgs + foca_plot_biomass_bgs,
            alpha_AP = plot_biomass + plot_biomass_bga,
            alpha_FP = plot_biomass + plot_biomass_bga + focs_plot_biomass + focs_plot_biomass_bga,
            alpha_PP = plot_biomass + foca_plot_biomass + plot_biomass_bga + foca_plot_biomass_bga,
            alpha_AA_fung = plot_biomass + fungicide_plot_biomass,
            alpha_FA_fung = plot_biomass + fungicide_plot_biomass + focs_plot_biomass + focs_fungicide_plot_biomass,
            alpha_PA_fung = plot_biomass + fungicide_plot_biomass + foca_plot_biomass + foca_fungicide_plot_biomass,
            alpha_AF_fung = plot_biomass + plot_biomass_bgs + fungicide_plot_biomass + fungicide_plot_biomass_bgs,
            alpha_FF_fung = plot_biomass + focs_plot_biomass + plot_biomass_bgs + focs_plot_biomass_bgs + fungicide_plot_biomass + focs_fungicide_plot_biomass + fungicide_plot_biomass_bgs + focs_fungicide_plot_biomass_bgs,
            alpha_PF_fung = plot_biomass + foca_plot_biomass + plot_biomass_bgs + foca_plot_biomass_bgs + fungicide_plot_biomass + foca_fungicide_plot_biomass + fungicide_plot_biomass_bgs + foca_fungicide_plot_biomass_bgs,
            alpha_AP_fung = plot_biomass + plot_biomass_bga + fungicide_plot_biomass + fungicide_plot_biomass_bga,
            alpha_FP_fung = plot_biomass + plot_biomass_bga + focs_plot_biomass + focs_plot_biomass_bga + fungicide_plot_biomass + fungicide_plot_biomass_bga + focs_fungicide_plot_biomass + focs_fungicide_plot_biomass_bga,
            alpha_PP_fung = plot_biomass + foca_plot_biomass + plot_biomass_bga + foca_plot_biomass_bga + fungicide_plot_biomass + foca_fungicide_plot_biomass + fungicide_plot_biomass_bga + foca_fungicide_plot_biomass_bga) 

# summarize for fig
alpha2 <- alpha %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  mutate(estimate = -1 * estimate) %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  ungroup() %>%
  transmute(Treatment = if_else(str_detect(parameter, "_fung") == T, "fungicide", "control"),
            Parameter = str_replace(parameter, "_fung", ""),
            Focal = if_else(str_starts(parameter, "alpha_A"), "Invader", "Competitor") %>%
              fct_relevel("Invader"),
            Age = if_else(str_starts(parameter, "alpha_P"), "adult", "1st year") %>%
              fct_relevel("1st year"),
            Competitor = case_when(str_ends(Parameter, "A") ~ "Invader",
                                   str_ends(Parameter, "F") ~ "1st year comp.",
                                   str_ends(Parameter, "P") ~ "Adult comp."),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper)

# separate by competitor
alpha_A <- alpha2 %>% filter(Competitor == "Invader")
alpha_F <- alpha2 %>% filter(Competitor == "1st year comp.")
alpha_P <- alpha2 %>% filter(Competitor == "Adult comp.")

# figs
alphaA_fig <- ggplot(alpha_A, aes(x = Focal, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "Invader competition") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

alphaF_fig <- ggplot(alpha_F, aes(x = Focal, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "1st yr comp. competition") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 7, color = "black", hjust = 0))

alphaP_fig <- ggplot(alpha_P, aes(x = Focal, y = Estimate, group = interaction(Treatment, Age))) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.3, position = position_dodge(0.5)) +
  geom_point(aes(fill = Treatment, shape = Age), size = 2, stroke = 0.3, 
             position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 22)) +
  labs(y = "Adult comp. competition") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 7, color = "black", hjust = 0))

# values for text
alpha %>%
  transmute(PP_eff = 100 * (alpha_PP_fung - alpha_PP) / alpha_PP,
            PF_eff = 100 * (alpha_FP_fung - alpha_FP) / alpha_FP,
            PA_eff = 100 * (alpha_FA_fung - alpha_FA) / alpha_FA) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate)

# table
write_csv(brms_SummaryTable(growthD2Mod2), "output/focal_growth_biomass_model_2019_density_exp.csv")


#### combine ####

# legend
leg <- get_legend(sev_fig)
top_row <- plot_grid(sev_fig + theme(legend.position = "none"), germ_fig, 
                     align = "h", axis = "l", ncol = 2,
                     labels = LETTERS[1:2],
                     label_size = 10)
mid_row1 <- plot_grid(est_fig, bio_fig,
                     align = "h", axis = "l", ncol = 2,
                     labels = LETTERS[3:4],
                     label_size = 10)
mid_row2 <- plot_grid(seed_fig, alphaA_fig, 
                        align = "h", axis = "l", ncol = 2,
                        labels = LETTERS[5:6],
                        label_size = 10)
bot_row <- plot_grid(alphaF_fig, alphaP_fig, 
                      align = "h", axis = "l", ncol = 2,
                      labels = LETTERS[7:8],
                      label_size = 10)
comb_fig <- plot_grid(top_row, mid_row1, mid_row2, bot_row, ncol = 1)

pdf("output/disease_treatment_effects_2018_2019_density_exp.pdf", width = 3.54, height = 6.29)
plot_grid(comb_fig, leg, nrow = 2, rel_heights = c(1, 0.07))
dev.off()
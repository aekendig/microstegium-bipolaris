##### outputs ####

# Figure 4


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(broom.mixed)
library(patchwork)

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

col_pal = c("black", "#238A8D")
shape_vals <- c(1, 2)

box_shade = "gray92"

dodge_width = 0.5


#### severity ####

# load model
load("output/focal_severity_model_aug_2019_dens_exp.rda")
load("output/focal_severity_model_jul_2019_dens_exp.rda")

# load data
sevD2Dat3_aug2 <- read_csv("output/focal_severity_model_data_aug_2019_dens_exp.csv")
sevD2Dat3_jul2 <- read_csv("output/focal_severity_model_data_jul_2019_dens_exp.csv")

# draws
sev_mv <- as_draws_df(sevD2Mod_sev_dens_aug) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(m_W = inv_logit_scaled(b_Intercept),
            m_F = inv_logit_scaled(b_Intercept + b_fungicide)) 

sev_ev <- as_draws_df(sevD2Mod_sev_dens_jul) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(s_W = inv_logit_scaled(b_Intercept + b_plant_groupEvS),
            s_F = inv_logit_scaled(b_Intercept + b_plant_groupEvS + b_fungicide + b_plant_groupEvS_fungicide),
            a_W = inv_logit_scaled(b_Intercept + b_plant_groupEvA),
            a_F = inv_logit_scaled(b_Intercept + b_plant_groupEvA + b_fungicide + b_plant_groupEvA_fungicide)) 

sev <- cbind(sev_mv, sev_ev)

# raw data
sev_raw_mv <- sevD2Dat3_aug2 %>%
  filter(plant_group == "Mv") %>%
  select(sp, age, site, treatment, next_severity_t)

sev_raw_ev <- sevD2Dat3_jul2 %>%
  filter(sp == "Ev") %>%
  select(sp, age, site, treatment, next_severity_t)

sev_raw <- sev_raw_mv %>%
  full_join(sev_raw_ev) %>%
  mutate(Species = if_else(sp == "Mv", "Invader (Mv)", "Competitor (Ev)") %>%
           fct_relevel("Invader (Mv)"),
         Age = fct_recode(age, 
                          "1st year" = "seedling") %>%
           fct_relevel("1st year"),
         Treatment = fct_recode(treatment,
                                "control" = "water") %>%
           fct_relevel("control"),
         SpeciesN = as.numeric(Species),
         Estimate = next_severity_t * 100)

# summarize for fig
sev2 <- sev %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Parameter = parameter, 
            Species = c(rep("Competitor (Ev)", 2), rep("Invader (Mv)", 2), rep("Competitor (Ev)", 2)) %>%
              fct_relevel("Invader (Mv)"),
            SpeciesN = as.numeric(Species),
            Age = c(rep("adult", 2), rep("1st year", 4)) %>%
                      fct_relevel("1st year"),
            Treatment = if_else(str_ends(parameter, "_W") == T, "control", "fungicide"),
            Estimate = estimate * 100,
            Lower = .lower * 100,
            Upper = .upper * 100)
  
# fig
sev_fig <- ggplot(sev2, aes(x = SpeciesN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = SpeciesN-0.5, xmax = Inf,
                fill = Species), alpha = 0.5) +
  geom_point(data = sev_raw, size = 0.25, alpha = 0.1,
             aes(color = Treatment), 
             position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75, 
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Disease severity (%)") +
  fig_theme +
  theme(axis.text.x = element_blank())

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
# did this in focal_severity_model_2018_2019.R too


#### germination ####

# load model and data
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
mvGermD1Dat <- read_csv("output/mv_germination_fungicide_model_data_2018_density_exp.csv")

# draws and summarize
mv_germ <- as_draws_df(mvGermD1Mod2)  %>%
  transmute(g_W = inv_logit_scaled(b_Intercept),
            g_F = inv_logit_scaled(b_Intercept + b_fungicide)) 

mv_germ2 <- mv_germ %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdi(estimate) %>%
  transmute(Species = "Invader (Mv)",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = 100 * estimate,
            Lower = 100 * .lower,
            Upper = 100 * .upper)

mv_germ_raw <- mvGermD1Dat %>%
  mutate(Estimate = 100 * (germination_final / seeds),
         Species = "Invader (Mv)",
         Treatment = treatment)

# load model
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")
evGermDat2 <- read_csv("output/ev_germination_fungicide_model_data_2018_2019_density_exp.rda")

# draws and summarize
ev_germ <- as_draws_df(evGermMod)  %>%
  transmute(g_W = inv_logit_scaled(b_Intercept),
            g_F = inv_logit_scaled(b_Intercept + b_fungicide)) 

ev_germ2 <- ev_germ %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Species = "Competitor (Ev)",
            Treatment = if_else(parameter == "g_W", "control", "fungicide"),
            Estimate = 100 * estimate,
            Lower = 100 * .lower,
            Upper = 100 * .upper)

ev_germ_raw <- evGermDat2 %>%
  mutate(Estimate = 100 * (germinants / seeds_planted),
         Species = "Competitor (Ev)",
         Treatment = treatment)

germ_raw <- mv_germ_raw %>%
  select(Estimate, Species, Treatment) %>%
  mutate(ID = 1:n()) %>%
  full_join(ev_germ_raw%>%
              select(Estimate, Species, Treatment) %>%
              mutate(ID = 1:n())) %>%
  mutate(Species = fct_relevel(Species, "Invader (Mv)"),
         SpeciesN = as.numeric(Species))

# combine
germ <- mv_germ2 %>%
  full_join(ev_germ2) %>%
  mutate(Species = fct_relevel(Species, "Invader (Mv)"),
         SpeciesN = as.numeric(Species))

# fig
germ_fig <- ggplot(germ, aes(x = SpeciesN, y = Estimate)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = SpeciesN-0.5, xmax = Inf,
                fill = Species), alpha = 0.5) +
  geom_point(data = germ_raw, size = 0.25, alpha = 0.1,
             aes(color = Treatment), 
             position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment), size = 2, shape = shape_vals[1], stroke = 0.75, 
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Germination (%)") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# values for text
mv_germ %>%
  transmute(eff = 100 * (g_F - g_W) / g_W) %>%
  mean_hdci(eff)

ev_germ %>%
  transmute(eff = 100 * (g_F - g_W) / g_W) %>%
  mean_hdci(eff)

# tables
write_csv(tidy(mvGermD1Mod2), "output/mv_germination_fungicide_model_2018_density_exp.csv")
write_csv(tidy(evGermMod), "output/ev_germination_fungicide_model_2018_2019_density_exp.csv")


#### establishment ####

# load model
load("output/survival_fungicide_model_2019_density_exp.rda")
survD2Dat2 <- read_csv("output/survival_fungicide_model_data_2019_density_exp.csv")

# draws
est <- as_draws_df(survFungD2Mod)  %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  transmute(e_A_W = inv_logit_scaled(b_Intercept),
            e_A_F = inv_logit_scaled(b_Intercept + b_fungicide),
            e_P_W = inv_logit_scaled(b_Intercept + b_focs),
            e_P_F = inv_logit_scaled(b_Intercept + b_focs + b_fungicide))

# summarize for fig
est2 <- est %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  mean_hdci(estimate) %>%
  transmute(Parameter = parameter,
            Species = rep(c("Invader (Mv)", "Competitor (Ev)"), each = 2),
            Treatment = if_else(str_detect(parameter, "W") == T, "control", "fungicide"),
            Estimate = 100 * estimate,
            Lower = 100 * .lower,
            Upper = 100 * .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader (Mv)"),
         SpeciesN = as.numeric(Species))

# raw data
est_raw <- survD2Dat2 %>%
  select(site, plot, treatment, sp, age, survival) %>%
  mutate(Estimate = 100 * survival,
         Species = fct_recode(sp,
                              "Invader (Mv)" = "Mv",
                              "Competitor (Ev)" = "Ev") %>%
           fct_relevel("Invader (Mv)"),
         Treatment = fct_recode(treatment, 
                                "control" = "water") %>%
           fct_relevel("control"),
         SpeciesN = as.numeric(Species))

# fig
est_fig <- ggplot(est2, aes(x = SpeciesN, y = Estimate)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = SpeciesN-0.5, xmax = Inf,
                fill = Species), alpha = 0.5) +
  geom_point(data = est_raw, size = 0.25, alpha = 0.1,
             aes(color = Treatment), 
             position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment), size = 2, shape = shape_vals[1], stroke = 0.75, 
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) +
  labs(y = "Establishment (%)") +
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
write_csv(tidy(survFungD2Mod), 
          "output/survival_fungicide_model_2019_density_exp.csv")


#### biomass ####

# load model
load("output/focal_growth_density_model_2019_density_exp.rda")
growthD2Dat2 <- read_csv("intermediate-data/focal_growth_density_data_2019_density_exp.csv")

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
            Species = rep(c("Invader (Mv)",
                            "Competitor (Ev)",
                            "Competitor (Ev)"), each = 2),
            Age = c(rep("1st year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader (Mv)"),
         SpeciesN = as.numeric(Species),
         Age = fct_relevel(Age, "1st year"))

# raw data
bio_raw <- growthD2Dat2 %>%
  mutate(Estimate = plant_growth,
         Species = fct_recode(sp,
                              "Invader (Mv)" = "Mv",
                              "Competitor (Ev)" = "Ev") %>%
           fct_relevel("Invader (Mv)"),
         Treatment = fct_recode(treatment, 
                                "control" = "water") %>%
           fct_relevel("control"),
         SpeciesN = as.numeric(Species),
         Age = fct_recode(age, 
                          "1st year" = "seedling") %>%
           fct_relevel("1st year"))

# fig
bio_fig <- ggplot(bio2, aes(x = SpeciesN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = SpeciesN-0.5, xmax = Inf,
                fill = Species), alpha = 0.5) +
  geom_point(data = bio_raw, size = 0.25, alpha = 0.1,
             aes(color = Treatment), 
             position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75, 
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Biomass (ln[g])") +
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
write_csv(tidy(growthD2Mod), 
          "output/focal_growth_density_model_2019_density_exp.csv")


#### seeds ####

# load model
load("output/focal_seed_density_model_2019_density_exp.rda")
seedD2Dat3 <- read_csv("intermediate-data/focal_seed_density_data_2019_density_exp.csv")

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
            Species = rep(c("Invader (Mv)",
                            "Competitor (Ev)",
                            "Competitor (Ev)"), each = 2),
            Age = c(rep("1st year", 4), rep("adult", 2)),
            Treatment = rep(c("control", "fungicide"), 3),
            Estimate = estimate,
            Lower = .lower,
            Upper = .upper) %>%
  mutate(Species = fct_relevel(Species, "Invader (Mv)"),
         SpeciesN = as.numeric(Species),
         Age = fct_relevel(Age, "1st year"))

# raw data
seed_raw <- seedD2Dat3 %>%
  mutate(Estimate = log_seeds,
         Species = fct_recode(sp,
                              "Invader (Mv)" = "Mv",
                              "Competitor (Ev)" = "Ev") %>%
           fct_relevel("Invader (Mv)"),
         Treatment = fct_recode(treatment, 
                                "control" = "water") %>%
           fct_relevel("control"),
         SpeciesN = as.numeric(Species),
         Age = fct_recode(age, 
                          "1st year" = "seedling") %>%
           fct_relevel("1st year"))

seed_fig <- ggplot(seed2, aes(x = SpeciesN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = SpeciesN-0.5, xmax = Inf,
                fill = Species), alpha = 0.5) +
  geom_point(data = seed_raw, size = 0.25, alpha = 0.1,
             aes(color = Treatment), 
             position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = dodge_width)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75,
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Seed production (ln[x+1])") +
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
write_csv(tidy(seedD2Mod), 
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
            Focal = if_else(str_starts(parameter, "alpha_A"), "Invader (Mv)", "Competitor (Ev)") %>%
              fct_relevel("Invader (Mv)"),
            FocalN = as.numeric(Focal),
            Age = if_else(str_starts(parameter, "alpha_P"), "adult", "1st year") %>%
              fct_relevel("1st year"),
            Competitor = case_when(str_ends(Parameter, "A") ~ "Invader",
                                   str_ends(Parameter, "F") ~ "1st year comp.",
                                   str_ends(Parameter, "P") ~ "Adult comp."),
            Estimate = -1 * estimate,
            Lower = -1 * .lower,
            Upper = -1 * .upper)

# separate by Competitor
alpha_A <- alpha2 %>% filter(Competitor == "Invader")
alpha_F <- alpha2 %>% filter(Competitor == "1st year comp.")
alpha_P <- alpha2 %>% filter(Competitor == "Adult comp.")

# figs
alphaA_fig <- ggplot(alpha_A, aes(x = FocalN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = FocalN-0.5, xmax = Inf,
                fill = Focal), alpha = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75,
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Invader (Mv) effect") +
  fig_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank())

alphaF_fig <- ggplot(alpha_F, aes(x = FocalN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = FocalN-0.5, xmax = Inf,
                fill = Focal), alpha = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75,
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "1st yr comp. (Ev) effect") +
  fig_theme +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 7, color = "black", hjust = 0))

alphaP_fig <- ggplot(alpha_P, aes(x = FocalN, y = Estimate, group = interaction(Treatment, Age))) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = FocalN-0.5, xmax = Inf,
                fill = Focal), alpha = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Treatment),
                width = 0, size = 0.5, position = position_dodge(dodge_width)) +
  geom_point(aes(color = Treatment, shape = Age), size = 2, stroke = 0.75,
             position = position_dodge(dodge_width)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("transparent", box_shade), guide = "none") +
  scale_shape_manual(values = shape_vals) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(sev2$Species),
                     limits = c(0.6, 2.4)) + 
  labs(y = "Adult comp. (Ev) effect") +
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

# find CI width that leads to sig competitive tolerance of Ev adult
alpha %>%
  transmute(PP_eff = 100 * (alpha_PP_fung - alpha_PP) / alpha_PP,
            PF_eff = 100 * (alpha_FP_fung - alpha_FP) / alpha_FP,
            PA_eff = 100 * (alpha_FA_fung - alpha_FA) / alpha_FA) %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") %>%
  group_by(parameter) %>%
  median_hdi(estimate, .width = 0.29)


# table
write_csv(tidy(growthD2Mod2), "output/focal_growth_biomass_model_2019_density_exp.csv")


#### combine ####

# legend
leg <- get_legend(sev_fig)

comb_fig <- (sev_fig + theme(legend.position = "none")) + 
  germ_fig +
  est_fig +
  bio_fig +
  seed_fig +
  alphaA_fig +
  alphaF_fig +
  alphaP_fig +
  plot_layout(nrow = 4) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

pdf("output/disease_treatment_effects_2018_2019_density_exp.pdf", width = 4.13, height = 6.29)
plot_grid(comb_fig, leg, nrow = 2, rel_heights = c(1, 0.07))
dev.off()
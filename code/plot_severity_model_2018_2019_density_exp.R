##### info ####

# file: plot_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 5/6/21
# goal: effects of within and outside of plot severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)

# import data
d1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("#472F7DFF", "#39568CFF", "#35B779FF")
shape_pal = c(23, 22, 24)


#### edit data ####

# # format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(!(month %in% c("sep"))) %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# one plot in one month is missing, use nearby plots
edgeSevD2Dat3 <- tibble(site = "D1", plot = 6, treatment = "fungicide", month = "early_aug") %>%
  mutate(edge_severity = edgeSevD2Dat2 %>%
           filter(site == "D1" & plot %in% c(4, 5) & treatment == "water" & month == "early_aug") %>%
           summarise(sev = mean(edge_severity)) %>%
           pull(sev)) %>%
  full_join(edgeSevD2Dat2) %>%
  filter(!is.na(edge_severity)) # use to check for missing values

# bg severity
bgSevD1Dat <- d1Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") | # select background measurements (includes focal)
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, prop_healthy))

bgSevD2Dat <- d2Dat %>%
  filter((plot %in% 2:4 & sp == "Mv") |
           (plot %in% 5:7 & sp == "Ev" & age == "seedling") |
           (plot %in% 8:10 & sp == "Ev" & age == "adult")) %>%
  rename(bg_lesions = lesions,
         bg_sp = sp,
         bg_age = age) %>%
  select(-c(severity, prop_healthy))

# focal dataset
focSevD1Dat <- d1Dat %>%
  filter(month != "jul") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  mutate(month = dplyr::recode(month, 
                               "late_aug" = "jul", # match prior month
                               "sep" = "late_aug")) %>%
  select(-c(severity, prop_healthy))

focSevD2Dat <- d2Dat %>%
  filter(month != "may") %>% # remove first month - no prior data
  rename(foc_sp = sp,
         foc_age = age,
         foc_lesions = lesions) %>%
  mutate(month = dplyr::recode(month,  # match prior month
                               "jun" = "may",
                               "jul" = "jun",
                               "early_aug" = "jul",
                               "late_aug" = "early_aug")) %>%
  select(-c(severity, prop_healthy))

# prior severity
priorSevD1Dat <- d1Dat %>%
  filter(month != "sep") %>% # remove last month
  mutate(foc_prior_disease = prop_healthy * lesions) %>%
  rename(foc_sp = sp,
         foc_age = age,
         foc_prior_prop_healthy = prop_healthy,
         foc_prior_lesions = lesions) %>%
  select(-severity)

priorSevD2Dat <- d2Dat %>%
  filter(month != "late_aug") %>% # remove last month
  mutate(foc_prior_disease = prop_healthy * lesions) %>%
  rename(foc_sp = sp,
         foc_age = age,
         foc_prior_prop_healthy = prop_healthy,
         foc_prior_lesions = lesions) %>%
  select(-severity)

# combine data
sevD1Dat <- focSevD1Dat %>%
  left_join(priorSevD1Dat) %>%
  left_join(bgSevD1Dat) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                           TRUE ~ bg_age),
         foc_sp_age = paste(foc_sp, foc_age, sep = "_"),
         bg_sp_age = paste(bg_sp, bg_age, sep = "_"),
         bg_lesions = case_when(foc_sp_age == bg_sp_age & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         foc_lesions_change = foc_lesions - foc_prior_lesions,
         foc_lesions_change = case_when(month == "jul" ~ foc_lesions_change / 2, # account for longer lapse between samples 
                                        TRUE ~ foc_lesions_change),
         bg_prior_disease = bg_lesions * foc_prior_prop_healthy) %>%
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_foc_ls_chg = mean(foc_lesions_change, na.rm = T),
         sd_foc_ls_chg = sd(foc_lesions_change, na.rm = T),
         mean_foc_pri_dis = mean(foc_prior_disease, na.rm = T),
         sd_foc_pri_dis = sd(foc_prior_disease, na.rm = T),
         mean_bg_pri_dis = mean(bg_prior_disease, na.rm = T),
         sd_bg_pri_dis = sd(bg_prior_disease, na.rm = T)) %>%
  ungroup() %>%
  mutate(foc_ls_chg_s = (foc_lesions_change - mean_foc_ls_chg) / sd_foc_ls_chg,
         foc_pr_dis_s = (foc_prior_disease - mean_foc_pri_dis) / sd_foc_pri_dis,
         bg_pr_dis_s = (bg_prior_disease - mean_bg_pri_dis) / sd_bg_pri_dis)

sevD2Dat <- focSevD2Dat %>%
  left_join(priorSevD2Dat) %>%
  left_join(bgSevD2Dat) %>%
  left_join(edgeSevD2Dat3) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         bg_sp = case_when(plot == 1 ~ foc_sp, # make the focal species the bg in 1 plots
                           TRUE ~ bg_sp),
         bg_age = case_when(plot == 1 ~ foc_age,
                            TRUE ~ bg_age),
         foc_sp_age = paste(foc_sp, foc_age, sep = "_"),
         bg_sp_age = paste(bg_sp, bg_age, sep = "_"),
         bg_lesions = case_when(foc_sp_age == bg_sp_age & plot == 1 ~ foc_lesions, # use focal lesions as bg in 1 plots
                                TRUE ~ bg_lesions),
         foc_lesions_change = foc_lesions - foc_prior_lesions,
         bg_prior_disease = bg_lesions * foc_prior_prop_healthy) %>%
  group_by(foc_sp, foc_age, bg_sp, bg_age) %>%
  mutate(mean_foc_ls_chg = mean(foc_lesions_change, na.rm = T),
         sd_foc_ls_chg = sd(foc_lesions_change, na.rm = T),
         mean_foc_pri_dis = mean(foc_prior_disease, na.rm = T),
         sd_foc_pri_dis = sd(foc_prior_disease, na.rm = T),
         mean_bg_pri_dis = mean(bg_prior_disease, na.rm = T),
         sd_bg_pri_dis = sd(bg_prior_disease, na.rm = T)) %>%
  ungroup() %>%
  mutate(foc_ls_chg_s = (foc_lesions_change - mean_foc_ls_chg) / sd_foc_ls_chg,
         foc_pr_dis_s = (foc_prior_disease - mean_foc_pri_dis) / sd_foc_pri_dis,
         bg_pr_dis_s = (bg_prior_disease - mean_bg_pri_dis) / sd_bg_pri_dis)

# modified bg data for figure
bgSevDat <- sevD1Dat %>%
  select(month, site, plot, treatment, bg_sp, bg_age, bg_sp_age, bg_lesions) %>%
  unique() %>%
  mutate(yearf = "2018") %>%
  full_join(sevD2Dat %>%
              select(month, site, plot, treatment, bg_sp, bg_age, bg_sp_age, bg_lesions) %>%
              unique() %>%
              mutate(yearf = "2019")) %>%
  filter(!is.na(bg_lesions)) %>%
  group_by(yearf, bg_sp_age) %>%
  mutate(mean_bg_ls = mean(bg_lesions),
         sd_bg_ls = sd(bg_lesions)) %>%
  ungroup() %>%
  mutate(bg_ls_s = (bg_lesions - mean_bg_ls) / sd_bg_ls,
         density_level = case_when(plot == 1 ~ "none",
                                   plot %in% c(2, 5, 8) ~ "low",
                                   plot %in% c(3, 6, 9) ~ "medium",
                                   plot %in% c(4, 7, 10) ~ "high") %>%
           fct_relevel("none", "low", "medium"),
         bg_sp_age = str_replace(bg_sp_age, "_", " "),
         bg_sp_age = recode(bg_sp_age, "Mv seedling" = "Mv"))


#### visualizations ####

# prior disease and lesions
ggplot(priorSevD1Dat, aes(foc_prior_lesions, foc_prior_disease)) +
  geom_point(aes(color = month)) +
  geom_smooth() +
  facet_wrap(foc_sp ~ foc_age)
# too correlated, have to pick one of these

# lesion change distributions
ggplot(sevD1Dat, aes(foc_lesions_change)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(foc_ls_chg_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# focal disease distributions
ggplot(sevD1Dat, aes(foc_prior_disease)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(foc_pr_dis_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# bg disease distributions
ggplot(sevD1Dat, aes(bg_prior_disease)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(bg_pr_dis_s)) +
  geom_density() +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# intraspecific
ggplot(sevD1Dat, aes(foc_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

ggplot(sevD2Dat, aes(foc_pr_dis_s, foc_ls_chg_s, color = treatment, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp ~ foc_age, scales = "free")

# interspecific
ggplot(sevD1Dat, aes(bg_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp_age ~ bg_sp_age, scales = "free")

ggplot(sevD2Dat, aes(bg_pr_dis_s, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(foc_sp_age ~ bg_sp_age, scales = "free")

# edge
ggplot(sevD2Dat, aes(edge_severity, foc_ls_chg_s, color = treatment)) +
  geom_point(aes(shape = month)) +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_wrap(~ foc_sp_age, scales = "free")


#### divide data ####

mvD1Dat <- sevD1Dat %>%
  filter(foc_sp == "Mv" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s)) %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Mv_seedling"))

evSD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_seedling" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s)) %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Ev_seedling"))

evAD1Dat <- sevD1Dat %>%
  filter(foc_sp_age == "Ev_adult" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s))

mvD2Dat <- sevD2Dat %>%
  filter(foc_sp == "Mv" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s) & !is.na(edge_severity)) %>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Mv_seedling"))

evSD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_seedling" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s) & !is.na(edge_severity))%>%
  mutate(bg_sp_age = fct_relevel(bg_sp_age, "Ev_seedling"))

evAD2Dat <- sevD2Dat %>%
  filter(foc_sp_age == "Ev_adult" & !is.na(foc_lesions_change) & !is.na(bg_pr_dis_s) & !is.na(edge_severity))


#### 2018 models ####

# Mv
mvSevD1Mod1 <- brm(foc_ls_chg_s ~ fungicide * bg_pr_dis_s * bg_sp_age + (1|plotf),
                   data = mvD1Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99))
summary(mvSevD1Mod1)

# Ev seedling
evSSevD1Mod1 <- update(mvSevD1Mod1, newdata = evSD1Dat)
summary(evSSevD1Mod1)

# Ev adult
evASevD1Mod1 <- update(mvSevD1Mod1, newdata = evAD1Dat)
summary(evASevD1Mod1)


#### 2019 models ####

# Mv
mvSevD2Mod1 <- brm(foc_ls_chg_s ~ fungicide * (bg_pr_dis_s * bg_sp_age + edge_severity) + (1|plotf),
                   data = mvD2Dat, family = gaussian,
                   prior <- c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99))
summary(mvSevD2Mod1)

# Ev seedling
evSSevD2Mod1 <- update(mvSevD2Mod1, newdata = evSD2Dat)
summary(evSSevD2Mod1)

# Ev adult
evASevD2Mod1 <- update(mvSevD2Mod1, newdata = evAD2Dat)
summary(evASevD2Mod1)

# no significant edge severity effects


#### extract effects ####

bg_eff_fun <- function(mod, foc, yr){
  
  if(foc == "Mv"){
    out <- posterior_samples(mod) %>%
      rename("Mv" = "b_bg_pr_dis_s", 
             "b_Ev_adult" = "b_bg_pr_dis_s:bg_sp_ageEv_adult",
             "b_Ev_seedling" = "b_bg_pr_dis_s:bg_sp_ageEv_seedling") %>%
      mutate(Ev_adult = Mv + b_Ev_adult,
             Ev_seedling = Mv + b_Ev_seedling)
  } else if(foc == "Ev seedling"){
    out <- posterior_samples(mod) %>%
      rename("Ev_seedling" = "b_bg_pr_dis_s", 
             "b_Ev_adult" = "b_bg_pr_dis_s:bg_sp_ageEv_adult",
             "b_Mv" = "b_bg_pr_dis_s:bg_sp_ageMv_seedling") %>%
      mutate(Ev_adult = Ev_seedling + b_Ev_adult,
             Mv = Ev_seedling + b_Mv)
  } else{
    out <- posterior_samples(mod) %>%
      rename("Ev_adult" = "b_bg_pr_dis_s", 
             "b_Ev_seedling" = "b_bg_pr_dis_s:bg_sp_ageEv_seedling",
             "b_Mv" = "b_bg_pr_dis_s:bg_sp_ageMv_seedling") %>%
      mutate(Ev_seedling = Ev_adult + b_Ev_seedling,
             Mv = Ev_adult + b_Mv)
  }

  
  out2 <- out %>%
    select(Mv, Ev_adult, Ev_seedling) %>%
    pivot_longer(cols = everything(),
                 names_to = "bg_sp_age",
                 values_to = "bg_eff") %>%
    group_by(bg_sp_age) %>%
    median_hdi(bg_eff) %>%
    mutate(focal = foc,
           background = str_replace(bg_sp_age, "_", " "),
           yearf = yr)
  
  return(out2)
  
}

bgEff <- bg_eff_fun(mvSevD1Mod1, "Mv", "2018") %>%
  full_join(bg_eff_fun(evSSevD1Mod1, "Ev seedling", "2018")) %>%
  full_join(bg_eff_fun(evASevD1Mod1, "Ev adult", "2018")) %>%
  full_join(bg_eff_fun(mvSevD2Mod1, "Mv", "2019")) %>%
  full_join(bg_eff_fun(evSSevD2Mod1, "Ev seedling", "2019")) %>%
  full_join(bg_eff_fun(evASevD2Mod1, "Ev adult", "2019")) %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"))


#### figure ####

bg_dens_figD1 <- ggplot(filter(bgSevDat, yearf == "2018"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.7)) +
  stat_summary(geom = "point", size = 2, fun = "mean", fill = "white", position = position_dodge(0.7), aes(shape = bg_sp_age)) +
  geom_text(x = 0.5, y = 1.03, label = "2018", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
  ylab("Background lesions") +
  xlab("Background plant density") +
  scale_color_manual(values = col_pal, name = "Plant group") +
  scale_shape_manual(values = shape_pal, name = "Plant group") +
  fig_theme +
  theme(legend.position = c(0.23, 0.75))

bg_dens_figD2 <- ggplot(filter(bgSevDat, yearf == "2019"), aes(density_level, bg_ls_s, color = bg_sp_age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(0.7)) +
  stat_summary(geom = "point", size = 2, fun = "mean", fill = "white", position = position_dodge(0.7), aes(shape = bg_sp_age)) +
  geom_text(x = 0.5, y = 1.3, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
  ylab("Background lesions") +
  xlab("Background plant density") +
  scale_color_manual(values = col_pal, name = "Plant group") +
  scale_shape_manual(values = shape_pal, name = "Plant group") +
  fig_theme

les_figD1 <- ggplot(filter(bgEff, yearf == "2018"), aes(background, bg_eff, color = focal)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, position = position_dodge(0.5)) +
  geom_point(aes(shape = focal, fill = sig), size = 2, position = position_dodge(0.5)) +
  geom_text(x = 0.48, y = 1.06, label = "2018", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
  xlab("Background plant group") +
  ylab("Change in focal lesions") +
  scale_color_manual(values = col_pal, name = "Focal plant group") +
  scale_shape_manual(values = shape_pal, name = "Focal plant group") +
  scale_fill_manual(values = c("white", "black"), guide = F) +
  fig_theme

les_figD2 <- ggplot(filter(bgEff, yearf == "2019"), aes(background, bg_eff, color = focal)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, position = position_dodge(0.5)) +
  geom_point(aes(shape = focal, fill = sig), size = 2, position = position_dodge(0.5)) +
  geom_text(x = 0.48, y = 1.18, label = "2019", color = "black", size = 2.5, check_overlap = T, hjust = 0) +
  xlab("Background plant group") +
  ylab("Change in focal lesions") +
  scale_color_manual(values = col_pal, name = "Focal plant group") +
  scale_shape_manual(values = shape_pal, name = "Focal plant group") +
  scale_fill_manual(values = c("white", "black"), guide = F) +
  fig_theme +
  theme(legend.position = c(0.6, 0.87))


# combine figure
pdf("output/plot_severity_effect_figure_2018_2019_density_exp.pdf", width = 4.7, height = 4.7)
plot_grid(bg_dens_figD1, bg_dens_figD2, les_figD1, les_figD2,
          nrow = 2,
          labels = LETTERS[1:4])
dev.off()
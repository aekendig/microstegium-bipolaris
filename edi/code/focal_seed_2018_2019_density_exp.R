##### outputs ####

# focal_seed_density_model_2018_density_exp.rda
# focal_seed_density_model_2019_density_exp.rda
# focal_seed_density_data_2018_density_exp.csv
# focal_seed_density_data_2019_density_exp.csv
# focal_seed_competition_coefficients_2018_2019_density_exp.csv (Table S5)

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(GGally)

# import plot information
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import seed data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")

# import survival data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")

# import tiller data
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### edit data ####

# plant group densities
plotDens <- plots %>%
  mutate(density = case_when(plot %in% 2:4 ~ background_density + 3,
                             plot %in% 5:7 ~ background_density + 3,
                             plot%in% 8:10 ~ background_density + 1,
                             plot == 1 ~ 0)) %>%
  select(plot, treatment, background, density_level, density)

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal)) %>%
  filter(!is.na(survival))

# 2019 list of all plants
# all dead plants were replaced
focD2Dat <- plots %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3))) %>%
  mutate(sp = "Mv",
         age = "seedling") %>%
  full_join(plots %>%
              select(plot, treatment) %>%
              expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
              expand_grid(tibble(ID = c("1", "2", "3", "A"),
                                 age = c(rep("seedling", 3), "adult"))) %>%
              mutate(sp = "Ev"))

# Ev seeds 2018
evSeedD1Dat2 <- evSeedD1Dat %>%
  filter(focal == 1 & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev" & survival == 1) %>%
              select(-survival)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# Ev seeds 2019
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(focD2Dat %>%
              filter(sp == "Ev")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# Mv seeds 2018
# quadrat was 0.49 x 0.25 m
mvSeedD1Dat2 <- mvSeedD1Dat %>% # none missing 
  left_join(tillerD1Dat %>%
              filter(sp == "Mv") %>%
              group_by(site, plot, treatment, sp) %>%
              summarise(tillers_jul = mean(tillers_jul, na.rm = T)) %>%
              ungroup()) %>%
  left_join(plotDens) %>%
  mutate(seeds_per_plant = seeds_per_stem * tillers_jul,
         seeds_per_m2 = (seeds_bio + seeds_soil) / (0.49 * 0.25),
         age = "seedling")

# Mv seeds 2019
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0))

# combine by year
seedD1Dat <- evSeedD1Dat2 %>%
  full_join(mvSeedD1Dat2 %>%
              rename(seeds = seeds_per_plant)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         seeds1 = seeds + 1,
         log_seeds = log(seeds1),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1)))

seedD2Dat <- evSeedD2Dat2 %>%
  full_join(mvSeedD2Dat2) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         seeds1 = seeds + 1,
         log_seeds = log(seeds1),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("s"),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1)))


#### fit models ####

# initial visualizations
ggplot(seedD1Dat %>% filter(seeds > 0), aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(1)) +
  facet_grid(focal ~ background, scales = "free")
# no non-zero seeds for Ev seedlings with low Ev seedling density

ggplot(seedD2Dat %>% filter(seeds > 0), aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se", position = position_dodge(1)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(1)) +
  facet_grid(focal ~ background, scales = "free")

# remove plot 1
seedD1Dat3 <- seedD1Dat %>%
  filter(plot > 1) %>%
  mutate(foc = fct_relevel(foc, "m"),
         bg = fct_relevel(bg, "m"))

seedD2Dat3 <- seedD2Dat %>%
  filter(plot > 1) %>%
  mutate(foc = fct_relevel(foc, "m"),
         bg = fct_relevel(bg, "m"))

# fit models
seedD1Mod <- brm(data = seedD1Dat3, family = gaussian,
                 log_seeds ~ foc * fungicide * (density + density:bg) + (1|plotf),
                 prior <- c(prior(normal(7, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedD1Mod)

seedD2Mod <- brm(data = seedD2Dat3, family = gaussian,
                 log_seeds ~ foc * fungicide * (density + density:bg) + (1|plotf),
                 prior <- c(prior(normal(7, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(seedD2Mod)

# save models and data
save(seedD1Mod, file = "output/focal_seed_density_model_2018_density_exp.rda")
save(seedD2Mod, file = "output/focal_seed_density_model_2019_density_exp.rda")
write_csv(seedD1Dat3, "intermediate-data/focal_seed_density_data_2018_density_exp.csv")
write_csv(seedD2Dat3, "intermediate-data/focal_seed_density_data_2019_density_exp.csv")

# load models
load("output/focal_seed_density_model_2018_density_exp.rda")
load("output/focal_seed_density_model_2019_density_exp.rda")


#### interaction coefficients (alphas) ####

# Mv background
mv_mv_ctrl_alpha = "density = 0"
mv_mv_fung_alpha = "density + fungicide:density = 0"
evS_mv_ctrl_alpha = "density + focs:density = 0"
evS_mv_fung_alpha = "density + fungicide:density + focs:density + focs:fungicide:density = 0"
evA_mv_ctrl_alpha = "density + foca:density = 0"
evA_mv_fung_alpha = "density + fungicide:density + foca:density + foca:fungicide:density = 0"

# EvS background
evS_evS_ctrl_alpha = "density + focs:density + density:bgs + focs:density:bgs = 0"
evS_evS_fung_alpha = "density + focs:density + density:bgs + focs:density:bgs + fungicide:density + focs:fungicide:density + fungicide:density:bgs + focs:fungicide:density:bgs = 0"
mv_evS_ctrl_alpha = "density +  density:bgs = 0"
mv_evS_fung_alpha = "density +  density:bgs + fungicide:density + fungicide:density:bgs = 0"
evA_evS_ctrl_alpha = "density +  density:bgs + foca:density + foca:density:bgs = 0"
evA_evS_fung_alpha = "density +  density:bgs + foca:density + foca:density:bgs + fungicide:density +  fungicide:density:bgs + foca:fungicide:density + foca:fungicide:density:bgs = 0"

# EvA background
evA_evA_ctrl_alpha = "density + foca:density + density:bga + foca:density:bga = 0"
evA_evA_fung_alpha = "density + foca:density + density:bga + foca:density:bga + fungicide:density + foca:fungicide:density + fungicide:density:bga + foca:fungicide:density:bga = 0"
mv_evA_ctrl_alpha = "density +  density:bga = 0"
mv_evA_fung_alpha = "density +  density:bga + fungicide:density + fungicide:density:bga = 0"
evS_evA_ctrl_alpha = "density + density:bga + focs:density + focs:density:bga = 0"
evS_evA_fung_alpha = "density + density:bga + focs:density + focs:density:bga + fungicide:density + fungicide:density:bga + focs:fungicide:density + focs:fungicide:density:bga = 0"

seedD1alphas <- hypothesis(seedD1Mod, 
                           c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                             evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                             evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                             evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                             mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                             evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                             evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                             mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                             evS_evA_ctrl_alpha, evS_evA_fung_alpha))

seedD2alphas <- hypothesis(seedD2Mod, 
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
alphaDat <- seedD1alphas[[1]] %>%
  mutate(year = "2018") %>%
  full_join(seedD2alphas[[1]] %>%
              mutate(year = "2019")) %>%
  mutate(foc_bg_trt = rep(c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                            "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                            "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                            "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                            "s_a_ctrl", "s_a_fung"), 2)) %>%
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
  left_join(seedD1Dat %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(treatment = fct_recode(treatment, "control" = "control (water)"),
         focal = fct_recode(focal, "Ev first-year" = "Ev seedling"),
         background = fct_recode(background, "Ev first-year" = "Ev seedling")) %>%
  select(year, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper, sig) %>%
  arrange(year, background, focal, treatment)

# save
write_csv(alphaDatSave, "output/focal_seed_competition_coefficients_2018_2019_density_exp.csv")

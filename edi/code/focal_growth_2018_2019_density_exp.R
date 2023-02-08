##### outputs ####

# focal_growth_density_data_2018_density_exp.csv
# focal_growth_density_data_2019_density_exp.csv
# focal_growth_density_model_2018_density_exp.rda
# focal_growth_density_model_2019_density_exp.rda
# focal_growth_biomass_model_2019_density_exp.rda
# focal_growth_competition_coefficients_supp_table_density_exp.csv (Table S3)
# focal_growth_pairwise_figure_2019_density_exp.pdf (Fig. S7)
# focal_growth_pairwise_figure_2018_density_exp.pdf (Fig. S8)


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)

# import plot information
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# import plot biomass data
plotBioD2Dat <- read_csv("intermediate-data/plot_biomass_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### edit data ####

# plant group densities
plotDens <- plots %>%
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1,
                             TRUE ~ background_density)) %>%
  select(plot, treatment, background, density)

# missing data
filter(tillerD1Dat, is.na(tillers_jul)) # 0
filter(tillerD1Dat, is.na(tillers_jun)) # 0
filter(tillerD1Dat, tillers_jul == 0 | tillers_jun == 0) # these plants are dead
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling
filter(plotBioD2Dat, is.na(biomass.g_m2)) # all no bg plots

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
         plotf = paste0(site, plot, str_sub(treatment, 1, 1))) %>%
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
         plotf = paste0(site, plot, str_sub(treatment, 1, 1))) %>%
  filter(!is.na(biomass_weight.g)) %>%
  left_join(plotBioD2Dat %>%
              select(site, treatment, plot, biomass.g_m2) %>%
              rename(plot_biomass = biomass.g_m2)) %>%
  mutate(plot_biomass = if_else(as.character(focal) == as.character(background), 
                                plot_biomass - biomass_weight.g, plot_biomass))


#### fit models ####

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
  filter(plot > 1) %>%
  mutate(background = fct_drop(background))

growthD2Dat2 <- growthD2Dat %>%
  filter(plot > 1) %>%
  mutate(background = fct_drop(background))

# biomass visualization
ggplot(growthD2Dat2, aes(plot_biomass, plant_growth, color = treatment)) +
  geom_point() +
  facet_grid(focal ~ background, scales = "free")

# fit models
growthD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (density + density:bg) + (1|plotf),
                   prior <- c(prior(normal(0.75, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.9999, max_treedepth = 15)) 
mod_check_fun(growthD1Mod)

growthD2Mod <- brm(data = growthD2Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (density + density:bg) + (1|plotf),
                   prior <- c(prior(normal(3, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99999, max_treedepth = 15)) 
mod_check_fun(growthD2Mod)

# biomass model
growthD2Mod2 <- brm(data = growthD2Dat2, family = gaussian,
                   plant_growth ~ foc * fungicide * (plot_biomass + plot_biomass:bg) + (1|plotf),
                   prior <- c(prior(normal(3, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                   control = list(adapt_delta = 0.99999, max_treedepth = 15)) 
mod_check_fun(growthD2Mod2)

# subset data to remove high EvA biomass
growthD2Dat2b <- growthD2Dat2 %>%
  filter(!(bg == "a" & plot_biomass > 150))

growthD2Mod2b <- update(growthD2Mod2, newdata = growthD2Dat2b)
summary(growthD2Mod2b)

# save models and data
save(growthD1Mod, file = "output/focal_growth_density_model_2018_density_exp.rda")
save(growthD2Mod, file = "output/focal_growth_density_model_2019_density_exp.rda")
save(growthD2Mod2, file = "output/focal_growth_biomass_model_2019_density_exp.rda")
save(growthD2Mod2b, file = "output/focal_growth_biomass_model_no_high_EvA_2019_density_exp.rda")
write_csv(growthD1Dat2, "intermediate-data/focal_growth_density_data_2018_density_exp.csv")
write_csv(growthD2Dat2, "intermediate-data/focal_growth_density_data_2019_density_exp.csv")

# load models
load("output/focal_growth_density_model_2018_density_exp.rda")
load("output/focal_growth_density_model_2019_density_exp.rda")
load("output/focal_growth_biomass_model_2019_density_exp.rda")
load("output/focal_growth_biomass_model_no_high_EvA_2019_density_exp.rda")


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

densityD1alphas <- hypothesis(growthD1Mod, 
                              c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                                evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                                evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                                evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                                mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                                evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                                evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                                mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                                evS_evA_ctrl_alpha, evS_evA_fung_alpha))

densityD2alphas <- hypothesis(growthD2Mod, 
                             c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                               evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                               evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                               evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                               mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                               evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                               evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                               mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                               evS_evA_ctrl_alpha, evS_evA_fung_alpha))

# Mv background
mv_mv_ctrl_alpha = "plot_biomass = 0"
mv_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass = 0"
evS_mv_ctrl_alpha = "plot_biomass + focs:plot_biomass = 0"
evS_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + focs:plot_biomass + focs:fungicide:plot_biomass = 0"
evA_mv_ctrl_alpha = "plot_biomass + foca:plot_biomass = 0"
evA_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + foca:plot_biomass + foca:fungicide:plot_biomass = 0"

# EvS background
evS_evS_ctrl_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs = 0"
evS_evS_fung_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs + fungicide:plot_biomass + focs:fungicide:plot_biomass + fungicide:plot_biomass:bgs + focs:fungicide:plot_biomass:bgs = 0"
mv_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs = 0"
mv_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + fungicide:plot_biomass + fungicide:plot_biomass:bgs = 0"
evA_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs = 0"
evA_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs + fungicide:plot_biomass +  fungicide:plot_biomass:bgs + foca:fungicide:plot_biomass + foca:fungicide:plot_biomass:bgs = 0"

# EvA background
evA_evA_ctrl_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga = 0"
evA_evA_fung_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga + fungicide:plot_biomass + foca:fungicide:plot_biomass + fungicide:plot_biomass:bga + foca:fungicide:plot_biomass:bga = 0"
mv_evA_ctrl_alpha = "plot_biomass +  plot_biomass:bga = 0"
mv_evA_fung_alpha = "plot_biomass +  plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga = 0"
evS_evA_ctrl_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga = 0"
evS_evA_fung_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga + focs:fungicide:plot_biomass + focs:fungicide:plot_biomass:bga = 0"

biomassD2alphas <- hypothesis(growthD2Mod2, 
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
alphaDat <- densityD1alphas[[1]] %>%
  mutate(year = "2018",
         response = "tillers",
         gradient = "density",
         data = "full") %>%
  full_join(densityD2alphas[[1]] %>%
              mutate(year = "2019",
                     response = "biomass",
                     gradient = "density",
                     data = "full")) %>%
  full_join(biomassD2alphas[[1]] %>%
              mutate(year = "2019",
                     response = "biomass",
                     gradient = "biomass",
                     data = "full")) %>%
  mutate(foc_bg_trt = rep(c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                            "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                            "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                            "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                            "s_a_ctrl", "s_a_fung"), 3)) %>%
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
  left_join(growthD2Dat2 %>%
              select(foc, focal, bg, background) %>%
              unique()) %>% 
  mutate(focal = case_when(focal == "Mv" ~ "Invader (Mv)",
                           focal == "Ev seedling" ~ "1st yr competitor (Ev)",
                           focal == "Ev adult" ~ "Adult competitor (Ev)") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         background = case_when(background == "Mv" ~ "Invader (Mv)",
                                background == "Ev seedling" ~ "1st yr competitor (Ev)",
                                background == "Ev adult" ~ "Adult competitor (Ev)") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)")) %>%
  select(year, response, gradient, data, focal, background, treatment, Estimate, Est.Error, CI.Lower, CI.Upper, sig) %>%
  arrange(year, response, gradient, data, background, focal, treatment)

# save
write_csv(filter(alphaDatSave, year == 2019 & gradient == "density"), 
          "output/focal_growth_competition_coefficients_supp_table_density_exp.csv")


#### figure ####

# density function
dens_fun <- function(f, b, trt, yr){
  
  if(yr == "2018"){
    dat <- growthD1Dat2 %>% filter(foc == f & bg == b & treatment == trt)
    density <- seq(min(dat$density), max(dat$density), length.out = 100)
  }else{
    dat <- growthD2Dat2 %>% filter(foc == f & bg == b & treatment == trt)
    density <- seq(min(dat$plot_biomass), max(dat$plot_biomass), length.out = 100)
  }
  
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
  rename(plot_biomass = density) %>%
  mutate(plant_growth = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q97.5"])

predDat <- predD1Dat %>%
  full_join(predD2Dat) %>%
  mutate(focal = fct_recode(foc, "Invader (Mv)" = "m", 
                            "1st yr competitor (Ev)" = "s", 
                            "Adult competitor (Ev)" = "a") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         background = fct_recode(bg, "Invader (Mv)" = "m", 
                                 "1st yr competitor (Ev)" = "s", 
                                 "Adult competitor (Ev)" = "a") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

# combine with preddat
predDat2 <- predDat %>%
  left_join(alphaDatSave %>%
              filter(data == "full" &  (year == "2018" | gradient == "biomass"))) %>%
  mutate(intra = case_when(foc == bg ~ "yes",
                           foc %in% c("a", "s") & bg %in% c("a", "s") ~ "yes",
                           TRUE ~ "no"))

# raw data
figDat <- growthD1Dat2 %>%
  select(site, plot, treatment, sp, age, ID, focal, background, density, plant_growth) %>%
  mutate(year = "2018") %>%
  full_join(growthD2Dat2 %>%
              select(site, plot, treatment, sp, age, ID, focal, background, plot_biomass, plant_growth) %>%
              mutate(year = "2019")) %>%
  left_join(plots %>%
              select(plot, density_level) %>%
              unique()) %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         focal = case_when(focal == "Mv" ~ "Invader (Mv)",
                           focal == "Ev seedling" ~ "1st yr competitor (Ev)",
                           focal == "Ev adult" ~ "Adult competitor (Ev)") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         background = case_when(background == "Mv" ~ "Invader (Mv)",
                                background == "Ev seedling" ~ "1st yr competitor (Ev)",
                                background == "Ev adult" ~ "Adult competitor (Ev)") %>%
           fct_relevel("Invader (Mv)", "1st yr competitor (Ev)"),
         density_level = fct_relevel(density_level, "low", "medium"),
         intra = case_when(focal ==  background ~ "yes",
                           str_detect(focal, "Ev") == T &  str_detect(background, "Ev") == T ~ "yes",
                           TRUE ~ "no"))

# alphas for figure
# sig beta values
alphaDat2 <- alphaDatSave %>%
  filter(data == "full" &  (year == "2018" | gradient == "biomass") & sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, foc, focal, bg, background) %>%
              summarise(density = max(density, na.rm = T),
                        plot_biomass = max(plot_biomass, na.rm = T))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(lower = min(CI.Lower)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(growth = min(plant_growth)) %>%
                          ungroup()) %>%
              rowwise() %>%
              mutate(plant_growth = min(c(growth, lower))) %>%
              ungroup() %>%
              select(-c(growth, lower))) %>%
  mutate(comp = as.character(sprintf("%.3f",round(Estimate, 3))))


# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        strip.placement = "outside")

col_pal = c("black", "#238A8D")

textSize = 2.5

# legend
legFig <- ggplot(predDat2, aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = figDat, aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.5) +
  facet_grid(rows = vars(focal),
             cols = vars(background)) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density treatment") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical") +
  guides(shape = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

leg <- get_legend(legFig)


# 2018
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat2, year == "2018"), 
            aes(label = paste("alpha", "==", comp, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 0.2, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Background density (", m^-2, ")", sep = ""))) +
  ylab(expression(paste("Focal growth (ln[", tillers[July], "/", tillers[June], "])", sep = ""))) +
  fig_theme

# print
pdf("output/focal_growth_pairwise_figure_2018_density_exp.pdf", width = 3.94, height = 4.33)
plot_grid(pairD1Fig, leg,
          nrow = 2, 
          rel_heights = c(1, 0.15))
dev.off()

# 2019
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = plot_biomass, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat2, year == "2019"), 
            aes(label = paste("alpha", "==", comp, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 0.2, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  xlab(expression(paste("Background biomass (g/", m^2, ")", sep = ""))) +
  ylab("Focal biomass (ln[g])") +
  fig_theme

# combine
pdf("output/focal_growth_pairwise_figure_2019_density_exp.pdf", width = 3.94, height = 4.33)
plot_grid(pairD2Fig, leg,
          nrow = 2, 
          rel_heights = c(1, 0.15))
dev.off()

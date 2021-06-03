##### info ####

# file: focal_seed_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/2/21
# goal: analyses of seeds as a function of density

#### approaches ####

# Ricker model:
# log(Nt+1/Nt) = r - alphaNN x Nt - alphaNM x Mt


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import seed data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
# ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R

# import survival data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")

# import tiller data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
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
focD2Dat <- plotsD %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3))) %>%
  mutate(sp = "Mv",
         age = "seedling") %>%
  full_join(plotsD %>%
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
  mutate(seeds = replace_na(seeds, 0),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_seeds = log(seeds + 1),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev())

# Ev seeds 2019
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(focD2Dat %>%
              filter(sp == "Ev")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_seeds = log(seeds + 1),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev())

# Mv seeds 2018
# quadrat was 0.49 x 0.25 m
mvSeedD1Dat2 <- mvSeedD1Dat %>% # none missing 
  left_join(growthD1Dat %>%
              filter(sp == "Mv") %>%
              group_by(site, plot, treatment, sp) %>%
              summarise(tillers_jul = mean(tillers_jul, na.rm = T)) %>%
              ungroup()) %>%
  left_join(plotDens) %>%
  mutate(seeds_per_plant = seeds_per_stem * tillers_jul,
         seeds_per_m2 = (seeds_bio + seeds_soil) / (0.49 * 0.25),
         log_seeds_per_plant = log(seeds_per_plant + 1),
         log_seeds_per_m2 = log(seeds_per_m2 + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         age = "seedling",
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev())

# Mv seeds 2019
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_seeds = log(seeds + 1),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev())


#### 2019 Mv models ####

# initial visualization
ggplot(mvSeedD2Dat2, aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
mvSeedD2Dat3 <- mvSeedD2Dat2 %>%
  filter(plot > 1) %>%
  mutate(background = fct_drop(background))

mvSeedD2Mod <- brm(data = mvSeedD2Dat3, family = gaussian,
                   log_seeds ~ fungicide * density * background + (1|site),
                   prior <- c(prior(normal(7, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.99)) 
mod_check_fun(mvSeedD2Mod)

# simulated data
mvSeedD2Sim <- mod_fit_fun2(dat = mvSeedD2Dat3, mod = mvSeedD2Mod, treatCol1 = fungicide, treatCol2 = background,
                         xCol = density, minX = 0, maxX = 67, yCol = log_seeds, f2t = T)
mvSeedD2Sim[[2]]

# save models
save(mvSeedD2Mod, file = "output/mv_seeds_Ricker_model_2019_density_exp.rda")


#### 2019 Ev seedling models ####

# initial visualization
ggplot(evSeedD2Dat2, aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
evSSeedD2Dat3 <- evSeedD2Dat2 %>%
  filter(plot > 1 & age == "seedling") %>%
  mutate(background = fct_drop(background))

# fit model
evSSeedD2Mod <- update(mvSeedD2Mod, newdata = evSSeedD2Dat3,
                       prior = set_prior("normal(2.5, 1)", class = "Intercept"),
                       control = list(adapt_delta = 0.999))
mod_check_fun(evSSeedD2Mod)
# bimodal distribution because there are so many zero seeds
# could do a different model that takes this into account

# simulated data
evSSeedD2Sim <- mod_fit_fun2(dat = evSSeedD2Dat3, mod = evSSeedD2Mod, treatCol1 = fungicide, treatCol2 = background,
                            xCol = density, minX = 0, maxX = 67, yCol = log_seeds, f2t = T)
evSSeedD2Sim[[2]]
# still provides reasonable estimates despite bimodal distribution

# save models
save(evSSeedD2Mod, file = "output/ev_seedling_seeds_Ricker_model_2019_density_exp.rda")


#### 2019 Ev adult models ####

# initial visualization
ggplot(evSeedD2Dat2, aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
evASeedD2Dat3 <- evSeedD2Dat2 %>%
  filter(plot > 1 & age == "adult") %>%
  mutate(background = fct_drop(background))

# fit model
evASeedD2Mod <- update(evSSeedD2Mod, newdata = evASeedD2Dat3,
                       prior = set_prior("normal(3.5, 1)", class = "Intercept"))
mod_check_fun(evASeedD2Mod)

# simulated data
evASeedD2Sim <- mod_fit_fun2(dat = evASeedD2Dat3, mod = evASeedD2Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = log_seeds, f2t = T)
evASeedD2Sim[[2]]

# save models
save(evASeedD2Mod, file = "output/ev_adult_seeds_Ricker_model_2019_density_exp.rda")


#### 2018 Mv models ####

# initial visualization
ggplot(mvSeedD1Dat2, aes(x = density, y = log_seeds_per_plant, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
mvSeedD1Dat3 <- mvSeedD1Dat2 %>%
  filter(plot > 1) %>%
  mutate(background = fct_drop(background),
         log_seeds = log_seeds_per_plant)

# fit model
mvSeedD1Mod <- update(evSSeedD2Mod, newdata = mvSeedD1Dat3,
                       prior = set_prior("normal(7, 1)", class = "Intercept"))
mod_check_fun(mvSeedD1Mod)

# simulated data
mvSeedD1Sim <- mod_fit_fun2(dat = mvSeedD1Dat3, mod = mvSeedD1Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = log_seeds_per_plant, f2t = T)
mvSeedD1Sim[[2]]

# save models
save(mvSeedD1Mod, file = "output/mv_seeds_Ricker_model_2018_density_exp.rda")


#### 2018 Ev seedling models ####

# initial visualization
ggplot(evSeedD1Dat2, aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
evSSeedD1Dat3 <- evSeedD1Dat2 %>%
  filter(plot > 1 & age == "seedling") %>%
  mutate(background = fct_drop(background))

# fit model
evSSeedD1Mod <- update(evSSeedD2Mod, newdata = evSSeedD1Dat3,
                       prior = set_prior("normal(1, 1)", class = "Intercept"))
mod_check_fun(evSSeedD1Mod)
# bimodal distribution because there are so many zero seeds
# could do a different model that takes this into account

# simulated data
evSSeedD1Sim <- mod_fit_fun2(dat = evSSeedD1Dat3, mod = evSSeedD1Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = log_seeds, f2t = T)
evSSeedD1Sim[[2]]
# still provides reasonable estimates despite bimodal distribution

# save models
save(evSSeedD1Mod, file = "output/ev_seedling_seeds_Ricker_model_2018_density_exp.rda")


#### 2018 Ev adult models ####

# initial visualization
ggplot(evSeedD1Dat2, aes(x = density, y = log_seeds, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_se") +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(age ~ background, scales = "free_x")

# remove none
# make background a factor
evASeedD1Dat3 <- evSeedD1Dat2 %>%
  filter(plot > 1 & age == "adult") %>%
  mutate(background = fct_drop(background))

# fit model
evASeedD1Mod <- update(evSSeedD2Mod, newdata = evASeedD1Dat3,
                       prior = set_prior("normal(3.5, 1)", class = "Intercept"))
mod_check_fun(evASeedD1Mod)

# simulated data
evASeedD1Sim <- mod_fit_fun2(dat = evASeedD1Dat3, mod = evASeedD1Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = log_seeds, f2t = T)
evASeedD1Sim[[2]]

# save models
save(evASeedD1Mod, file = "output/ev_adult_seeds_Ricker_model_2018_density_exp.rda")


#### figure ####

# prediction function
pair_pred_fun <- function(year, foc, bg, trt){
  
  # model and data
  if(year == "2018" & foc == "Mv_seedling"){
    mod <- mvSeedD1Mod
    dat <- mvSeedD1Dat3
  }else if(year == "2018" & foc == "Ev_seedling"){
    mod <- evSSeedD1Mod
    dat <- evSSeedD1Dat3
  }else if(year == "2018" & foc == "Ev_adult"){
    mod <- evASeedD1Mod
    dat <- evASeedD1Dat3
  }else if(year == "2019" & foc == "Mv_seedling"){
    mod <- mvSeedD2Mod
    dat <- mvSeedD2Dat3
  }else if(year == "2019" & foc == "Ev_seedling"){
    mod <- evSSeedD2Mod
    dat <- evSSeedD2Dat3
  }else if(year == "2019" & foc == "Ev_adult"){
    mod <- evASeedD2Mod
    dat <- evASeedD2Dat3
  }
  
  # select for background
  dat2 <- dat %>%
    filter(background == bg & treatment == trt)
  
  # sequence of density values
  # other variables
  # predicted values
  simDat <- tibble(density = seq(min(dat2$density), max(dat2$density), length.out = 100)) %>%
    mutate(background = bg,
           site = "A",
           fungicide = case_when(trt == "water" ~ 0,
                                 trt == "fungicide" ~ 1)) %>%
    mutate(log_seeds = fitted(mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(mod, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
    select(-background)
  
  # output
  return(simDat)
  
}

# predicted data
predDat <- tibble(focal = c("Ev_adult", "Ev_seedling", "Mv_seedling")) %>%
  expand_grid(tibble(background = c("Ev_adult", "Ev_seedling", "Mv_seedling"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  expand_grid(tibble(year = c("2018", "2019"))) %>%
  mutate(pred = pmap(list(year, focal, background, treatment), pair_pred_fun)) %>%
  unnest(pred) %>%
  mutate(focal = fct_recode(focal, "Mv" = "Mv_seedling") %>%
           str_replace("_", " "),
         background = fct_recode(background, "Mv" = "Mv_seedling") %>%
           str_replace("_", " "),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# combine alphas
alphasDat <- posterior_samples(mvSeedD1Mod) %>%
  mutate(focal = "Mv", year = "2018") %>%
  full_join(posterior_samples(evSSeedD1Mod) %>%
              mutate(focal = "Ev seedling", year = "2018")) %>%
  full_join(posterior_samples(evASeedD1Mod) %>%
              mutate(focal = "Ev adult", year = "2018")) %>%
  full_join(posterior_samples(mvSeedD2Mod) %>%
              mutate(focal = "Mv", year = "2019")) %>%
  full_join(posterior_samples(evSSeedD2Mod) %>%
              mutate(focal = "Ev seedling", year = "2019")) %>%
  full_join(posterior_samples(evASeedD2Mod) %>%
              mutate(focal = "Ev adult", year = "2019")) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  mutate(mv_water = b_density,
         evS_water = b_density + b_density_backgroundEv_seedling,
         evA_water = b_density + b_density_backgroundEv_adult,
         mv_fungicide = b_density + b_fungicide_density,
         evS_fungicide = evS_water + b_fungicide_density + b_fungicide_density_backgroundEv_seedling,
         evA_fungicide = evA_water + b_fungicide_density + b_fungicide_density_backgroundEv_adult) %>%
  select(year, focal, mv_water, evS_water, evA_water, mv_fungicide, evS_fungicide, evA_fungicide) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(year, focal),
               names_to = c(".value", "treatment"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_longer(cols = -c(year, focal, treatment),
               names_to = "background",
               values_to = "alpha") %>%
  mutate(alpha = alpha) %>%
  group_by(year, treatment, focal, background) %>%
  mean_hdi(alpha) %>%
  ungroup() %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"),
         background = fct_recode(background, "Ev adult" = "evA", "Ev seedling" = "evS", "Mv" = "mv"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# add to predicted dataset
predDat2 <- predDat %>%
  inner_join(alphasDat %>%
              select(year, treatment, focal, background, sig))

# raw data
figDat <- mvSeedD1Dat3 %>%
  select(site, plot, treatment, sp, age, plant_group, background, density, log_seeds) %>%
  mutate(year = "2018") %>%
  full_join(evSSeedD1Dat3 %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, log_seeds) %>%
              mutate(year = "2018"))%>%
  full_join(evASeedD1Dat3 %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, log_seeds) %>%
              mutate(year = "2018"))%>%
  full_join(mvSeedD2Dat3 %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, log_seeds) %>%
              mutate(year = "2019")) %>%
  full_join(evSSeedD2Dat3 %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, log_seeds) %>%
              mutate(year = "2019"))%>%
  full_join(evASeedD2Dat3 %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, log_seeds) %>%
              mutate(year = "2019"))%>%
  mutate(focal = str_replace(plant_group, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# sig alpha values
alphasDat2 <- alphasDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, background) %>%
              summarise(density = max(density))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(upper = max(upper)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(raw = max(log_seeds))) %>%
              rowwise() %>%
              mutate(log_seeds = max(c(raw, upper))) %>%
              ungroup() %>%
              select(-c(raw, upper))) %>%
  rename(param = alpha) %>%
  mutate(param = round(param, 2))

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
                   density = c(3, 3),
                   log_seeds = c(5.7, 6),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# water figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = log_seeds)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2018"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Seeds") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# fungicide figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = density, y = log_seeds)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2019"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Density (", m^-2, ")", sep = ""))) +
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
pdf("output/focal_seeds_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()

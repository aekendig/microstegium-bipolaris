##### outputs ####

# Figure S6

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)
library(tidybayes)

# import data
bgBioD2Dat <- read_csv("data/bg_biomass_2019_density_exp.csv") 
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### seeds ####

# planted density
plots
evPlantsD2Dat <- tibble(plot = c(1, 5:7, 1, 8:10),
                        planted = c(3, 7, 11, 19, 
                                    1, 3, 5, 9),
                        age = rep(c("seedling", "adult"), each = 4))

# format seeds
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>% # complete list of Ev
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  filter((plot %in% c(1, 5:7) & age == "seedling") | (plot %in% c(1, 8:10) & age == "adult")) %>% # only use the type planted as background
  group_by(site, plot, treatment, age) %>%
  summarise(seeds_per_plant = mean(seeds)) %>%
  ungroup() %>%
  left_join(evPlantsD2Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-seeds_per_plant)

# check
filter(evSeedD2Dat2, is.na(seeds))

# figures
ggplot(evSeedD2Dat2, aes(treatment, seeds)) +
  geom_boxplot() +
  facet_wrap(~ age, scales = "free")

ggplot(evSeedD2Dat2, aes(x = plot, y = seeds)) +
  geom_point() +
  facet_grid(age ~ treatment, scales = "free")


##### biomass ####

# zero background biomass
bgZero <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 2),
               treatment = rep(c("water", "fungicide"), 4)) %>%
  mutate(plot = 1,
         biomass_bg = 0)

# background biomass: combine plots
bgBioD2Dat2 <- bgBioD2Dat %>%
  group_by(site, plot, treatment) %>%
  summarise(biomass_bg = sum(biomass.g)) %>%
  ungroup() %>%
  full_join(bgZero)

# use average of others in plot if plant is missing biomass
evBioD2Dat2 <- bgBioD2Dat2 %>%
  inner_join(evBioD2Dat %>%
               filter(ID %in% c("1", "2", "3") & plot %in% c(1, 5:7)) %>%
               group_by(site, plot, treatment) %>%
               mutate(weight_adj = mean(weight, na.rm = T)) %>%
               ungroup() %>%
               mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                         TRUE ~ weight),
                      age = "seedling") %>%
               group_by(site, plot, treatment, age) %>%
               summarise(biomass_foc = sum(weight)) %>%
               ungroup() %>%
               full_join(evBioD2Dat %>%
                           filter(ID == "A" & plot %in% c(1, 8:10)) %>%
                           select(site, plot, treatment, weight) %>%
                           mutate(age = "adult") %>%
                           rename(biomass_foc = weight))) %>%
  mutate(biomass = biomass_bg + biomass_foc)

# figuress
ggplot(evBioD2Dat2, aes(treatment, biomass)) +
  geom_boxplot() +
  facet_wrap(~ age, scales = "free")

ggplot(evBioD2Dat2, aes(x = plot, y = biomass)) +
  geom_point() +
  facet_grid(age ~ treatment, scales = "free")


#### combine ####

# combine data
plotD2Dat <- evSeedD2Dat2 %>%
  full_join(evBioD2Dat2) %>%
  rename(density = planted) %>%
  mutate(fungicide = if_else(treatment == "fungicide", 1, 0),
         treatment = fct_recode(treatment, control = "water") %>%
           fct_relevel("control"))

# split by age
seedlingD2Dat <- plotD2Dat %>%
  filter(age == "seedling")

adultD2Dat <- plotD2Dat %>%
  filter(age == "adult")

# data for figure
figD2Dat <- plotD2Dat %>%
  select(site, plot, treatment, age, density, seeds, biomass) %>%
  pivot_longer(cols = c(seeds, biomass),
               names_to = "response",
               values_to = "value") %>%
  mutate(plant_group = if_else(age == "adult", "Adult competitor (Ev)",
                               "1st yr competitor (Ev)") %>%
           fct_relevel("1st yr competitor (Ev)"))


#### fit models ####

# initial visualization
ggplot(plotD2Dat, aes(x = density, y = biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ age, scales = "free")

ggplot(plotD2Dat, aes(x = density, y = seeds, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ age, scales = "free")

# priors
plotD2Dat %>%
  filter(plot == 1) %>%
  group_by(treatment, age) %>%
  summarise(b0 = mean(biomass/density),
            s0 = mean(seeds/density))

x <- seq(-1, 20, length.out = 100)
y <- dgamma(x, shape = 5, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# models
evSBioDensMod <- brm(data = seedlingD2Dat, family = gaussian,
                    bf(biomass ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(1, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 

evSSeedDensMod <- brm(data = seedlingD2Dat, family = gaussian,
                     bf(seeds ~ (density * b0)/(1 + alpha * density), 
                        b0 ~ 0 + treatment + (1|site), 
                        alpha ~ 0 + treatment, 
                        nl = T),
                     prior <- c(prior(gamma(5, 1), nlpar = "b0", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 

evABioDensMod <- brm(data = adultD2Dat, family = gaussian,
                     bf(biomass ~ (density * b0)/(1 + alpha * density), 
                        b0 ~ 0 + treatment + (1|site), 
                        alpha ~ 0 + treatment, 
                        nl = T),
                     prior <- c(prior(gamma(5, 1), nlpar = "b0", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 

evASeedDensMod <- brm(data = adultD2Dat, family = gaussian,
                      bf(seeds ~ (density * b0)/(1 + alpha * density), 
                         b0 ~ 0 + treatment + (1|site), 
                         alpha ~ 0 + treatment, 
                         nl = T),
                      prior <- c(prior(normal(80, 10), nlpar = "b0"),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                      iter = 6000, warmup = 1000, chains = 3,
                      control = list(adapt_delta = 0.999)) 

mod_check_fun(evSBioDensMod)
mod_check_fun(evSSeedDensMod)
mod_check_fun(evABioDensMod)
mod_check_fun(evASeedDensMod)

# save models
save(evSBioDensMod, file = "output/evS_plot_biomass_density_model_2019_density_exp.rda")
save(evSSeedDensMod, file = "output/evS_plot_seed_density_model_2019_density_exp.rda")
save(evABioDensMod, file = "output/evA_plot_biomass_density_model_2019_density_exp.rda")
save(evASeedDensMod, file = "output/evA_plot_seed_density_model_2019_density_exp.rda")

# load models
load("output/evS_plot_biomass_density_model_2019_density_exp.rda")
load("output/evS_plot_seed_density_model_2019_density_exp.rda")
load("output/evA_plot_biomass_density_model_2019_density_exp.rda")
load("output/evA_plot_seed_density_model_2019_density_exp.rda")

# output tables
write_csv(tidy(evSBioDensMod), "output/evS_plot_biomass_density_model_2019_dens_exp.csv")
write_csv(tidy(evSSeedDensMod), "output/evS_plot_seed_density_model_2019_dens_exp.csv")
write_csv(tidy(evABioDensMod), "output/evA_plot_biomass_density_model_2019_dens_exp.csv")
write_csv(tidy(evASeedDensMod), "output/evA_plot_seed_density_model_2019_dens_exp.csv")


#### predicted values ####

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}

# prediction dataset
plotPredDatTemplate <- plotD2Dat %>%
  group_by(treatment, fungicide, age) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

evSSeedPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "seedling") %>%
  mutate(value = fitted(evSSeedDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evSSeedDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evSSeedDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evASeedPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "adult") %>%
  mutate(value = fitted(evASeedDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evASeedDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evASeedDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evSBioPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "seedling") %>%
  mutate(value = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evABioPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "adult") %>%
  mutate(value = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evPredD2Dat <- evSSeedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(evASeedPredD2Dat %>%
              mutate(response = "seeds")) %>%
  full_join(evSBioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  full_join(evABioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  mutate(plant_group = if_else(age == "adult", "Adult competitor (Ev)",
                               "1st yr competitor (Ev)") %>%
           fct_relevel("1st yr competitor (Ev)"))


#### figure ####

# theme
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

# colors
col_pal = c("black", "#238A8D")

# dodge size
dodge_width <- 0.5

# labels
plot_labels <- c(biomass = "Biomass~(g~m^-2)",
                 seeds = "Seed~production~(m^-2)")

# figure
pdf("output/ev_density_figure_2019_density_exp.pdf", width = 3.54, height = 3.54)
ggplot(evPredD2Dat, aes(x = density, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment)) +
  geom_point(data = figD2Dat, aes(color = treatment), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  facet_grid(rows = vars(response),
             cols = vars(plant_group),
             scales = "free",
             switch = "y",
             labeller = labeller(response = as_labeller(plot_labels, label_parsed))) +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  labs(x = expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(axis.title.y = element_blank())
dev.off()


#### treatment effects ####

# treatment effects
b0_eff = "b0_treatmentfungicide - b0_treatmentcontrol = 0"
alpha_eff = "alpha_treatmentfungicide - alpha_treatmentcontrol = 0"

evSBioEff <- hypothesis(evSBioDensMod, c(b0_eff, alpha_eff))
evSSeedEff <- hypothesis(evSSeedDensMod, c(b0_eff, alpha_eff))
evABioEff <- hypothesis(evABioDensMod, c(b0_eff, alpha_eff))
evASeedEff <- hypothesis(evASeedDensMod, c(b0_eff, alpha_eff))

# combine treatment effects
evSBioEff[[1]] %>%
  mutate(response = "biomass",
         age = "seedling") %>%
  full_join(evSSeedEff[[1]] %>%
              mutate(response = "seeds",
                     age = "seedling")) %>%
  full_join(evABioEff[[1]] %>%
              mutate(response = "biomass",
                     age = "adult")) %>%
  full_join(evASeedEff[[1]] %>%
              mutate(response = "seeds",
                     age = "adult")) %>%
  mutate(parameter = rep(c("b0", "alpha"), 4)) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob)) %>%
  relocate(age, response, parameter)




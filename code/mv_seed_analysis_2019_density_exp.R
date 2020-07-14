##### info ####

# file: mv_seed_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/11/20
# goal: analyze Mv seeds


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)
library(glmmTMB)

# import data
plant_dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
plot_dat <- read_csv("intermediate-data/mv_plot_level_seeds_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sev19 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

# severity
sev19b <- sev19 %>%
  mutate(plant_severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, plot, treatment, sp, ID, plant_severity) %>%
  pivot_wider(names_from = "month", values_from = "plant_severity", names_glue = "severity_{month}") %>%
  filter(sp == "Mv") %>%
  mutate(plant = as.numeric(ID)) %>%
  select(-ID)

# add density to data for visualization
plot_dat2 <- left_join(plot_dat, plots)

# data for stats
plant_dat2 <- left_join(plant_dat, treat) %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0") %>% as.numeric(),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide"),
         exp_plot = paste(plot, toupper(substr(treatment, 1, 1)), sep = ""),
         mv_density = case_when(plot %in% c(1, 5:10) ~ 3,
                                TRUE ~ background_density + 3),
         ev_seedling_density = case_when(plot %in% c(5:7) ~ background_density + 3,
                                         TRUE ~ 3),
         ev_adult_density = case_when(plot %in% c(8:10) ~ background_density + 1,
                                         TRUE ~ 1)) %>%
  left_join(sev19b)

# no background data
no_back_plant_dat <- filter(plant_dat2, plot == 1)

# individual-level seeds and plots
plant_dat3 <- left_join(plant_dat, plots) %>%
  mutate(fungicide = recode(treatment, fungicide = "1", water = "0") %>% as.numeric(),
         Treatment = recode(treatment, water = "control (water)", fungicide = "fungicide"))


#### visualize ####

# template theme
temp_theme <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside")

# treatment figure
seed_viz <- ggplot(plot_dat2, aes(x = background_density, y = seeds, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  xlab("Background density") +
  temp_theme
seed_viz
# seeds drop with Microstegium density, but are pretty unaffected by Elymus density

# flower seeds
seed_viz %+%
  aes(y = flower_seeds)
# similar pattern

# stem seeds
seed_viz %+%
  aes(y = stem_seeds)
# similar pattern

# individual seeds
ggplot(mv_dat, aes(x = mv_density, y = seeds, fill = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3)) +
  scale_fill_manual(values = c("black", "white"), name = "Treatment") +
  temp_theme

# no background
mv_dat %>%
  filter(plot == 1) %>%
  ggplot(aes(x = Treatment, y = seeds, fill = treatment)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  temp_theme


#### direct disease models ####

# control seeds
filter(no_back_plant_dat, fungicide == 0) %>%
  summarise(seeds = mean(seeds))

# model
mv_no_back_seed_mod <- brm(data = no_back_plant_dat, family = gaussian,
                          seeds ~ fungicide + (1|site),
                          prior <- c(prior(normal(1200, 100), class = Intercept),
                                     prior(normal(0, 100), class = b),
                                     prior(cauchy(0, 1), class = sd),
                                     prior(cauchy(0, 1), class = sigma)),
                          iter = 6000, warmup = 1000, chains = 1,
                          control = list(adapt_delta = 0.99))

# check model and add chains
summary(mv_no_back_seed_mod)
prior_summary(mv_no_back_seed_mod)
mv_no_back_seed_mod <- update(mv_no_back_seed_mod, chains = 3)
plot(mv_no_back_seed_mod)
pp_check(mv_no_back_seed_mod, nsamples = 100)

# reduction due to fungicide
-0.12*1182

# model with direct fungicide effects
mv_no_back_seed_fung_mod <- brm(data = no_back_plant_dat, family = gaussian,
                               seeds ~ Treatment + fungicide + (1|site),
                               prior <- c(prior(normal(1182, 100), class = Intercept),
                                          prior(normal(0, 100), class = b),
                                          prior(normal(-142, 0.001), class = b, coef = "Treatmentfungicide"),
                                          prior(cauchy(0, 1), class = sd),
                                          prior(cauchy(0, 1), class = sigma)),
                               iter = 6000, warmup = 1000, chains = 3,
                               control = list(adapt_delta = 0.99))

# check model
summary(mv_no_back_seed_fung_mod)
prior_summary(mv_no_back_seed_fung_mod)
plot(mv_no_back_seed_fung_mod)
pp_check(mv_no_back_seed_fung_mod, nsamples = 100)


#### disease and density models ####

# split by background
mv_mv_seed_dat <- plant_dat2 %>%
  filter(background == "Mv seedling" & !is.na(seeds))

mv_evs_seed_dat <- plant_dat2 %>%
  filter(background == "Ev seedling")

mv_eva_seed_dat <- plant_dat2 %>%
  filter(background == "Ev adult" & !is.na(seeds))

# check for complete values
sum(is.na(mv_mv_seed_dat$seeds))
sum(is.na(mv_evs_seed_dat$seeds))
sum(is.na(mv_eva_seed_dat$seeds))

# Beverton-Holt function
bh_fun <- function(dat_in, a, y0){
  
  # extract values
  xmin = 0
  xmax = max(dat_in$background_density)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean", aes(color = treatment)) +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, aes(color = treatment)) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# try values of a
mean(filter(mv_mv_seed_dat, density_level == "low")$seeds)
bh_fun(mv_mv_seed_dat, 0.05, 1300)
mean(filter(mv_evs_seed_dat, density_level == "low")$seeds)
bh_fun(mv_evs_seed_dat, -0.01, 1000)
mean(filter(mv_eva_seed_dat, density_level == "low")$seeds)
bh_fun(mv_eva_seed_dat, 0.1, 2000)

# distributions
x = 0:10
y = dexp(x, 0.5)
plot(x, y, type = "l")

# Mv background model
mv_mv_seed_mod <- brm(data = mv_mv_seed_dat, family = gaussian,
                        bf(seeds ~ y0 / (1 + alpha * background_density),
                           y0 ~ 0 + Treatment + (1|site),
                           alpha ~ 0 + Treatment,
                           nl = T),
                      prior <- c(prior(normal(1300, 100), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_mv_seed_mod)
prior_summary(mv_mv_seed_mod)
mv_mv_seed_mod <- update(mv_mv_seed_mod, chains = 3)
plot(mv_mv_seed_mod)
pp_check(mv_mv_seed_mod, nsamples = 100)    

# Ev seedling background model
mv_evs_seed_mod <- brm(data = mv_evs_seed_dat, family = gaussian,
                      bf(seeds ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(1000, 100), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 1), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 3,
                      control = list(adapt_delta = 0.99))

# check model and increase chains
summary(mv_evs_seed_mod)
prior_summary(mv_evs_seed_mod)
plot(mv_evs_seed_mod)
pp_check(mv_evs_seed_mod, nsamples = 100)    

# Ev adult background model
mv_eva_seed_mod <- brm(data = mv_eva_seed_dat, family = gaussian,
                       bf(seeds ~ y0 / (1 + alpha * background_density),
                          y0 ~ 0 + Treatment + (1|site),
                          alpha ~ 0 + Treatment,
                          nl = T),
                       prior <- c(prior(normal(2000, 100), nlpar = "y0", class = "b", lb = 0),
                                  prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                  prior(cauchy(0, 1), nlpar = "y0", class = "sd"),
                                  prior(cauchy(0, 1), class = "sigma")),
                       iter = 8000, warmup = 3000, chains = 1,
                       control = list(adapt_delta = 0.99999))

# check model and increase chains
summary(mv_eva_seed_mod)
mv_eva_seed_mod <- update(mv_eva_seed_mod, chains = 3)
mv_eva_seed_mod <- update(mv_eva_seed_mod, iter = 10000, warmup = 5000)
# issue with bulk ESS despite higher chain iterations and r-hat for site is over 1
# adjust priors for sd
mv_eva_seed_mod <-brm(data = mv_eva_seed_dat, family = gaussian,
                      bf(seeds ~ y0 / (1 + alpha * background_density),
                         y0 ~ 0 + Treatment + (1|site),
                         alpha ~ 0 + Treatment,
                         nl = T),
                      prior <- c(prior(normal(2000, 100), nlpar = "y0", class = "b", lb = 0),
                                 prior(exponential(0.5), nlpar = "alpha", lb = 0),
                                 prior(cauchy(0, 10), nlpar = "y0", class = "sd"),
                                 prior(cauchy(0, 10), class = "sigma")),
                      iter = 6000, warmup = 1000, chains = 1,
                      control = list(adapt_delta = 0.99))
# can accept R-hat values less than 1.05 with 4 chains: https://mc-stan.org/rstan/reference/Rhat.html
mv_eva_seed_mod <- update(mv_eva_seed_mod, chains = 4)
plot(mv_eva_seed_mod)
# site still has two peaks
pp_check(mv_eva_seed_mod, nsamples = 100)   


#### model fit figure ####

# model fits over simulated data
seed_mv_sim_dat <- tibble(background_density = rep(seq(0, 64, length.out = 300), 2),
                         fungicide = rep(c(0, 1), each = 300)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_mv_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seed_lower = fitted(mv_mv_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seed_upper = fitted(mv_mv_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"])

seed_evs_sim_dat <- tibble(background_density = rep(seq(0, 16, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_evs_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seed_lower = fitted(mv_evs_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seed_upper = fitted(mv_evs_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"])

seed_eva_sim_dat <- tibble(background_density = rep(seq(0, 8, length.out = 100), 2),
                          fungicide = rep(c(0, 1), each = 100)) %>%
  mutate(site = NA,
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide")) %>%
  mutate(seeds = fitted(mv_eva_seed_mod, newdata = ., re_formula = NA)[, "Estimate"],
         seed_lower = fitted(mv_eva_seed_mod, newdata = ., re_formula = NA)[, "Q2.5"],
         seed_upper = fitted(mv_eva_seed_mod, newdata = ., re_formula = NA)[, "Q97.5"])

# combine simulated data
seed_sim_dat <- seed_mv_sim_dat %>%
  mutate(background = "Mv seedling") %>%
  full_join(seed_evs_sim_dat %>%
              mutate(background = "Ev seedling")) %>%
  full_join(seed_eva_sim_dat %>%
              mutate(background = "Ev adult"))

# colors
col_pal = c("#a6611a", "#018571")

# figure
pdf("output/mv_seed_analysis_model_fits_2019_density_exp.pdf")
ggplot(plant_dat3, aes(x = background_density, y = seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), aes(group = Treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 3, shape = 21, position = position_dodge(0.3), aes(fill = Treatment)) +
  geom_ribbon(data = seed_sim_dat, alpha = 0.5, aes(ymin = seed_lower, ymax = seed_upper, fill = Treatment)) +
  geom_line(data = seed_sim_dat, aes(color = Treatment)) +
  facet_wrap(~ background, scales = "free_x", strip.position = "bottom") +
  xlab("Background density") +
  ylab(expression(paste(italic(Microstegium), " seedmass (g ", m^-1, ")", sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  temp_theme
dev.off()


#### severity and biomass ####

# water treatment
mv_sev_dat <- filter(plant_dat2, treatment == "water") %>%
  as.data.frame()

# check for density-severity correlations
for(i in 33:37){
  mv_sev_dat$severity = mv_sev_dat[,i]
  print(cor.test(mv_sev_dat$background_density, mv_sev_dat$severity))
}
# none significantly correlated

# seed-severity relationships
mv_seed_sev_jun_19_mod <- glmmTMB(seeds ~ severity_jun + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_seed_sev_jun_19_mod) # not sig
mv_seed_sev_jul_19_mod <- glmmTMB(seeds ~ severity_jul + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_seed_sev_jul_19_mod) # not sig
mv_seed_sev_eau_19_mod <- glmmTMB(seeds ~ severity_early_aug + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_seed_sev_eau_19_mod) # not sig
mv_seed_sev_lau_19_mod <- glmmTMB(seeds ~ severity_late_aug + (1|site/plot), family = gaussian, data = mv_sev_dat)
summary(mv_seed_sev_lau_19_mod) # not sig


#### output ####
save(mv_no_back_seed_mod, file = "output/mv_seeds_no_background_model_2019_density_exp.rda")
save(mv_no_back_seed_fung_mod, file = "output/mv_seeds_no_background_model_greenhouse_fungicide_2019_density_exp.rda")
save(mv_mv_seed_mod, file = "output/mv_seeds_mv_background_model_2019_density_exp.rda")
save(mv_evs_seed_mod, file = "output/mv_seeds_ev_seedling_background_model_2019_density_exp.rda")
save(mv_eva_seed_mod, file = "output/mv_seeds_ev_adult_background_model_2019_density_exp.rda")

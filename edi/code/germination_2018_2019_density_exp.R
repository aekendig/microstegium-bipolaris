##### outputs ####

# mv_germination_fungicide_model_2018_density_exp.rda
# mv_germination_fungicide_model_data_2018_density_exp.csv
# ev_germination_fungicide_model_2018_2019_density_exp.rda
# ev_germination_fungicide_model_data_2018_2019_density_exp.rda
# mv_germination_infection_model_2018_density_exp.csv (Table S18)
# mv_seed_infection_dark_fungicide_model_2018_density_exp.csv (Table S17)
# mv_germination_infection_figure_2018_density_exp.pdf (Fig. S5)


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)
library(car)
library(tidybayes)
library(broom.mixed)

# import data
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}

# convert logit to probability
logit2prob <- function(x){
  exp(x)/(1 + exp(x))
}


#### edit data ####

# notes
unique(mvGermD1Dat1$notes_check_1) # seedlings with lesions in notes
unique(mvGermD1Dat1$notes_check_2) # may be a contaminate on plate
unique(mvGermD1Dat1$notes_check_3)
unique(mvGermD1Dat1$notes_germination_check_1) # some plates were put into fridge during one day
unique(mvGermD1Dat1$notes_germination_final) 
unique(mvGermD1Dat2$notes) # may be a contaminate on plate

# Mv data
# average across trials
mvGermD1Dat <- mvGermD1Dat1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final),
         seeds_dark_check_1 = rowSums(cbind(seeds_dark_check_1, seeds_pink_check_1, seeds_red_check_1, seeds_green_check_1), na.rm = T),
         seeds_dark_check_2 = rowSums(cbind(seeds_dark_check_2, seeds_pink_check_2, seeds_red_check_2, seeds_green_check_2), na.rm = T),
         seeds_seeds_dark_check_3 = rowSums(cbind(seeds_dark_check_3, seeds_red_check_3, seeds_green_check_3), na.rm = T),
         seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, seeds_dark_check_3, na.rm = T),
         seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, seeds_light_check_3, na.rm = T)) %>%
  select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light) %>%
  full_join(mvGermD1Dat2 %>%
              mutate(seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, na.rm = T),
                     seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, na.rm = T)) %>%
              select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light)) %>%
  mutate(seeds_infect = case_when(seeds_dark > 0 | seeds_light > 0 ~ 1,
                                  TRUE ~ 0), # a single seed can be infected with dark or light
         prop_germ = germination_final / seeds,
         prop_dark = seeds_dark / seeds,
         prop_light = seeds_light / seeds,
         site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           dplyr::recode("F" = "fungicide", "W" = "control") %>%
           fct_rev(),
         site = ifelse(site == "P1", "D1", site),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)))

# infection variable
unique(mvGermD1Dat$seeds_infect) # all samples have some infection

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 ) %>%
  data.frame()

# trials per plot
mvGermD1Dat %>%
  count(site_plot) %>%
  rename(trials = "n") %>%
  count(trials)

# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
evGermDat2 <- evGermDat %>%
  filter(seeds_planted > 0) %>%
  mutate(week_4_emerg = case_when(week_4_emerg < week_3_emerg ~ week_3_emerg,
                                  TRUE ~ week_4_emerg),
         week_3_cut_tops = case_when(week_4_cut_tops > week_3_cut_tops & week_4_cut_tops <= week_2_emerg ~ week_4_cut_tops,
                                     TRUE ~ week_3_cut_tops),
         week_4_cut_tops = case_when(week_4_cut_tops > week_2_emerg ~ week_3_cut_tops,
                                     week_4_cut_tops < week_3_cut_tops ~ week_3_cut_tops,
                                     TRUE ~ week_4_cut_tops),
         week_3_new_emerg = week_3_emerg - week_3_cut_tops,
         week_4_new_emerg = week_4_emerg - week_4_cut_tops,
         germinants = week_2_emerg + week_4_new_emerg + week_4_soil_germ,
         germination = germinants/seeds_planted,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         yearf = ifelse(year == 2018, 0, 1),
         treatment = fct_recode(treatment, "control" = "water") %>%
           fct_rev())
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()
# 40 in 2018
# 56 in 2019

# split by year
evGermD1Dat <- evGermDat2 %>% filter(year == 2018)
evGermD2Dat <- evGermDat2 %>% filter(year == 2019)


#### Mv models ####

# initial visualization
mvGermD1Dat %>%
  select(prop_dark, prop_light, prop_germ) %>%
  ggpairs()

ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# light/dark infection correlation
cor.test(~ prop_dark + prop_light, data = mvGermD1Dat) # not correlated

# model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ prop_dark + prop_light + (1|plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermD1Mod)
# prop dark decreases germination
# prop light increases germination

mvGermD1Mod2 <- brm(data = mvGermD1Dat, family = binomial,
                    germination_final | trials(seeds) ~ fungicide + (1|plotf),
                    prior <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3)
summary(mvGermD1Mod2)
# fungicide doesn't affect germination

mvPropDarkMod2 <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ fungicide + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod2)
# fungicide decreases prop dark

mvPropLightMod2 <- brm(data = mvGermD1Dat, family = binomial,
                      seeds_light | trials(seeds) ~ fungicide + (1|plotf),
                      prior <- c(prior(normal(0, 10), class = "Intercept"),
                                 prior(normal(0, 10), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropLightMod2)
# fungicide doesn't affect prop light

# save
save(mvGermD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")
save(mvGermD1Mod2, file = "output/mv_germination_fungicide_model_2018_density_exp.rda")
save(mvPropDarkMod2, file = "output/mv_seed_infection_dark_fungicide_model_2018_density_exp.rda")
save(mvPropLightMod2, file = "output/mv_seed_infection_light_fungicide_model_2018_density_exp.rda")

# load
load("output/mv_germination_infection_model_2018_density_exp.rda")
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
load("output/mv_seed_infection_dark_fungicide_model_2018_density_exp.rda")
load("output/mv_seed_infection_light_fungicide_model_2018_density_exp.rda")

# tables
write_csv(tidy(mvGermD1Mod), "output/mv_germination_infection_model_2018_density_exp.csv")
write_csv(tidy(mvPropDarkMod2), "output/mv_seed_infection_dark_fungicide_model_2018_density_exp.csv")

# save corresponding data
write_csv(mvGermD1Dat, "output/mv_germination_fungicide_model_data_2018_density_exp.csv")


#### Ev models ####

# initial visualization
ggplot(evGermDat2, aes(year, germination, color = age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(0.2))

ggplot(evGermDat2, aes(treatment, germination)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                 germinants | trials(seeds_planted) ~ fungicide + (1|site) + (1|yearf),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, 
                 control = list(adapt_delta = 0.999, max_treedepth = 15))
mod_check_fun(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# save corresponding data
write_csv(evGermDat2, "output/ev_germination_fungicide_model_data_2018_2019_density_exp.rda")


#### figure ####

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.box = "vertical",
        plot.title = element_text(size = 7, hjust = 0.5))

# simulated data
mvFungSim <- tibble(fungicide = c(0, 1)) %>%
  mutate(plotf = "A",
         treatment = c("control", "fungicide"),
         seeds = round(mean(mvGermD1Dat$seeds))) %>%
  mutate(prop_dark = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

mvGermSim <- tibble(prop_dark = seq(min(mvGermD1Dat$prop_dark), max(mvGermD1Dat$prop_dark), length.out = 100)) %>%
  mutate(plotf = "A",
         seeds = round(mean(mvGermD1Dat$seeds)),
         prop_light = mean(mvGermD1Dat$prop_light)) %>%
  mutate(prop_germ = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

# figures
mvFungFig <- ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
    # geom_point(alpha = 0.3, size = 0.7, position = position_jitter(width = 0.35)) +
  # geom_violin() +
  ggdist::stat_dots(dotsize = .4, 
                    binwidth = 0.01,
                    color = "#238A8DFF",
                    fill = "#238A8DFF",
                    alpha = 0.7,
                    side = "both") +
    geom_errorbar(data = mvFungSim, width = 0, size = 0.25, aes(ymin = lower, ymax = upper)) +
    geom_point(data = mvFungSim, size = 1.25) +
  fig_theme +
  theme(axis.title.x = element_blank()) +
  ylab("Infected fraction")

mvGermFig <- ggplot(mvGermD1Dat, aes(prop_dark, prop_germ)) +
  geom_point(alpha = 0.3, size = 0.7, color = "#238A8DFF") +
  geom_line(data = mvGermSim) +
  geom_ribbon(data = mvGermSim, aes(ymin = lower, ymax = upper), alpha = 0.4) +
  fig_theme +
  xlab("Infected fraction") +
  ylab("Germination fraction")

pdf("output/mv_germination_infection_figure_2018_density_exp.pdf", width = 3.54, height = 1.57)
plot_grid(mvFungFig, mvGermFig,
          nrow = 1,
          labels = LETTERS[1:2],
          label_size = 10)
dev.off()


#### values for text ####

# Mv fungicide effect on infection
hypothesis(mvPropDarkMod2, "logit2prob(Intercept + fungicide) - logit2prob(Intercept) = 0")

# Mv infection effect on germination
filter(mvGermSim, prop_dark == min(prop_dark))
filter(mvGermSim, prop_dark == max(prop_dark))

# Mv infection effect on germination
set.seed(184)
posterior_predict(mvGermD1Mod,
                  newdata = filter(mvGermD1Dat, 
                                   prop_dark %in% c(max(mvGermD1Dat$prop_dark), min(mvGermD1Dat$prop_dark))) %>%
                    select(prop_dark) %>%
                    unique() %>%
                    mutate(prop_light = 0,
                           seeds = 30,
                           plotf = "A"),
                  allow_new_levels = T) %>%
  as_tibble(.name_repair = ~ c("min_inf", "max_inf")) %>% # min_inf is first row
  mutate(min_inf_prop = min_inf / 30,
         max_inf_prop = max_inf / 30,
         inf_eff = 100 * (max_inf_prop - min_inf_prop) / min_inf_prop) %>%
  select(min_inf_prop, max_inf_prop, inf_eff) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)

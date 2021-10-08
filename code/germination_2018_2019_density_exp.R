##### info ####

# file: germination_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/5/21
# goal: analyses of germination


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)
library(car)

# import data
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")
sevD1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R

# model functions
source("code/brms_model_fitting_functions.R")

# convert logit to probability
logit2prob <- function(x){
  exp(x)/(1 + exp(x))
}


#### edit data ####

# severity
sevD1Dat2 <- sevD1Dat %>%
  select(month, site, treatment, plot, sp, age, severity) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity") %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

sevD2Dat2 <- sevD2Dat %>%
  select(month, site, treatment, plot, sp, age, severity) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity") %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

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
           dplyr::recode("F" = "fungicide", "W" = "control (water)") %>%
           fct_rev(),
         site = ifelse(site == "P1", "D1", site),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  left_join(sevD1Dat2 %>%
              filter(sp == "Mv"))

# infection variable
unique(mvGermD1Dat$seeds_infect) # all samples have some infection

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 | prop_infect > 1) %>%
  data.frame()
# don't use prop_infect
filter(mvGermD1Dat, is.na(jul_severity)) # 1 plot
filter(mvGermD1Dat, is.na(late_aug_severity)) # 3 plots
filter(mvGermD1Dat, is.na(sep_severity)) # 3 plots (one missing late aug)

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
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev()) %>%
  left_join(sevD1Dat2 %>%
              filter(sp == "Ev") %>%
              mutate(year = 2018) %>%
              full_join(sevD2Dat2 %>%
                          filter(sp == "Ev") %>%
                          mutate(year = 2019))) %>%
  mutate(final_severity = case_when(year == 2018 ~ sep_severity,
                                    year == 2019 ~ late_aug_severity))
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()
# 40 in 2018
# 56 in 2019

# number missing?
filter(evGermDat2, is.na(jul_severity)) %>% data.frame() # 3 all 2019
filter(evGermDat2, is.na(late_aug_severity)) %>% data.frame() # 27
filter(evGermDat2, is.na(final_severity)) %>% data.frame() # 31
filter(evGermDat2, year == 2018 & is.na(late_aug_severity)) %>% data.frame() # 12
filter(evGermDat2, year == 2018 & is.na(sep_severity)) %>% data.frame() # 16
filter(evGermDat2, year == 2019 & is.na(early_aug_severity)) %>% data.frame() # 7
filter(evGermDat2, year == 2019 & is.na(late_aug_severity)) %>% data.frame() # 15

# split by year
evGermD1Dat <- evGermDat2 %>% filter(year == 2018)
evGermD2Dat <- evGermDat2 %>% filter(year == 2019)


#### Mv models ####

# initial visualization
mvGermD1Dat %>%
  select(jul_severity, late_aug_severity, sep_severity,
         prop_dark, prop_light, prop_germ) %>%
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
                   germination_final | trials(seeds) ~ sep_severity + (1|plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
summary(mvGermD1Mod2)
# no effect of severity

mvGermD1Mod3 <- brm(data = mvGermD1Dat, family = binomial,
                    germination_final | trials(seeds) ~ fungicide + (1|plotf),
                    prior <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3)
summary(mvGermD1Mod3)
# fungicide doesn't affect germination

mvPropDarkMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ sep_severity + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod)
# severity increases prop dark

mvPropDarkMod2 <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ fungicide + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod2)
# fungicide decreases prop dark

mvPropLightMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_light | trials(seeds) ~ sep_severity + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropLightMod)
# severity increases prop light

mvPropLightMod2 <- brm(data = mvGermD1Dat, family = binomial,
                      seeds_light | trials(seeds) ~ fungicide + (1|plotf),
                      prior <- c(prior(normal(0, 10), class = "Intercept"),
                                 prior(normal(0, 10), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropLightMod2)
# fungicide doesn't affect prop light

# save
save(mvGermD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")
save(mvPropDarkMod, file = "output/mv_seed_infection_dark_model_2018_density_exp.rda")
save(mvPropDarkMod2, file = "output/mv_seed_infection_dark_fungicide_model_2018_density_exp.rda")
save(mvPropLightMod, file = "output/mv_seed_infection_light_model_2018_density_exp.rda")
save(mvPropLightMod2, file = "output/mv_seed_infection_light_fungicide_model_2018_density_exp.rda")


#### Ev model ####

# initial visualization
evGermD1Dat %>%
  select(jul_severity, late_aug_severity, sep_severity, germination) %>%
  ggpairs()

evGermD2Dat %>%
  select(jul_severity, early_aug_severity, late_aug_severity, germination) %>%
  ggpairs()

evGermDat2 %>%
  select(late_aug_severity, final_severity, germination) %>%
  ggpairs()

ggplot(evGermDat2, aes(year, germination, color = age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(0.2))

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                   germinants | trials(seeds_planted) ~ final_severity*yearf + (1|site),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, 
                 control = list(adapt_delta = 0.999, max_treedepth = 15))
mod_check_fun(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_severity_model_2018_2019_density_exp.rda")


#### figure ####

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

# simulated data
mvFungSim <- tibble(fungicide = c(0, 1)) %>%
  mutate(plotf = "A",
         treatment = c("control (water)", "fungicide"),
         seeds = round(mean(mvGermD1Dat$seeds))) %>%
  mutate(prop_dark = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

mvPropDarkSim <- tibble(sep_severity = seq(min(mvGermD1Dat$sep_severity, na.rm = T), max(mvGermD1Dat$sep_severity, na.rm = T), length.out = 100)) %>%
  mutate(plotf = "A",
         seeds = round(mean(mvGermD1Dat$seeds))) %>%
  mutate(prop_dark = fitted(mvPropDarkMod, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvPropDarkMod, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvPropDarkMod, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

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
    geom_errorbar(data = mvFungSim, width = 0, size = 0.75, aes(ymin = lower, ymax = upper)) +
    geom_point(data = mvFungSim, size = 2.5) +
  geom_text(x = 1.5, y = 0.75, label = "*", size = 6, check_overlap = T) +
  fig_theme +
  xlab("Plot treatment") +
  ylab("Proportion seeds infected")

mvDarkFig <- ggplot(mvGermD1Dat, aes(sep_severity, prop_dark)) +
  geom_point(alpha = 0.3, size = 0.7, color = "#238A8DFF") +
  geom_line(data = mvPropDarkSim) +
  geom_ribbon(data = mvPropDarkSim, aes(ymin = lower, ymax = upper), alpha = 0.4) +
  fig_theme +
  xlab("Leaf disease severity") +
  ylab("Proportion seeds infected")

mvGermFig <- ggplot(mvGermD1Dat, aes(prop_dark, prop_germ)) +
  geom_point(alpha = 0.3, size = 0.7, color = "#238A8DFF") +
  geom_line(data = mvGermSim) +
  geom_ribbon(data = mvGermSim, aes(ymin = lower, ymax = upper), alpha = 0.4) +
  fig_theme +
  xlab("Proportion seeds infected") +
  ylab("Proportion seeds germinated") +
  theme(axis.title.y = element_text(size = 10, hjust = -0.05))

tiff("output/mv_germination_infection_figure_2018_density_exp.tiff", width = 18, height = 6, units = "cm", res = 300, compression = "lzw")
plot_grid(mvFungFig, mvDarkFig, mvGermFig,
          nrow = 1,
          labels = LETTERS[1:3],
          vjust = 1.1)
dev.off()


#### values for text ####

# Mv fungicide effect on infection
hypothesis(mvPropDarkMod2, "logit2prob(Intercept + fungicide) - logit2prob(Intercept) = 0")

# Mv severity effect on infection
filter(mvPropDarkSim, sep_severity == min(sep_severity))
filter(mvPropDarkSim, sep_severity == max(sep_severity))

# Mv infection effect on germination
filter(mvGermSim, prop_dark == min(prop_dark))
filter(mvGermSim, prop_dark == max(prop_dark))

# Ev severity hypotheses
ev18 = "final_severity = 0"
ev19 = "final_severity + final_severity:yearf = 0"
hypothesis(evGermMod, c(ev18, ev19))

# Mv model table
mod_table_fun <- function(mod){
  fixef(mod) %>%
    as_tibble() %>%
    mutate(Predictor = rownames(fixef(mod)),
           N = summary(mod)$nobs,
           SD = summary(mod)$random$plotf$Estimate) %>%
    filter(Predictor != "Intercept")
  }

mv_mod_table <- bind_rows(mod_table_fun(mvPropDarkMod2) %>%
                            mutate(Response = "prop. infected (dark)"),
                          mod_table_fun(mvPropLightMod2) %>%
                            mutate(Response = "prop. infected (light)"),
                          mod_table_fun(mvPropDarkMod) %>%
                            mutate(Response = "prop. infected (dark)"),
                          mod_table_fun(mvPropLightMod) %>%
                            mutate(Response = "prop. infected (light)"),
                          mod_table_fun(mvGermD1Mod) %>%
                            mutate(Response = "prop. germinated")
) %>%
  relocate(Response, Predictor)

write_csv(mv_mod_table, "output/mv_germination_infection_table_2018_density_exp.csv")

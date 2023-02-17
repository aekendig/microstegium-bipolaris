##### outputs ####

# survival_fungicide_model_2019_density_exp.rda
# survival_fungicide_model_data_2019_density_exp.csv
# ev_adult_survival_fungicide_model_2018_2019_density_exp.rda


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}


#### edit data ####

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  left_join(plots) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  select(-c(month, field_notes, seeds_produced)) %>%
  filter(!is.na(survival))

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
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

# 2019 focal survival
# 2018 survival starts in June because plants were replaced through May 24
survD2Dat2 <- survD2Dat %>%
  filter(focal == 1 & replace_date > 20190531) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(plantings = length(unique(replace_date)) + 1) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plots) %>%
  mutate(plantings = replace_na(plantings, 1),
         survival = ifelse(plantings > 1, 0, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""))

# winter survival 2018-2019
winSurvD1Dat <- survD1Dat %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  left_join(plots) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         focal = paste(sp, age, sep = " ") %>%
           fct_recode(Mv = "Mv seedling"),
         foc = fct_recode(focal, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         bg = fct_recode(background, m = "Mv", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April)

# annual adult survival 2018-2019
adultSurvD1Dat <- survD1Dat %>%
  filter(month == "April" & age == "adult" & !is.na(survival)) %>%
  select(-field_notes) %>%
  left_join(plots) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0))
# includes non-focal


#### fit models ####

survFungD1Mod <- brm(data = survD1Dat2, family = bernoulli,
                    survival ~ fungicide * foc + (1|plotf),
                    prior <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3) 
mod_check_fun(survFungD1Mod)

survFungD2Mod <- brm(data = survD2Dat2, family = bernoulli,
                    survival ~ fungicide * foc + (1|plotf),
                    prior <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(survFungD2Mod)

winSurvFungD1Mod <- brm(data = winSurvD1Dat, family = bernoulli,
                       survival ~ fungicide * foc + (1|plotf),
                       prior <- c(prior(normal(0, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")), # use default for sigma
                       iter = 6000, warmup = 1000, chains = 3, cores = 3,
                       control = list(adapt_delta = 0.99)) 
mod_check_fun(winSurvFungD1Mod)

# save models
save(survFungD1Mod, file = "output/survival_fungicide_model_2018_density_exp.rda")
save(survFungD2Mod, file = "output/survival_fungicide_model_2019_density_exp.rda")
save(winSurvFungD1Mod, file = "output/winter_survival_fungicide_model_2018_density_exp.rda")

# load
load("output/survival_fungicide_model_2018_density_exp.rda")
load("output/survival_fungicide_model_2019_density_exp.rda")
load("output/winter_survival_fungicide_model_2018_density_exp.rda")

# save corresponding data
write_csv(survD2Dat2, "output/survival_fungicide_model_data_2019_density_exp.csv")


#### adult survival ####

mean(adultSurvD1Dat$survival)

adultSurvD1Mod <- brm(data = adultSurvD1Dat, family = bernoulli,
                      survival ~ 1 + (1|site),
                      prior <- prior(normal(0, 1), class = "Intercept"), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.99)) 
mod_check_fun(adultSurvD1Mod)

adultSurvD1Mod2 <- brm(data = adultSurvD1Dat, family = bernoulli,
                      survival ~ fungicide + (1|site),
                      prior <- c(prior(normal(0, 10), class = "Intercept"), # use default for sigma
                                 prior(normal(0, 10), class = "b")),
                      iter = 6000, warmup = 1000, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.999)) 
mod_check_fun(adultSurvD1Mod2)

# save
save(adultSurvD1Mod, file = "output/ev_adult_survival_model_2018_2019_density_exp.rda")
save(adultSurvD1Mod2, file = "output/ev_adult_survival_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_adult_survival_model_2018_2019_density_exp.rda")
load("output/ev_adult_survival_fungicide_model_2018_2019_density_exp.rda")


#### fungicide effects ####

mv_fung_eff = "fungicide = 0"
evS_fung_eff = "fungicide + fungicide:focs = 0"
evA_fung_eff = "fungicide + fungicide:foca = 0"

hypothesis(survFungD1Mod,
                          c(mv_fung_eff, evS_fung_eff, evA_fung_eff)) [[1]] %>%
  mutate(Year = "2018", Season = "growing season") %>%
  full_join(hypothesis(survFungD2Mod,
                       c(mv_fung_eff, evS_fung_eff, evA_fung_eff)) [[1]] %>%
              mutate(Year = "2019", Season = "growing season")) %>%
  full_join(hypothesis(winSurvFungD1Mod,
                       c(mv_fung_eff, evS_fung_eff)) [[1]] %>%
              mutate(Year = "2018-2019", Season = "winter")) %>%
  mutate(Focal = c(rep(c("Mv", "Ev first-year", "Ev adult"), 2), "Ev adult", "Ev first-year")) %>%
  select(Year, Season, Focal, Estimate:CI.Upper) %>%
  arrange(Season, Year, Focal)


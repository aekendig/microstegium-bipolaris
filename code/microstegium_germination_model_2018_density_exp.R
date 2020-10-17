##### info ####

# file: microstegium_germination_model_2018_density_exp
# author: Amy Kendig
# date last edited: 10/16/20
# goal: analyze Microstegium germination


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
set1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
set2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")


#### edit data ####

# notes
unique(set1$notes_check_1) # seedlings with lesions in notes
unique(set1$notes_check_2) # may be a contaminate on plate
unique(set1$notes_check_3)
unique(set1$notes_germination_check_1) # some plates were put into fridge during one day
unique(set1$notes_germination_final) 
unique(set2$notes) # may be a contaminate on plate

# check columns
unique(set1$germination_check_1)
unique(set1$germination_final)

# condense info in set 1
mvGermD1Dat <- set1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final)) %>%
  select(site_plot, trial, seeds, germination_final) %>%
  full_join(set2 %>%
              select(site_plot, trial, seeds, germination_final)) %>%
  mutate(site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           recode("F" = "fungicide", "W" = "water"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot))


#### initial visualizations ####

ggplot(mvGermD1Dat, aes(x = as.factor(plot), y = germination_final/seeds, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(mvGermD1Dat, aes(x = treatment, y = germination_final/seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### model ####

# initial fit
mvGermD1Mod1 <- brm(germination_final | trials(seeds) ~ fungicide + (1|site/plotr),
                     data = mvGermD1Dat, family = binomial,
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 6000, warmup = 1000, chains = 1)
# 82 divergent transitions
summary(mvGermD1Mod1)

# increase chains and adapt delta
mvGermD1Mod2 <- update(mvGermD1Mod1, chains = 3,
                        control = list(adapt_delta = 0.999, max_treedepth = 15))
summary(mvGermD1Mod2)
plot(mvGermD1Mod2)
pp_check(mvGermD1Mod2, nsamples = 50)

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), seeds = c(30, 30)) %>%
  mutate(germination_final = fitted(mvGermD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Estimate"],
         germination_final_lower = fitted(mvGermD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q2.5"],
         germination_final_upper = fitted(mvGermD1Mod2, newdata = ., re_formula = NA, type = "response")[, "Q97.5"],
         treatment = ifelse(fungicide == 0, "water", "fungicide"))

# fit figure
ggplot(mvGermD1Dat, aes(treatment, germination_final/seeds)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  geom_point(data = fitDat, size = 2, color = "red", alpha = 0.5) +
  geom_errorbar(data = fitDat, aes(ymin = germination_final_lower/seeds, ymax = germination_final_upper/seeds), color = "red", width = 0, alpha = 0.5) +
  theme_bw()


#### output ####

save(mvGermD1Mod2, file = "output/microstegium_germination_model_2018_density_exp.rda")

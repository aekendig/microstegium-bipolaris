##### info ####

# file: mv-degredation-by-treatment-2018
# author: Amy Kendig
# date last edited: 7/1/19
# goal: see how fungicide treatment and disease severity affected Mv biomass degredation


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

# run leaf scan files
source("./code/mv-leaf-scans-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf")))

# run covariate file
source("./code/covariate-data-processing-2018.R")
rm(list = setdiff(ls(), c("mleaf", "covar", "esurv", "msurv", "bgd")))

# import data
mb <- read_csv("./data/mv-biomass-oct-2018-density-exp.csv")
ml <- read_csv("./data/mv-litter-biomass-apr-2019-density-exp.csv")
til <- read_csv("./data/mv-disease-sep-2018-density-exp.csv")

# plot parameters
colpal = c("#3CBB75FF", "#39568CFF")
sm_txt = 12
lg_txt = 14
an_txt = 3

# functions
plot_intervals <- function(data) {
  ggplot(data, aes(y = parameter, yend = parameter)) + 
    geom_vline(xintercept = 0) +
    geom_segment(aes(x = ll, xend = hh), size = 1) + 
    geom_segment(aes(x = l, xend = h), size = 2) +
    geom_point(aes(x = m), size = 3, color = "red") +
    theme_bw() +
    facet_wrap(~model)
}


#### edit data ####

# edit leaf scan data for focals
ml_f <- mleaf %>%
  ungroup() %>%
  filter(month == "September" & focal == 1) %>%
  select(site, plot, treatment, ID, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
  mutate(plot = as.numeric(plot)) %>%
  full_join(filter(til, ID %in% c("1", "2", "3"))) %>%
  mutate(leaves_infec2 = ifelse(leaves_infec == 0 & !is.na(leaf_area.pix), 1, leaves_infec),
         ind_severity = (lesion_area.pix - green_area.pix) * leaves_infec2 / (leaf_area.pix * leaves_tot)) %>%
  group_by(site, plot, treatment) %>%
  summarise(severity = mean(ind_severity, na.rm = T),
            n = sum(!is.na(ind_severity)),
            severity_se = sd(ind_severity, na.rm = T)/n) %>%
  filter(plot == 4)

# edit leaf scan data for background
ml_b <- mleaf %>%
  ungroup() %>%
  mutate(plot = as.numeric(plot)) %>%
  filter(month == "September" & focal == 0 & plot == 4) %>%
  select(site, plot, treatment, lesion_area.pix, green_area.pix, leaf_area.pix) %>%
  full_join(til %>% 
              filter(!(ID %in% c("1", "2", "3")) & plot == 4) %>%
              group_by(site, plot, treatment) %>%
              summarize(n = sum(!is.na(leaves_tot)),
                leaves_infec = sum(leaves_infec, na.rm = T),
                        leaves_tot = sum(leaves_tot, na.rm = T))) %>%
  mutate(severity = (lesion_area.pix - green_area.pix) * leaves_infec / (leaf_area.pix * leaves_tot))

# check notes
filter(mb, plot == 4) %>%
  select(processing_notes) %>%
  unique()

# edit values and merge with other data
dat <- ml %>%
  select(-c(date, mv_weigh_date)) %>%
  rename(litter.g = mv.g) %>%
  full_join(filter(mb, plot == 4)) %>%
  left_join(covar) %>%
  mutate(log_bio.g = log(bio.g),
         log_litter.g = log(litter.g),
         disease = ifelse(treatment == "fungicide", 0, 1),
         bio_diff.g = litter.g - bio.g)


#### visualize ####

dat %>%
  ggplot(aes(x = bio.g, y = litter.g, color = treatment)) +
  geom_point()

dat %>%
  ggplot(aes(x = treatment, y = bio_diff.g, shape = site)) +
  geom_point()

dat %>%
  ggplot(aes(x = treatment, y = bio_diff.g)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2)

summary(dat$bio_diff.g)
summary(dat$bio.g)
summary(dat$litter.g)


#### make data long ####

dat2 <- dat %>%
  select(site, plot, treatment, sm_adj, cc_adj, pm_adj, bio.g, litter.g) %>%
  gather(month, biomass, -c(site, plot, treatment, sm_adj, cc_adj, pm_adj)) %>%
  mutate(month = recode(month, bio.g = "October", litter.g = "April"),
         month = factor(month, levels = c("October", "April")),
         log_bio = log(biomass))


#### stats ####

# log-transformed, uninformative priors
mod <- brm(data = dat2, family = gaussian,
           log_bio ~ month * treatment + sm_adj + cc_adj + pm_adj,
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b)),
           iter = 6000, warmup = 1000, chains = 3, cores = 2)

summary(mod)
plot(mod)

# non-transformed, uninformative priors
mod2 <- brm(data = dat2, family = gaussian,
           biomass ~ month * treatment + sm_adj + cc_adj + pm_adj,
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b)),
           iter = 6000, warmup = 1000, chains = 3, cores = 2)

summary(mod2)
plot(mod2)

# compare transformed and non
pp_check(mod, nsamples = 50)
pp_check(mod2, nsamples = 50)

dat2 <- dat2 %>%
  mutate(pred_tu = fitted(mod)[,1],
         pred_nu = fitted(mod2)[,1])

cor.test(dat2$log_bio, dat2$pred_tu) # 0.71
cor.test(dat2$biomass, dat2$pred_nu) # 0.66

# transformed, informative priors
# Kourtev et al. 2002: change in biomass from September to March
mean(c(0.58 * log(6) - 0.59 * log(6), 0.45 * log(6) - 0.41 * log(6)))
# Stricker et al. 2016: 39% increase in biomass with fungicide
-log(1.39)

mod3 <- brm(data = dat2, family = gaussian,
            log_bio ~ month * treatment + sm_adj + cc_adj + pm_adj,
            prior <- c(prior(normal(0, 100), class = Intercept),
                       prior(normal(0.03, 10), class = b, coef = monthApril),
                       prior(normal(-0.33, 10), class = b, coef = treatmentwater),
                       prior(normal(0, 10), class = b, coef = monthApril:treatmentwater),
                       prior(normal(0, 10), class = b, coef = sm_adj),
                       prior(normal(0, 10), class = b, coef = cc_adj),
                       prior(normal(0, 10), class = b, coef = pm_adj)),
            iter = 6000, warmup = 1000, chains = 3, cores = 2)

summary(mod3)
plot(mod3)

imod <- as.array(mod) %>%
  mcmc_intervals_data(pars = paste("b_", rownames(fixef(mod)), sep = "")) %>%
  mutate(model = "uninformative")

imod3 <- as.array(mod3) %>%
  mcmc_intervals_data(pars = paste("b_", rownames(fixef(mod3)), sep = "")) %>%
  mutate(model = "informative")

plot_intervals(full_join(imod, imod3)) # basically no effect

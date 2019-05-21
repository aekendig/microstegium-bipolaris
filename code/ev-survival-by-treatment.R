##### info ####

# file: ev-survival-by-treatment-2018
# author: Amy Kendig
# date last edited: 4/30/19
# goal: evaluate treatment effects on Ev survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# run survival data files
source("./code/ev-survival-data-processing-2018.R")

# clear everything except survival data
rm(list = setdiff(ls(), "esurv"))

# run covariate data files
source("./code/covariate-data-processing-2018.R")

# clear everything except survival data
rm(list = setdiff(ls(), c("esurv", "covar")))

# import data
trt <- read_csv("./data/plot-treatments-2018-density-exp.csv")
# import covariates


#### edit data ####

# merge with treatment and covariates
d <- full_join(esurv, trt) %>%
  full_join(filter(covar, !(site == "D4" & treatment == "fungicide" & plot %in% c(8, 10))))

# select September data (summer survival)
ds <- d %>%
  filter(month == "September" & !is.na(survival) & !is.na(treatment))

# select April data (winter survival given summer survival)
da <- d %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April)


#### visualize data ####

# counts
ds %>%
  group_by(age, treatment, background, background_density) %>%
  summarise(n = length(survival), s = sum(survival)) %>%
  gather(key = "type", value = "count", -c(age, treatment, background, background_density)) %>%
  ggplot(aes(x = background_density, y = count, colour = treatment, shape = type)) +
  facet_grid(age~background, scales = "free") +
  geom_point()

# proportion
ds %>%
  group_by(age, treatment, background, background_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = background_density, y = p, colour = treatment)) +
  facet_grid(age~background, scales = "free") +
  geom_point(position = position_dodge(1))

# proportion April
da %>%
  group_by(age, treatment, background, background_density) %>%
  summarise(p = sum(survival) / length(survival)) %>%
  ggplot(aes(x = background_density, y = p, colour = treatment)) +
  facet_grid(age~background, scales = "free") +
  geom_point(position = position_dodge(1))


#### statistical models ####

# seedling by type, summer
dsst <- ds %>%
  filter(age == "seedling", !(background %in% "none")) %>%
  mutate(sm = scale(soil_moisture_oct),
         cc = scale(canopy_cover),
         mb = scale(mv_bio.g))

msst0 <- brm(data = dsst, family = bernoulli,
             survival ~ 1,
             prior = prior(normal(0, 10), class = Intercept),
             iter = 6000, warmup = 1000, chains = 1, cores = 1)
summary(msst0)
inv_logit_scaled(fixef(msst0)[1])
mean(dsst$survival)

msst1 <- brm(data = dsst, family = bernoulli,
             survival ~ 1 + (1|site),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, chains = 1, cores = 1)


## older code

# adult with adult, summer
maas <- brm(data = filter(ds, background %in% c("none", "Ev adult") & age == "adult"), 
            family = binomial,
            survival|trials(1) ~ background_density * treatment + soil_moisture_oct + canopy_cover + mv_bio.g + (1|site),
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
            iter = 6000, warmup = 1000, chains = 1, cores = 1,
            control = list(adapt_delta = 0.95))
# many divergent transitions - probably because most of the values are 1
summary(maas)
filter(ds, background %in% c("none", "Ev adult") & age == "adult") %>% summarise(s = sum(survival), n = length(survival)) # only 2 deaths

# adult with seedling, summer - will be same issue as above

# adult with Mv, summer
mams <- brm(data = filter(ds, background %in% c("none", "Mv seedling") & age == "adult"), 
            family = binomial,
            survival|trials(1) ~ background_density * treatment + soil_moisture_oct + canopy_cover + mv_bio.g + (1|site),
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
            iter = 6000, warmup = 1000, chains = 1, cores = 1,
            control = list(adapt_delta = 0.95))
summary(mams)


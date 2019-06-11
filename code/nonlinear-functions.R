##### info ####

# file: nonlinear-functions
# author: Amy Kendig
# date last edited: 6/5/19
# goal: demonstrations of nonlinear functions that can be used to describe species interactions


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)


#### define functions ####

# References
# Law and Watkinson 1987
# Levine and HilleRisLambers 2009
# Hart et al. 2018

# exponential
ex_fun <- function(y0, x, a) {
  y = y0 * exp(a * x)
  return(y)
}

# Beverton-Holt
bh_fun <- function(y0, x, a) {
  y = y0 / (1 + a * x)
  return(y)
}

# Beverton-Holt with exponent
bhe_fun <- function(y0, x, a, b) {
  y = y0 / (1 + a * x)^b
  return(y)
}

# non-linear isoclines
nli_fun <- function(y0, x, a) {
  y = y0 / (1 + x^a)
  return(y)
}

# linear
ln_fun <- function(y0, x, a) {
  y = y0 - a * x
  return(y)
}

# exponential/log
el_fun <- function(y0, x, a) {
  y = y0 * exp(a * log(x + 1))
  return(y)
}

# another Beverton-Holt with exponent
abhe_fun <- function(y0, x, a, b) {
  y = y0 / (1 + (a * x)^b)
  return(y)
}

# another linear
aln_fun <- function(y0, x, a) {
  y = 1 + y0 * (1 - a * x)
  return(y)
}


#### simulate data ####

# y0 is based on biomass
y0 = 14

# functions
funs = tibble(funs = c("abhe", "aln", "bh", "bhe", "el", "ex", "ln", "nli"))

# a is the rate of increase or decrease
sim_dat = expand.grid(0:64, c(-2:2, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1)) %>%
  merge(funs, all = T)

# name columns
colnames(sim_dat) = c("x", "a", "func")

# simulate y values
sim_dat <- sim_dat %>%
  mutate(y = case_when(func == "abhe" ~ abhe_fun(y0, x, a, 0.2),
                       func == "aln" ~ aln_fun(y0, x, a),
                       func == "bh" ~ bh_fun(y0, x, a),
                       func == "bhe" ~ bhe_fun(y0, x, a, 0.2),
                       func == "el" ~ el_fun(y0, x, a),
                       func == "ex" ~ ex_fun(y0, x, a),
                       func == "ln" ~ ln_fun(y0, x, a),
                       func == "nli" ~ nli_fun(y0, x, a)))

#### Figure ####

# all results
sim_dat %>%
  ggplot(aes(x = x, y = y, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~func, scales = "free")

# remove extreme values
sim_dat %>%
  filter(y < 50 & y > 0) %>%
  ggplot(aes(x = x, y = y, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~func, scales = "free")

# look at the nli and el models (saturating postiive effects)
sim_dat %>%
  filter(func %in% c("el", "nli")) %>%
  ggplot(aes(x = x, y = y, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~func, scales = "free")

# el does poorly with high a values
sim_dat %>%
  filter((func == "el" & a < 1) | func == "nli") %>%
  ggplot(aes(x = x, y = y, colour = as.factor(a))) + 
  geom_line() +
  facet_wrap(~func, scales = "free")

# nli has some weird behavior near zero
sim_dat %>%
  filter(func == "nli" & x < 5) %>%
  ggplot(aes(x = x, y = y, colour = as.factor(a))) + 
  geom_line()
# 0 raised to negative numbers or 0 isn't calculated as 0


#### conclusions ####

# the exponential log function works well as long as a is less than 1
# the non-linear isocline function works well as long as x is greater than 0
# the nli was analyzed in Hart et al. 2018 (competitive ability value available)
# but I don't think it makes sense that one neighbor automatically decreases the biomass to half of when the plant is alone. it seems to assume that the focal counts as 1 and when it is alone, it is at half of its carrying capacity (given facilitation)
# I'll use the exponential log function
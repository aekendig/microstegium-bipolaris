##### info ####

# file: annual_perennial_simulation_parameter_distributions_2019_density_exp
# author: Amy Kendig
# date last edited: 10/17/20
# goal: demonstrate sampling from statistical model distributions for simulation model


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(gtools)

# number of iterations
n_samps <- 1000 

# call model script
source("code/annual_perennial_simulation_model_2019_density_exp.R")


#### conditions ####

# simulation time
simtime <- 2

# initial conditions
A0 <- 100 # initial annual population size
S0 <- 10 # initial perennial seedling population size
P0 <- 100 # initial perennial adult population size
L0 <- 0 # initial annual litter amount


#### parameter distributions ####

# samples from parameter distributions
samps <- n_samps

# initiate list
params <- list()

# save parameter values from each sample
for(iter in 1:samps){
  
  params[[iter]] <- sim_fun(A0, S0, P0, L0, simtime, disease = 0, iter)[[2]]
  
}

# convert to dataframe
param_dat <- do.call(rbind, params)

# visualize
pdf("output/simulation_parameter_distributions.pdf")
ggplot(param_dat, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ description, 
             scales = "free",
             ncol = 2) +
  theme_bw()
dev.off()


#### notes ####

# for proportion variables, the distributions of the logit-transformed variable (sampled from the posterior distribution) look like normal distributions, doing the inverse transformation to obtain a proportion leads to bimodal distributions around zero and one

# log-transformed annual biomass and seeds produced are sampled from the functions, which is in a normal distribution, but when they are exponentiated (to biomass in g or number of seeds), the distribution is skewed

# perennial plants are more likely to produce seeds than not only when their biomass exceeds ~ 1 g

# the back-transformation for proportion variables doesn't work for values of 1 and taking the log of other variables doesn't work for values of 0


#### alternate figure ####

# transform variables
param_dat2 <- param_dat %>%
  mutate(value_adj = case_when(parameter %in% c("H.A", "H.S", "H.P", "W.S", "W.P", "G.A", "G.S") & value == 1 ~ 0.9999,
                               parameter %in% c("V.A", "V.S", "V.P", "Y.A", "Y.S", "Y.P") & value == 0 ~ 0.0001,
                               TRUE ~ value),
         transformed_value = case_when(parameter %in% c("H.A", "H.S", "H.P", "W.S", "W.P", "G.A", "G.S") ~ logit(value_adj),
                                       parameter %in% c("V.A", "V.S", "V.P", "Y.A", "Y.S", "Y.P") ~ log(value_adj),
                                       TRUE ~ value_adj))


# visualize
pdf("output/simulation_transformed_parameter_distributions.pdf")
ggplot(param_dat2, aes(x = transformed_value)) +
  geom_histogram() +
  facet_wrap(~ description, 
             scales = "free",
             ncol = 2) +
  theme_bw()
dev.off()

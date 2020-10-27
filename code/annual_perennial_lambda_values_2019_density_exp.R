##### info ####

# file: annual_perennial_lambda_values_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: lifetime reproduction


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidybayes)

# number of iterations
n_samps <- 1000

# call model script
source("code/annual_perennial_kortessis_model_2019_density_exp.R")


#### conditions ####

# simulation time
simtime <- 2

# initial conditions
A0 <- 0 # initial annual population size
S0 <- 0 # initial perennial seedling population size
P0 <- 0 # initial perennial adult population size
L0 <- 0 # initial annual litter amount

# samples from parameter distributions
samps <- n_samps


#### no disease simulations ####

# initiate list
params <- list()

# save parameter values from each sample
for(iter in 1:samps){
  
  mod <- sim_fun(A0, S0, P0, L0, simtime, disease = 0, iter)
  params[[iter]] <- mod[[2]]
  
}

# convert to dataframe
param_dat <- do.call(rbind, params)


#### disease simulations ####

# initiate list
params2 <- list()

# save parameter values from each sample
for(iter in 1:samps){
  
  mod2 <- sim_fun(A0, S0, P0, L0, simtime, disease = 1, iter)
  params2[[iter]] <- mod2[[2]]
  
}

# convert to dataframe
param_dat2 <- do.call(rbind, params2)


#### figure ####

# relative fitness
rel_param_dat <- param_dat %>%
  mutate(Disease = "without disease") %>%
  full_join(param_dat2 %>%
              mutate(Disease = "with disease")) %>%
  select(-description) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(rel_fit = l.A / l.P,
         Disease = fct_relevel(Disease, "without disease"))

# figure
pdf("output/annual_perennial_lambda_values_2019_density_exp.pdf", width = 3, height = 3)
ggplot(rel_param_dat, aes(Disease, rel_fit, color = Disease)) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  stat_pointinterval(point_interval = mean_hdci, .width = 0.95, interval_size = 1, point_size = 3) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  ylab(expression(paste("Relative fitness (", italic("M. vimineum"), "/", italic("E. virginicus"), ")", sep = ""))) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 10),
        strip.background = element_blank())
dev.off()

# test different
rel_param_dat %>%
  select(time:Disease, rel_fit) %>%
  mutate(Disease = recode(Disease, "without disease" = "no", "with disease" = "yes")) %>%
  pivot_wider(names_from = Disease, values_from = rel_fit) %>%
  mutate(dis_diff = no - yes) %>%
  mean_hdi(dis_diff)

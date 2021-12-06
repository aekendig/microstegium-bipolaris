#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(deSolve)
library(cowplot)

# import parameters
parms_all <- read_csv("output/model_parameters_2018_2019_density_exp.csv")
# model_parameters_2018_2019_density_exp

# continuous models
source("code/continuous_AFP_model.R")
source("code/continuous_AF_model.R")
source("code/continuous_AP_model.R")
source("code/continuous_FP_model.R")
source("code/continuous_A_model.R")
source("code/continuous_F_model.R")
source("code/continuous_P_model.R")

# discrete model
source("code/discrete_AFP_no_seed_infection_model.R")


#### parameters ####

# growing season days
gs_days <- 160
gs_days_seq <- seq(0, gs_days, by = 1)

# simulation years
years <- 100

# control parameters
parms_ctrl <- parms_all %>%
  filter(Treatment == "control" | is.na(Treatment))

# fungicide parameters
parms_fung <- parms_all %>%
  filter(Treatment == "fungicide" | is.na(Treatment))

# conidia per 1 g litter
conidia_init <- as.numeric(parms_all[parms_all$Parameter == "h", "Estimate"])

# initial plant sizes
bio_A_init <- as.numeric(parms_all[parms_all$Parameter == "b_A", "Estimate"])


#### annual alone ####

# continuous input values
xstart <- c(LogB_A = log(bio_A_init * 10), I_A = 0, D = 0, C = conidia_init)

# apply continuous function
modA_cont_ctrl <- ode(y = xstart, times = gs_days_seq, 
                      func = cont_A_mod, parms = parms_ctrl) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(B_A = exp(LogB_A),
         S_A = B_A - I_A) %>%
  select(-LogB_A) %>%
  rename("I_C" = "C", "I_D" = "D") %>%
  pivot_longer(cols = -time,
               names_to = c("state", "plant"),
               names_pattern = "(.)_(.)",
               values_to = "pop")

# figure of continuous results
modA_cont_ctrl %>%
  ggplot(aes(x = time, y = pop, color = state)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ plant, scales = "free_y")

# discrete input values
modA_disc_ctrl <- disc_AFP_mod(A0 = 10, F0 = 0, P0 = 0, L0 = 0, C0 = conidia_init, IA0 = 0, IF0 = 0, IP0 = 0, BP0 = 0, BF0 = 0, simtime = years, gs_time = gs_days, disc_parms = parms_ctrl, con_parms = parms_ctrl)

# create figure
modA_disc_ctrl %>%
  filter(species %in% c("Annual", "Litter")) %>%
  ggplot(aes(x = time, y = N)) +
  geom_line() +
  theme_classic() +
  facet_wrap(~ species, scales = "free_y")

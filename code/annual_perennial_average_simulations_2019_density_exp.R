##### info ####

# file: annual_perennial_simulation_model_2019_density_exp
# author: Amy Kendig
# date last edited: 3/29/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)

# load models
load("output/evA_growing_season_survival_model_2019_density_exp.rda")
load("output/evS_winter_survival_model_2018_density_exp.rda")
load("output/evA_winter_survival_model_2018_density_exp.rda")
load("output/mv_germination_model_2018_density_exp.rda")
load("output/ev_germination_model_2018_2019_density_exp.rda")
load("output/mv_growing_season_survival_model_2019_density_exp.rda")
load("output/evS_growing_season_survival_model_2019_density_exp.rda")
load("output/microstegium_litter_establishment_bh_model_2018_greenhouse_exp.rda")
load("output/elymus_litter_establishment_bh_model_2018_greenhouse_exp.rda")
load("output/mv_biomass_mv_density_model_2019_density_exp.rda")
load("output/evS_biomass_mv_density_model_2019_density_exp.rda")
load("output/evA_biomass_mv_density_model_2019_density_exp.rda")
load("output/mv_biomass_evS_density_model_2019_density_exp.rda")
load("output/evS_biomass_evS_density_model_2019_density_exp.rda")
load("output/evA_biomass_evS_density_model_2019_density_exp.rda")
load("output/mv_biomass_evA_density_model_2019_density_exp.rda")
load("output/evS_biomass_evA_density_model_2019_density_exp.rda")
load("output/evA_biomass_evA_density_model_2019_density_exp.rda")
load("output/mv_seeds_per_biomass_model_2019_density_exp.rda")
load("output/evS_seeds_per_biomass_model_2019_density_exp.rda")
load("output/evA_seeds_per_biomass_model_2019_density_exp.rda")
load("output/mv_plot_biomass_density_model_2019_density_exp.rda")
load("output/focal_growth_no_background_model_2019_density_exp.rda")
load("output/evS_tillers_mv_density_model_2018_density_exp.rda")


#### parameters from models ####

# survival
u.P <- posterior_samples(evASurvD2Mod) %>%
  transmute(u.P = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(u.P) %>%
  mean()

w.S <- posterior_samples(evSWinSurvD1Mod) %>%
  transmute(w.S = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(w.S) %>%
  mean()
w.P <- posterior_samples(evAWinSurvD1Mod) %>%
  transmute(w.P = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(w.P) %>%
  mean()

# germination
g.A <- posterior_samples(mvGermD1Mod) %>%
  transmute(g.A = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(g.A) %>%
  mean()
g.S <- posterior_samples(evGermMod) %>%
  transmute(g.S = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(g.S) %>%
  mean()

# establishment
e.A <- posterior_samples(mvSurvD2Mod) %>%
  transmute(e.A = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(e.A) %>%
  mean()
e.S <- posterior_samples(evSSurvD2Mod) %>%
  transmute(e.S = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide))) %>%
  pull(e.S) %>%
  mean()

# litter effects
beta.A <- posterior_samples(mvLitEstBhMod2) %>%
  transmute(beta.A = b_betaL_Intercept) %>%
  pull(beta.A) %>%
  mean()
beta.S <- posterior_samples(evLitEstBhMod2) %>%
  transmute(beta.S = b_betaL_Intercept) %>%
  pull(beta.S) %>%
  mean()

# biomass
b.A <- posterior_samples(mvMvD2Mod) %>%
  transmute(b.A = b_b0_treatmentfungicide) %>%
  pull(b.A) %>%
  mean()
b.S <- posterior_samples(evSMvD2Mod) %>%
  transmute(b.S = b_b0_treatmentfungicide) %>%
  pull(b.S) %>%
  mean()
b.P <- posterior_samples(evAMvD2Mod) %>%
  transmute(b.P = b_b0_treatmentfungicide) %>%
  pull(b.P) %>%
  mean()

# competition coefficients
alpha.AA <- posterior_samples(mvMvD2Mod) %>%
  transmute(alpha.AA = b_alpha_treatmentfungicide) %>%
  pull(alpha.AA) %>%
  mean()
alpha.SA <- posterior_samples(evSMvD2Mod) %>%
  transmute(alpha.SA = b_alpha_treatmentfungicide) %>%
  pull(alpha.SA) %>%
  mean()
alpha.PA <- posterior_samples(evAMvD2Mod) %>%
  transmute(alpha.PA = b_alpha_treatmentfungicide) %>%
  pull(alpha.PA) %>%
  mean()

alpha.AS <- posterior_samples(mvEvSD2Mod) %>%
  transmute(alpha.AS = b_alpha_treatmentfungicide) %>%
  pull(alpha.AS) %>%
  mean()
alpha.SS <- posterior_samples(evSEvSD2Mod) %>%
  transmute(alpha.SS = b_alpha_treatmentfungicide) %>%
  pull(alpha.SS) %>%
  mean()
alpha.PS <- posterior_samples(evAEvSD2Mod) %>%
  transmute(alpha.PS = b_alpha_treatmentfungicide) %>%
  pull(alpha.PS) %>%
  mean()

alpha.AP <- posterior_samples(mvEvAD2Mod) %>%
  transmute(alpha.AP = b_alpha_treatmentfungicide) %>%
  pull(alpha.AP) %>%
  mean()
alpha.SP <- posterior_samples(evSEvAD2Mod) %>%
  transmute(alpha.SP = b_alpha_treatmentfungicide) %>%
  pull(alpha.SP) %>%
  mean()
alpha.PP <- posterior_samples(evAEvAD2Mod) %>%
  transmute(alpha.PP = b_alpha_treatmentfungicide) %>%
  pull(alpha.PP) %>%
  mean()

# seeds
y.A <- posterior_samples(mvSeedsBioD2Mod) %>%
  transmute(y.A = b_Intercept + b_fungicide) %>%
  pull(y.A) %>%
  mean()
y.S <- posterior_samples(evSSeedsBioD2Mod) %>%
  transmute(y.S = b_Intercept + b_fungicide) %>%
  pull(y.S) %>%
  mean()
y.P <- posterior_samples(evASeedsBioD2Mod) %>%
  transmute(y.P = b_Intercept + b_fungicide) %>%
  pull(y.P) %>%
  mean()

yb.A <- posterior_samples(mvSeedsBioD2Mod) %>%
  rename("b_fungi_log_bio" = "b_fungicide:log_bio") %>%
  transmute(yb.A = b_log_bio + b_fungi_log_bio) %>%
  pull(yb.A) %>%
  mean()
yb.S <- posterior_samples(evSSeedsBioD2Mod) %>%
  rename("b_fungi_log_bio" = "b_fungicide:log_bio") %>%
  transmute(yb.S = b_log_bio + b_fungi_log_bio) %>%
  pull(yb.S) %>%
  mean()
yb.P <- posterior_samples(evASeedsBioD2Mod) %>%
  rename("b_fungi_log_bio" = "b_fungicide:log_bio") %>%
  transmute(yb.P = b_log_bio + b_fungi_log_bio) %>%
  pull(yb.P) %>%
  mean()

# perennial lifespan
l.P <- ifelse(w.P * u.P > 0.99, 0.99, w.P * u.P)


#### disease effects ####

# increases perennial germination
evGermDE <- posterior_samples(evGermMod) %>%
  transmute(fung = exp(b_Intercept + b_fungicide)/(1 + exp(b_Intercept + b_fungicide)),
            water = exp(b_Intercept)/(1 + exp(b_Intercept)),
            DE = water / fung) %>%
  median_hdi(DE) %>% pull(DE)

# decreases annual biomass
mvBioDE <- posterior_samples(mvBioDensMod) %>%
  transmute(fung = (67 * b_b0_treatmentfungicide) / (1 + b_alpha_treatmentfungicide * 67),
            water = (67 * b_b0_treatmentwater) / (1 + b_alpha_treatmentwater * 67),
            DE = water / fung) %>%
  median_hdi(DE) %>% pull(DE)

# decreases perennial biomass
evBioDE <- posterior_samples(growthD2Mod) %>%
  rename(b_fungicide_Ev_seedling = "b_fungicide:plant_groupEv_seedling") %>%
  transmute(water = exp(b_Intercept + b_plant_groupEv_seedling),
            fung = exp(b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_Ev_seedling),
            DE = water / fung) %>%
  median_hdi(DE) %>% pull(DE)

# increases the competitive effect of annual on perennial
alphaSADE <- posterior_samples(evSMvD1Mod) %>%
  transmute(fung = b_alpha_treatmentfungicide,
            water = b_alpha_treatmentwater,
            DE = water / fung) %>%
  median_hdi(DE) %>% pull(DE)


#### constant parameters ####

# annual seed survival
s.A <- s.S <- 0.1 
# median of Redwood et al. 2018, Mv = 0.15; Garrison and Stier 2010, Ev = 0.05

# litter decomposition 
d <- 0.35
# median of DeMeester and Richter 2010, Mv = 0.6; Hobbie et al. 2008, min = 0.1


#### model ####

sim_fun = function(A0, S0, P0, L0, simtime, g.S.DE, b.A.DE, b.S.DE, alpha.SA.DE){
  
  # initialize populations
  A <- rep(NA,simtime)
  S <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  
  A[1] <- A0
  S[1] <- S0
  P[1] <- P0
  L[1] <- L0
  
  
  # germination
  g.S <- g.S * g.S.DE
  g.S <- ifelse(g.S > 1, 1, g.S) # restrict germination fraction to maximum 1
  
  # maximum growth
  b.A <- b.A * b.A.DE
  b.S <- b.S * b.S.DE
  
  # competitive effect of annual
  alpha.SA <- alpha.SA * alpha.SA.DE
  
  # initialize parameter vectors
  E.A_vec <- rep(NA, simtime-1)
  E.S_vec <- rep(NA, simtime-1)
  B.A_vec <- rep(NA, simtime-1)
  B.S_vec <- rep(NA, simtime-1)
  B.P_vec <- rep(NA, simtime-1)
  Y.A_vec <- rep(NA, simtime-1)
  Y.S_vec <- rep(NA, simtime-1)
  Y.P_vec <- rep(NA, simtime-1)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # seedling establishment
    E.A <- e.A / (1 + beta.A * L[t])
    E.S <- e.S / (1 + beta.S * L[t])
    
    # reduce growth due to competition
    B.A <- b.A / (1 + alpha.AA * g.A * E.A * A[t] + alpha.AS * g.S * E.S * S[t] + alpha.AP * P[t])
    B.S <- b.S / (1 + alpha.SA * g.A * E.A * A[t] + alpha.SS * g.S * E.S * S[t] + alpha.SP * P[t])
    B.P <- b.P / (1 + alpha.PA * g.A * E.A * A[t] + alpha.PS * g.S * E.S * S[t] + alpha.PP * P[t])
    
    # seed production (log-log regression)
    Y.A <- exp(y.A + yb.A * log(B.A)) - 1
    Y.S <- exp(y.S + yb.S * log(B.S)) - 1
    Y.P <- exp(y.P + yb.P * log(B.P)) - 1
    
    # population size
    A[t+1] = s.A * (1-g.A) * A[t] + g.A * E.A * Y.A * A[t]  
    L[t+1] = g.A * E.A * B.A * A[t] + g.S * E.S * B.S * S[t] + B.P * P[t] + (1 - d) * L[t]    
    S[t+1] = s.S * (1-g.S) * S[t] + g.S * E.S * Y.S * S[t] + u.P * Y.P * P[t]
    P[t+1] = l.P * P[t] + g.S * E.S * w.S * S[t]  
    
    # correct to prevent negative numbers
    A[t+1] = ifelse(A[t+1] < 1, 0, A[t+1])
    L[t+1] = ifelse(L[t+1] < 1, 0, L[t+1])
    S[t+1] = ifelse(S[t+1] < 1, 0, S[t+1])
    P[t+1] = ifelse(P[t+1] < 1, 0, P[t+1])
    
    # save parameter value
    E.A_vec[t] <- E.A
    E.S_vec[t] <- E.S
    B.A_vec[t] <- B.A
    B.S_vec[t] <- B.S
    B.P_vec[t] <- B.P
    Y.A_vec[t] <- Y.A
    Y.S_vec[t] <- Y.S
    Y.P_vec[t] <- Y.P
  }
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 4),
                species = rep(c("Perennial seedling", 
                                "Perennial adult", 
                                "Annual", 
                                "Litter"), 
                              each = simtime),
                N = c(S, P, A, L))
  
  # parameter data
  dfP <- tibble(time = rep(1:(simtime-1), 8),
                parameter = rep(c("E.A", "E.S", "B.A", "B.S", "B.P", "Y.A", "Y.S", "Y.P"), 
                                each = simtime - 1),
                description = rep(c("annual seedling establishment",
                                    "perennial seedling establishment",
                                    "annual biomass production",
                                    "perennial seedling biomass production",
                                    "perennial adult biomass production",
                                    "annual seed production",
                                    "perennial seedling seed production",
                                    "perennial adult seed production"), 
                                  each = simtime - 1),
                value = c(E.A_vec, E.S_vec, B.A_vec, B.S_vec, B.P_vec, Y.A_vec, Y.S_vec, Y.P_vec))
  
  # return
  return(list(dfN, dfP))
}


#### simulations ####

# simulation time
simtimeR <- 600
simtimeI <- 200

# initial conditions
A0E <- 0 # initial annual population size
S0E <- 0 # initial perennial seedling population size
P0E <- 10 # initial perennial adult population size
L0E <- 0 # initial annual litter amount

# function to simulate invasions
inv_fun <- function(g_S_DE, b_A_DE, b_S_DE, alpha_SA_DE){
  
  # establish residents
  mod_est <- sim_fun(A0 = A0E, S0 = S0E, P0 = P0E, L0 = L0E, simtime = simtimeR, 
                     g.S.DE = g_S_DE, b.A.DE = b_A_DE, b.S.DE = b_S_DE, alpha.SA.DE = alpha_SA_DE)
  
  # extract final abundances
  S0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Perennial seedling") %>% pull(N)
  P0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Perennial adult") %>% pull(N)
  L0Ei <- mod_est[[1]] %>% filter(time == simtimeR & species == "Litter") %>% pull(N)
  
  # invasion
  mod_inv <- sim_fun(A0 = 10, S0 = S0Ei, P0 = P0Ei, L0 = L0Ei, simtime = simtimeI, 
                     g.S.DE = g_S_DE, b.A.DE = b_A_DE, b.S.DE = b_S_DE, alpha.SA.DE = alpha_SA_DE)
  
  # output
  mod_out <- mod_est[[1]] %>% 
    full_join(mod_inv[[1]] %>%
                mutate(time = time + simtimeR))
  return(mod_out)
}

# no disease
mod_0 <- inv_fun(g_S_DE = 1, b_A_DE = 1, b_S_DE = 1, alpha_SA_DE = 1) %>%
  mutate(DE = 0)

# observed disease
mod_1 <- inv_fun(g_S_DE = 1, b_A_DE = mvBioDE, b_S_DE = evBioDE, alpha_SA_DE = alphaSADE) %>%
  mutate(DE = 1)

# observed disease
mod_10 <- inv_fun(g_S_DE = 1, b_A_DE = mvBioDE/10, b_S_DE = evBioDE/10, alpha_SA_DE = alphaSADE*10) %>%
  mutate(DE = 10)


#### combine data ####

abund <- mod_0 %>%
  full_join(mod_1) %>%
  full_join(mod_10)



#### figures ####

ggplot(abund, aes(time, N, linetype = as.factor(DE))) +
  geom_line() +
  facet_wrap(~ species, scales = "free")

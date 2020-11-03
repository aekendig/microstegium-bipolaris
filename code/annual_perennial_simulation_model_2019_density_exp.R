##### info ####

# file: annual_perennial_simulation_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# number of samples
# n_samps <- 10 # specify in scripts that source this one

# load packages
library(tidyverse)
library(brms)

# load scripts
source("code/microstegium_establishment_parameter_2019_density_exp.R")
source("code/elymus_seedling_establishment_parameter_2019_density_exp.R")

source("code/elymus_adult_gs_survival_parameter_2019_density_exp.R")

source("code/elymus_seedling_ngs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_ngs_survival_parameter_2019_density_exp.R")

source("code/microstegium_biomass_parameter_2019_density_exp.R")
source("code/elymus_seedling_biomass_parameter_2019_density_exp.R")
source("code/elymus_adult_biomass_parameter_2019_density_exp.R")

source("code/microstegium_seed_production_parameter_2019_density_exp.R")
source("code/elymus_seedling_seed_production_parameter_2019_density_exp.R")
source("code/elymus_adult_seed_production_parameter_2019_density_exp.R")

source("code/microstegium_germination_parameter_2018_density_exp.R")
source("code/elymus_germination_parameter_2019_density_exp.R")

# constant parameters
s.A <- 0.15 # annual seed survival (Redwood et al. 2018)
s.S <- 0.05 # perennial seed survival (Garrison and Stier 2010)
d <- 0.61 # annual litter decomposition (DeMeester and Richter 2010)


#### model ####

sim_fun = function(A0, S0, P0, L0, simtime, disease, iter){
  
  # initialize populations
  A <- rep(NA,simtime)
  S <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  
  A[1] <- A0
  S[1] <- S0
  P[1] <- P0
  L[1] <- L0
  
  # initialize parameter vectors
  E.A_vec <- rep(NA, simtime-1)
  E.S_vec <- rep(NA, simtime-1)
  u.P_vec <- rep(NA, simtime-1)
  w.S_vec <- rep(NA, simtime-1)
  w.P_vec <- rep(NA, simtime-1)
  B.A_vec <- rep(NA, simtime-1)
  B.S_vec <- rep(NA, simtime-1)
  B.P_vec <- rep(NA, simtime-1)
  Y.A_vec <- rep(NA, simtime-1)
  Y.S_vec <- rep(NA, simtime-1)
  Y.P_vec <- rep(NA, simtime-1)
  g.A_vec <- rep(NA, simtime-1)
  g.S_vec <- rep(NA, simtime-1)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # perennial adult growing season survival
    u.P <- u_P_fun(disease, iter)
    
    # perennial non-growing season survival
    w.S <- w_S_fun(disease, iter)
    w.P <- w_P_fun(disease, iter)
    
    # germination depends on disease and decreases due to litter
    g.A <- g_A_fun(disease, iter)
    g.S <- g_S_fun(disease, iter)
    
    # seedling establishment
    E.A <- E_A_fun(disease, g.A, A[t], g.S, S[t], P[t], L[t], iter)
    E.S <- E_S_fun(disease, g.A, A[t], g.S, S[t], P[t], L[t], iter)
    
    # reduce growth due to competition
    B.A <- B_A_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter)
    B.S <- B_S_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter)
    B.P <- B_P_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter) 
    
    # reduce seed production due to competition
    Y.A <- Y_A_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter)
    Y.S <- Y_S_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter)
    Y.P <- Y_P_fun(disease, g.A, E.A, A[t], g.S, E.S, S[t], P[t], iter)
    
    # perennial lifespan
    l.P <- ifelse(w.P * u.P > 0.99, 0.99, w.P * u.P)
    
    # population size
    A[t+1] = s.A * (1-g.A) * A[t] + g.A * E.A * Y.A * A[t]  
    L[t+1] = g.A * E.A * B.A * A[t] + g.S * E.S * B.S * S[t] + u.P * B.P * P[t] + (1 - d) * L[t]    
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
    u.P_vec[t] <- u.P
    w.S_vec[t] <- w.S
    w.P_vec[t] <- w.P
    B.A_vec[t] <- B.A
    B.S_vec[t] <- B.S
    B.P_vec[t] <- B.P
    Y.A_vec[t] <- Y.A
    Y.S_vec[t] <- Y.S
    Y.P_vec[t] <- Y.P
    g.A_vec[t] <- g.A
    g.S_vec[t] <- g.S
  }
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 4),
                species = rep(c("Perennial seedling", 
                                "Perennial adult", 
                                "Annual", 
                                "Litter"), 
                              each = simtime),
                N = c(S, P, A, L),
                iteration = iter)

  # parameter data
  dfP <- tibble(time = rep(1:(simtime-1), 13),
                parameter = rep(c("E.A", "E.S", "u.P", "w.S", "w.P", "B.A", "B.S", "B.P", "Y.A", "Y.S", "Y.P", "g.A", "g.S"), 
                                each = simtime - 1),
                description = rep(c("annual seedling establishment",
                                    "perennial seedling establishment",
                                    "perennial adult growing season survival",
                                    "perennial seedling non-growing season survival",
                                    "perennial adult non-growing season survival",
                                    "annual biomass production",
                                    "perennial seedling biomass production",
                                    "perennial adult biomass production",
                                    "annual seed production",
                                    "perennial seedling seed production",
                                    "perennial adult seed production",
                                    "annual germination",
                                    "perennial germination"), 
                                  each = simtime - 1),
                value = c(E.A_vec, E.S_vec, u.P_vec, w.S_vec, w.P_vec, B.A_vec, B.S_vec, B.P_vec, Y.A_vec, Y.S_vec, Y.P_vec, g.A_vec, g.S_vec),
                iteration = iter)

  # return
  return(list(dfN, dfP))
}




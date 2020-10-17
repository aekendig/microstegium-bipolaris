##### info ####

# file: annual_perennial_simulation_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/17/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# number of samples
# n_samps <- 10 # specify in scripts that source this one

# load packages
library(tidyverse)
library(brms)

# load scripts
source("code/microstegium_gs_survival_parameter_2019_density_exp.R")
source("code/elymus_seedling_gs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_gs_survival_parameter_2019_density_exp.R")

source("code/elymus_seedling_ngs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_ngs_survival_parameter_2019_density_exp.R")

source("code/microstegium_biomass_parameter_2019_density_exp.R")
source("code/elymus_seedling_biomass_parameter_2019_density_exp.R")
source("code/elymus_adult_biomass_parameter_2019_density_exp.R")

source("code/microstegium_seeds_per_biomass_parameter_2019_density_exp.R")
source("code/elymus_seeds_per_biomass_parameter_2019_density_exp.R")

source("code/microstegium_germination_parameter_2018_density_exp.R")
source("code/elymus_germination_parameter_2019_density_exp.R")

# constant parameters
s.A <- 0.05 # annual seed survival (Redwood et al. 2018)
s.S <- 0.09 # perennial seed survival (Garrison and Stier 2010)
b <- 0.33 # annual litter decomposition (DeMeester and Richter 2010)


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
  H.A_vec <- rep(NA, simtime-1)
  H.S_vec <- rep(NA, simtime-1)
  H.P_vec <- rep(NA, simtime-1)
  W.S_vec <- rep(NA, simtime-1)
  W.P_vec <- rep(NA, simtime-1)
  V.A_vec <- rep(NA, simtime-1)
  V.S_vec <- rep(NA, simtime-1)
  V.P_vec <- rep(NA, simtime-1)
  Y.A_vec <- rep(NA, simtime-1)
  Y.S_vec <- rep(NA, simtime-1)
  Y.P_vec <- rep(NA, simtime-1)
  G.A_vec <- rep(NA, simtime-1)
  G.S_vec <- rep(NA, simtime-1)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # growing season survival
    H.A <- H_A_fun(disease, A[t], S[t], P[t], iter)
    H.S <- H_S_fun(disease, A[t], S[t], P[t], iter)
    H.P <- H_P_fun(disease, A[t], S[t], P[t], iter)
    
    # non-growing season survival
    W.S <- W_S_fun(disease, A[t], S[t], P[t], iter)
    W.P <- W_P_fun(iter)
    
    # reduce growth due to competition
    V.A <- V_A_fun(disease, A[t], S[t], P[t], iter)
    V.S <- V_S_fun(disease, A[t], S[t], P[t], iter)
    V.P <- V_P_fun(disease, A[t], S[t], P[t], iter) 
    
    # seed production based on biomass
    Y.A <- Y_A_fun(disease, V.A, iter)
    Y.S <- Y_SP_fun(disease, V.S, iter)
    Y.P <- Y_SP_fun(disease, V.P, iter)
    
    # germination depends on disease and decreases due to litter
    G.A <- G_A_fun(disease, L[t], iter)
    litter_pres <- ifelse(L[t] > 0, 1, 0)
    G.S <- G_S_fun(disease, litter_pres, iter)
    
    # population size
    A[t+1] = s.A * (1-G.A) * A[t] + G.A * H.A * Y.A * A[t]  
    L[t+1] = G.A * H.A * V.A * A[t] + b * L[t]    
    S[t+1] = s.S * (1-G.S) * S[t] + G.S * H.S * Y.S * S[t] + H.P * Y.P * P[t]
    P[t+1] = H.P * W.P * P[t] + G.S * H.S * W.S * S[t]  
    
    # correct to prevent negative numbers
    A[t+1] = ifelse(A[t+1] < 1, 0, A[t+1])
    L[t+1] = ifelse(L[t+1] < 1, 0, L[t+1])
    S[t+1] = ifelse(S[t+1] < 1, 0, S[t+1])
    P[t+1] = ifelse(P[t+1] < 1, 0, P[t+1])
    
    # save parameter value
    H.A_vec[t] <- H.A
    H.S_vec[t] <- H.S
    H.P_vec[t] <- H.P
    W.S_vec[t] <- W.S
    W.P_vec[t] <- W.P
    V.A_vec[t] <- V.A
    V.S_vec[t] <- V.S
    V.P_vec[t] <- V.P
    Y.A_vec[t] <- Y.A
    Y.S_vec[t] <- Y.S
    Y.P_vec[t] <- Y.P
    G.A_vec[t] <- G.A
    G.S_vec[t] <- G.S
  }
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 4),
                species = rep(c("Elymus seedling", 
                                "Elymus adult", 
                                "Microstegium", 
                                "Microstegium litter"), 
                              each = simtime),
                N = c(S, P, A, L))

  # parameter data
  dfP <- tibble(time = rep(1:(simtime-1), 13),
                parameter = rep(c("H.A", "H.S", "H.P", "W.S", "W.P", "V.A", "V.S", "V.P", "Y.A", "Y.S", "Y.P", "G.A", "G.S"), 
                                each = simtime - 1),
                description = rep(c("annual growing season survival",
                                    "perennial seedling growing season survival",
                                    "perennial adult growing season survival",
                                    "perennial seedling non-growing season survival",
                                    "perennial adult non-growing season survival",
                                    "annual vigor",
                                    "perennial seedling vigor",
                                    "perennial adult vigor",
                                    "annual seed production",
                                    "perennial seedling seed production",
                                    "perennial adult seed production",
                                    "annual germination",
                                    "perennial germination"), 
                                  each = simtime - 1),
                value = c(H.A_vec, H.S_vec, H.P_vec, W.S_vec, W.P_vec, V.A_vec, V.S_vec, V.P_vec, Y.A_vec, Y.S_vec, Y.P_vec, G.A_vec, G.S_vec))

  # return
  return(list(dfN, dfP))
}




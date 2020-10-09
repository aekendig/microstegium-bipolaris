##### info ####

# file: 
# author: Amy Kendig
# date last edited: 9/30/20
# goal: simulate populations with data


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# number of iterations
n_samps = 10

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

# simulation time
simtime = 500

# initial conditions
A0 = 1 # initial annual population size
S0 = 1 # initial perennial seedling population size
P0 = 0 # initial perennial adult population size
L0 = 0 # initial annual litter amount


#### model ####

simFun = function(params, A0, S0, P0, L0, simtime, disease, iter){
  
  # initialize populations
  A <- rep(NA,simtime)
  S <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  
  A[1] <- A0
  S[1] <- S0
  P[1] <- P0
  L[1] <- L0
  
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

    #### start here ####
    
    # reduce germination due to litter
    G.A = g.A / (1 + beta.A * L[t])
    G.S = g.S / (1 + beta.S * L[t])
    
    # population size
    N.A[t+1] = s.A * (1-G.A) * N.A[t] + G.A * h.A * Y.A * N.A[t]  
    L[t+1] = a * G.A * h.A * V.A * N.A[t] + (1-b) * L[t]    
    N.S[t+1] = s.S * (1-G.S) * N.S[t] + G.S * h.S * Y.S * N.S[t] + h.P * Y.P * N.P[t]
    N.P[t+1] = h.P * w.P * N.P[t] + G.S * h.S * w.S * N.S[t]  
    
    # correct to prevent negative numbers
    N.A[t+1] = ifelse(N.A[t+1] < 1, 0, N.A[t+1])
    L[t+1] = ifelse(L[t+1] < 1, 0, L[t+1])
    N.S[t+1] = ifelse(N.S[t+1] < 1, 0, N.S[t+1])
    N.P[t+1] = ifelse(N.P[t+1] < 1, 0, N.P[t+1])
  }
  
  # save data
  dfN = data.frame(time = rep(1:simtime, 4), N = c(N.S, N.P, N.A, L), species = rep(c("Elymus seedling", "Elymus adult", "Microstegium", "Microstegium litter"), each = simtime))
  
  # return
  return(dfN)
}



# define parameters
g.A = filter(params, symbol == "g.A")$value
beta.A = filter(params, symbol == "beta.A")$value
s.A = filter(params, symbol == "s.A")$value
h.A = filter(params, symbol == "h.A")$value
v.A = filter(params, symbol == "v.A")$value
alpha.AA = filter(params, symbol == "alpha.AA")$value
alpha.AS = filter(params, symbol == "alpha.AS")$value
alpha.AP = filter(params, symbol == "alpha.AP")$value
y.A = filter(params, symbol == "y.A")$value
y.Aint = filter(params, symbol == "y.Aint")$value
a = filter(params, symbol == "a")$value
b = filter(params, symbol == "b")$value
g.S = filter(params, symbol == "g.S")$value
beta.S = filter(params, symbol == "beta.A")$value
s.S = filter(params, symbol == "s.S")$value
v.S = filter(params, symbol == "v.S")$value
v.P = filter(params, symbol == "v.P")$value
alpha.SA = filter(params, symbol == "alpha.SA")$value
alpha.SS = filter(params, symbol == "alpha.SS")$value
alpha.SP = filter(params, symbol == "alpha.SP")$value
alpha.PA = filter(params, symbol == "alpha.PA")$value
alpha.PS = filter(params, symbol == "alpha.PS")$value
alpha.PP = filter(params, symbol == "alpha.PP")$value
y.S = filter(params, symbol == "y.S")$value
y.Sint = filter(params, symbol == "y.Sint")$value
y.P = filter(params, symbol == "y.P")$value
y.Pint = filter(params, symbol == "y.Pint")$value
w.S = filter(params, symbol == "w.S")$value
w.P = filter(params, symbol == "w.P")$value
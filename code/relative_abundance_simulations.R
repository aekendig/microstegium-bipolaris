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

# load scripts
source("code/microstegium_gs_survival_parameter_2019_density_exp.R")
source("code/elymus_seedling_gs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_gs_survival_parameter_2019_density_exp.R")

source("code/elymus_seedling_ngs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_ngs_survival_parameter_2019_density_exp.R")

source("code/microstegium_biomass_parameter_2019_density_exp.R")
source("code/elymus_seedling_biomass_parameter_2019_density_exp.R")
source("code/elymus_adult_biomass_parameter_2019_density_exp.R")

# simulation time
simtime = 10000

# invasion time
invtime = 6000

# initial conditions
N0.A = 1 # initial annual population size
N0.S = 1 # initial perennial seedling population size
N0.P = 0 # initial perennial adult population size
L0 = 0 # initial annual litter amount
Ni.A = 1 # introduction of annual 
Ni.S = 1 # introduction of perennial


#### model ####

simFun = function(params, A0, S0, P0, L0, simtime, disease, iter){
  
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
    
    
    # reduce germination due to litter
    G.A = g.A / (1 + beta.A * L[t])
    G.S = g.S / (1 + beta.S * L[t])
    
    # growing season survival
    H.A <- H_A_fun(disease, A[t], S[t], P[t], iter)
    H.S <- H_S_fun(disease, A[t], S[t], P[t], iter)
    H.P <- H_P_fun(disease, A[t], S[t], P[t], iter)
    
    # reduce growth due to competition
    V.A = v.A / (1 + alpha.AA * G.A * N.A[t] + alpha.AP * N.P[t] + alpha.AS * G.S * N.S[t])    
    V.S = v.S / (1 + alpha.SA * G.A * N.A[t] + alpha.SP * N.P[t] + alpha.SS * G.S * N.S[t]) 
    V.P = v.P / (1 + alpha.PA * G.A * N.A[t] + alpha.PP * N.P[t] + alpha.PS * G.S * N.S[t])  
    
    # seed production based on biomass
    Y.A = V.A^y.A * y.Aint
    # Y.S = V.S^y.S * y.Sint 
    Y.S = V.S^y.S * y.Sint - 1
    Y.S = ifelse(Y.S < 0, 0, Y.S)
    # Y.P = V.P^y.P * y.Pint
    Y.P = V.P^y.P * y.Pint - 1
    Y.P = ifelse(Y.P < 0, 0, Y.P)
    
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
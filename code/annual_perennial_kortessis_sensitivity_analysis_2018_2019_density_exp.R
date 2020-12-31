##### info ####

# file: annual_perennial_kortessis_sensitivity_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/19/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# number of samples
# n_samps <- 10 # specify in scripts that source this one

# load packages
library(tidyverse)
library(brms)

# load scripts
source("code/microstegium_establishment_bh_parameter_2018_2019_density_exp.R")
source("code/elymus_seedling_establishment_bh_parameter_2018_2019_density_exp.R")

source("code/elymus_adult_gs_survival_parameter_2018_2019_density_exp.R")
source("code/elymus_adult_ngs_survival_parameter_2019_density_exp.R")

source("code/microstegium_biomass_fung_parameter_2018_2019_density_exp.R")
source("code/elymus_adult_biomass_fung_parameter_2019_density_exp.R")

source("code/microstegium_seed_production_parameter_2018_2019_density_exp.R")
source("code/elymus_seedling_seed_production_parameter_2018_2019_density_exp.R")
source("code/elymus_adult_seed_production_parameter_2018_2019_density_exp.R")

source("code/microstegium_germination_parameter_2018_density_exp.R")
source("code/elymus_germination_parameter_2018_2019_density_exp.R")

# constant parameters
s.A <- 0.15 # annual seed survival (Redwood et al. 2018)
s.P <- 0.05 # perennial seed survival (Garrison and Stier 2010)
d <- 0.61 # annual litter decomposition (DeMeester and Richter 2010)
h <- 0.29 # seedling survival from germination to establishment (Emery et al. 2013)


#### model ####

sens_fun = function(disease, params)

sim_fun = function(A0, S0, P0, L0, simtime, disease, iter, d_in){
  
  # initialize populations
  A <- rep(NA,simtime)
  S <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  
  A[1] <- A0
  S[1] <- S0
  P[1] <- P0
  L[1] <- L0
  
  # constant parameters
  
  # decomposition rate
  d <- d_in
  
  # perennial adult growing season survival
  u.P <- u_P_fun(disease, iter)
  
  # perennial non-growing season survival
  w.P <- w_P_fun(disease, iter)
  
  # perennial annual survival
  s <- ifelse(w.P * u.P > 0.99, 0.99, w.P * u.P)
  
  # biomass production
  b.A <- b_A_fung_fun(disease, iter)
  b.P <- b_P_fung_fun(disease, iter) 
  
  # germination
  g.A <- g_A_fun(disease, iter)
  g.P <- g_S_fun(disease, iter)
  
  # maximum establishment
  e.A <- ifelse(disease == 1, 
                exp(E_A_df$int_wat[iter])/(1 + exp(E_A_df$int_wat[iter])), 
                exp(E_A_df$int_fun[iter])/(1 + exp(E_A_df$int_fun[iter]))) * h
  e.P <- ifelse(disease == 1, 
                exp(E_S_df$int_wat[iter])/(1 + exp(E_S_df$int_wat[iter])), 
                exp(E_S_df$int_fun[iter])/(1 + exp(E_S_df$int_fun[iter]))) * h
  
  # maximum seed production
  y.A <- ifelse(disease == 1, 
                exp(Y_A_dens$int_wat[iter]) - 1, 
                exp(Y_A_dens$int_wat[iter] * Y_A_fun_eff$fun_eff[iter]) - 1)
  y.A <- ifelse(y.A < 0, 0, y.A)
  
  y.P <- ifelse(disease == 1, 
                exp(Y_P_dens$int_wat[iter]) - 1, 
                exp(Y_P_dens$int_wat[iter] * Y_P_fun_eff$fun_eff[iter]) - 1)
  y.P <- ifelse(y.P < 0, 0, y.P)
  
  y.S <- ifelse(disease == 1, 
                exp(Y_S_dens$int_wat[iter]) - 1, 
                exp(Y_S_dens$int_wat[iter] * Y_S_fun_eff$fun_eff[iter]) - 1)
  y.S <- ifelse(y.S < 0, 0, y.S)
  
  y.1 <- ifelse(y.P > 0, y.S / y.P, 0)
  
  # competition coefficients
  alpha.A.A <- ifelse(disease == 1, 
                      Y_A_dens$mv_dens_wat, 
                      Y_A_dens$mv_dens_fun)
  alpha.A.P <- ifelse(disease == 1, 
                      Y_A_dens$evA_dens_wat, 
                      Y_A_dens$evA_dens_fun)
  gamma.A <- ifelse(disease == 1, 
                    Y_A_dens$evS_dens_wat/alpha.A.P, 
                    Y_A_dens$evS_dens_fun/alpha.A.P)
  
  alpha.P.A <- ifelse(disease == 1, 
                      Y_P_dens$mv_dens_wat, 
                      Y_P_dens$mv_dens_fun)
  alpha.P.P <- ifelse(disease == 1, 
                      Y_P_dens$evA_dens_wat, 
                      Y_P_dens$evA_dens_fun)
  gamma.P <- ifelse(disease == 1, 
                    Y_P_dens$evS_dens_wat/alpha.P.P, 
                    Y_P_dens$evS_dens_fun/alpha.P.P)
  
  alpha.S.A <- ifelse(disease == 1, 
                      Y_S_dens$mv_dens_wat, 
                      Y_S_dens$mv_dens_fun)
  alpha.S.P <- ifelse(disease == 1, 
                      Y_S_dens$evA_dens_wat, 
                      Y_S_dens$evA_dens_fun)
  gamma.S <- ifelse(disease == 1, 
                    Y_S_dens$evS_dens_wat/alpha.S.P, 
                    Y_S_dens$evS_dens_fun/alpha.S.P)
  
  alpha.A <- mean(c(alpha.A.A, alpha.P.A, alpha.S.A))
  alpha.P <- mean(c(alpha.A.P, alpha.P.P, alpha.S.P))
  gamma <- mean(c(gamma.A, gamma.P, gamma.S))
  
  alpha.i.P <- alpha.P * (1 + gamma * (1 - s))
  
  # litter effect values
  beta.A <- ifelse(disease == 1,
                   E_A_beta_wat[iter],
                   E_A_beta_fun$beta[iter])
  
  beta.P <- ifelse(disease == 1,
                   E_S_beta_wat[iter],
                   E_S_beta_fun$beta[iter])
  
  # lambda values
  l.A <- g.A * e.A * y.A / (1 - s.A * (1 - g.A))
  l.P <- y.P * e.P * g.P * (y.1 + 1 / (1 - s))/(1 - s.P * (1 - g.P))
  
  # equilibrium litter density
  L.P.1 = -0.5 * (b.P / (d * alpha.i.P) + 1 / beta.P)
  L.P.2 = sqrt(0.25 * (b.P / (d * alpha.i.P) + 1 / beta.P)^2 + (b.P / d) * (l.P - 1) / (alpha.i.P * beta.P))
  L.P = ifelse((L.P.1 + L.P.2) > (L.P.1 - L.P.2), L.P.1 + L.P.2, L.P.1 - L.P.2)
  
  L.A.1 = -0.5 * (b.A / (d * alpha.A) + 1 / beta.A)
  L.A.2 = sqrt(0.25 * (b.A / (d * alpha.A) + 1 / beta.A)^2 + (b.A / d) * (l.A - 1) / (alpha.A * beta.A))
  L.A = ifelse((L.A.1 + L.A.2) > (L.A.1 - L.A.2), L.A.1 + L.A.2, L.A.1 - L.A.2)
  
  # initialize parameter vectors
  E.A_vec <- rep(NA, simtime-1)
  E.P_vec <- rep(NA, simtime-1)
  e.A_vec <- rep(e.A, simtime-1)
  e.P_vec <- rep(e.P, simtime-1)
  s_vec <- rep(s, simtime-1)
  b.A_vec <- rep(b.A, simtime-1)
  b.P_vec <- rep(b.P, simtime-1)
  Y.A_vec <- rep(NA, simtime-1)
  y.1_vec <- rep(y.1, simtime-1)
  Y.P_vec <- rep(NA, simtime-1)
  g.A_vec <- rep(g.A, simtime-1)
  g.P_vec <- rep(g.P, simtime-1)
  l.A_vec <- rep(l.A, simtime-1)
  l.P_vec <- rep(l.P, simtime-1)
  L.A_vec <- rep(L.A, simtime-1)
  L.P_vec <- rep(L.P, simtime-1)
  beta.A_vec <- rep(beta.A, simtime-1)
  beta.P_vec <- rep(beta.P, simtime-1)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # seedling establishment
    E.A <- E_A_bh_fun(disease, L[t], iter) * h
    E.P <- E_S_bh_fun(disease, L[t], iter) * h
    
    # reduce seed production due to competition
    Y.A <- y.A / (1 + alpha.A * g.A * E.A * A[t] + alpha.P * P[t] + alpha.P * gamma * g.P * E.P * S[t])
    Y.P <- y.P / (1 + alpha.A * g.A * E.A * A[t] + alpha.P * P[t] + alpha.P * gamma * g.P * E.P * S[t])
    
    # population size
    A[t+1] = s.A * (1 - g.A) * A[t] + g.A * E.A * Y.A * A[t]  
    L[t+1] = g.A * E.A * b.A * A[t] + b.P * P[t] + (1 - d) * L[t]    
    S[t+1] = s.P * (1 - g.P) * S[t] + g.P * E.P * Y.P * y.1 * S[t] + Y.P * P[t]
    P[t+1] = s * P[t] + g.P * E.P * S[t]  
    
    # correct to prevent negative numbers
    A[t+1] = ifelse(A[t+1] < 0, 0, A[t+1])
    L[t+1] = ifelse(L[t+1] < 0, 0, L[t+1])
    S[t+1] = ifelse(S[t+1] < 0, 0, S[t+1])
    P[t+1] = ifelse(P[t+1] < 0, 0, P[t+1])
    
    # save parameter value
    E.A_vec[t] <- E.A
    E.P_vec[t] <- E.P
    Y.A_vec[t] <- Y.A
    Y.P_vec[t] <- Y.P
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
  dfP <- tibble(time = rep(1:(simtime-1), 18),
                parameter = rep(c("E.A", "E.P", "e.A", "e.P", "s", "b.A", "b.P", "Y.A", "y.1", "Y.P", "g.A", "g.P", "l.A", "l.P", "L.A", "L.P", "beta.A", "beta.P"), 
                                each = simtime - 1),
                value = c(E.A_vec, E.P_vec, e.A_vec, e.P_vec, s_vec, b.A_vec, b.P_vec, Y.A_vec, y.1_vec, Y.P_vec, g.A_vec, g.P_vec, l.A_vec, l.P_vec, L.A_vec, L.P_vec, beta.A_vec, beta.P_vec),
                iteration = iter)

  # return
  return(list(dfN, dfP))
}




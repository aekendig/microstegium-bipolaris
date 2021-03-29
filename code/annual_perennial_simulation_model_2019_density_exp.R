##### info ####

# file: annual_perennial_simulation_model_2019_density_exp
# author: Amy Kendig
# date last edited: 3/28/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# number of samples
# n_samps <- 10 # specify in scripts that source this one

# load packages
library(tidyverse)
library(brms)

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


# posterior samples
u.P.samps <- posterior_samples(evASurvD2Mod) %>%
  transmute(u.P = exp(b_Intercept)/(1 + exp(b_Intercept)))

w.S.samps <- posterior_samples(evSWinSurvD1Mod) %>%
  transmute(w.S = exp(b_Intercept)/(1 + exp(b_Intercept)))
w.P.samps <- posterior_samples(evAWinSurvD1Mod) %>%
  transmute(w.P = exp(b_Intercept)/(1 + exp(b_Intercept)))

g.A.samps <- posterior_samples(mvGermD1Mod) %>%
  transmute(g.A = exp(b_Intercept)/(1 + exp(b_Intercept)))
g.S.samps <- posterior_samples(evGermMod) %>%
  transmute(g.S = exp(b_Intercept)/(1 + exp(b_Intercept)))

e.A.samps <- posterior_samples(mvSurvD2Mod) %>%
  transmute(e.A = exp(b_Intercept)/(1 + exp(b_Intercept)))
e.S.samps <- posterior_samples(evSSurvD2Mod) %>%
  transmute(e.S = exp(b_Intercept)/(1 + exp(b_Intercept)))

beta.A.samps <- posterior_samples(mvLitEstBhMod2) %>%
  transmute(beta.A = b_betaL_Intercept)
beta.S.samps <- posterior_samples(evLitEstBhMod2) %>%
  transmute(beta.S = b_betaL_Intercept)

b.A.samps <- posterior_samples(mvMvD2Mod) %>%
  transmute(b.A = b_b0_treatmentwater)
b.S.samps <- posterior_samples(evSMvD2Mod) %>%
  transmute(b.S = b_b0_treatmentwater)
b.P.samps <- posterior_samples(evAMvD2Mod) %>%
  transmute(b.P = b_b0_treatmentwater)

alpha.AA.samps <- posterior_samples(mvMvD2Mod) %>%
  transmute(alpha.AA = b_alpha_treatmentwater)
alpha.SA.samps <- posterior_samples(evSMvD2Mod) %>%
  transmute(alpha.SA = b_alpha_treatmentwater)
alpha.PA.samps <- posterior_samples(evAMvD2Mod) %>%
  transmute(alpha.PA = b_alpha_treatmentwater)

alpha.AS.samps <- posterior_samples(mvEvSD2Mod) %>%
  transmute(alpha.AS = b_alpha_treatmentwater)
alpha.SS.samps <- posterior_samples(evSEvSD2Mod) %>%
  transmute(alpha.SS = b_alpha_treatmentwater)
alpha.PS.samps <- posterior_samples(evAEvSD2Mod) %>%
  transmute(alpha.PS = b_alpha_treatmentwater)

alpha.AP.samps <- posterior_samples(mvEvAD2Mod) %>%
  transmute(alpha.AP = b_alpha_treatmentwater)
alpha.SP.samps <- posterior_samples(evSEvAD2Mod) %>%
  transmute(alpha.SP = b_alpha_treatmentwater)
alpha.PP.samps <- posterior_samples(evAEvAD2Mod) %>%
  transmute(alpha.PP = b_alpha_treatmentwater)

y.A.samps <- posterior_samples(mvSeedsBioD2Mod) %>%
  transmute(y.A = b_Intercept)
y.S.samps <- posterior_samples(evSSeedsBioD2Mod) %>%
  transmute(y.S = b_Intercept)
y.P.samps <- posterior_samples(evASeedsBioD2Mod) %>%
  transmute(y.P = b_Intercept)

yb.A.samps <- posterior_samples(mvSeedsBioD2Mod) %>%
  transmute(yb.A = b_log_bio)
yb.S.samps <- posterior_samples(evSSeedsBioD2Mod) %>%
  transmute(yb.S = b_log_bio)
yb.P.samps <- posterior_samples(evASeedsBioD2Mod) %>%
  transmute(yb.P = b_log_bio)

# constant parameters
s.A <- s.P <- 0.1 # annual seed survival 
# median of Redwood et al. 2018, Mv = 0.15; Garrison and Stier 2010, Ev = 0.05
d <- 0.35 # litter decomposition 
# median of DeMeester and Richter 2010, Mv = 0.6; Hobbie et al. 2008, min = 0.1


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
  
  # perennial adult growing season survival
  u.P <- u.P.samps$u.P[iter]
  
  # perennial non-growing season survival
  w.S <- w.P.samps$w.S[iter]
  w.P <- w.P.samps$w.P[iter]
  
  # perennial lifespan
  l.P <- ifelse(w.P * u.P > 0.99, 0.99, w.P * u.P)
  
  # germination
  g.A <- g.A.samps$g.A[iter]
  g.S <- ifelse(disease == T, g.S.samps$g.S[iter], g.S.samps$g.S[iter] * 0.7)
  
  # maximum establishment
  e.A <- e.A.samps$e.A[iter]
  e.S <- e.S.samps$e.S[iter]
  
  # litter effect
  beta.A <- beta.A.samps$beta.A[iter]
  beta.S <- beta.S.samps$beta.S[iter]
  
  # maximum growth
  b.A <- b.A.samps$b.A[iter]
  b.S <- b.S.samps$b.S[iter]
  b.P <- b.P.samps$b.P[iter]
  
  # competitive effect of annual
  alpha.AA <- alpha.AA.samps$alpha.AA[iter]
  alpha.SA <- ifelse(disease == T, alpha.SA.samps$alpha.SA[iter], alpha.SA.samps$alpha.SA[iter] * 0.11)
  alpha.PA <- alpha.PA.samps$alpha.PA[iter]
  
  # competitive effect of perennial seedling
  alpha.AS <- alpha.AS.samps$alpha.AS[iter]
  alpha.SS <- alpha.SS.samps$alpha.SS[iter]
  alpha.PS <- alpha.PS.samps$alpha.PS[iter]
  
  # competitive effect of perennial adult
  alpha.AP <- alpha.AP.samps$alpha.AP[iter]
  alpha.SP <- alpha.SP.samps$alpha.SP[iter]
  alpha.PP <- alpha.PP.samps$alpha.PP[iter]
  
  # minimum seeds
  y.A <- y.A.samps$y.A[iter]
  y.S <- y.S.samps$y.S[iter]
  y.P <- y.P.samps$y.P[iter]
  
  # biomass-seed conversion
  yb.A <- yb.A.samps$yb.A[iter]
  yb.S <- yb.S.samps$yb.S[iter]
  yb.P <- yb.P.samps$yb.P[iter]
  
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




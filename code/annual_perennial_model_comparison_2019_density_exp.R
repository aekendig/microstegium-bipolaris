##### info ####

# file: annual_perennial_model_comparison_2019_density_exp
# author: Amy Kendig
# date last edited: 10/29/20
# goal: compare Kortessis model to experimental model


#### set-up ####

# clear environment
rm(list = ls())

# number of samples
n_samps <- 300

# load packages
library(tidyverse)
library(brms)

# load scripts
source("code/microstegium_establishment_bh_parameter_2019_density_exp.R")
source("code/elymus_seedling_establishment_bh_parameter_2019_density_exp.R")
source("code/microstegium_establishment_parameter_2019_density_exp.R")
source("code/elymus_seedling_establishment_parameter_2019_density_exp.R")

source("code/elymus_adult_gs_survival_parameter_2019_density_exp.R")
source("code/elymus_adult_ngs_survival_parameter_2019_density_exp.R")

source("code/elymus_seedling_ngs_survival_parameter_2019_density_exp.R")

source("code/microstegium_biomass_fung_parameter_2019_density_exp.R")
source("code/elymus_adult_biomass_fung_parameter_2019_density_exp.R")

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
s.P <- 0.05 # perennial seed survival (Garrison and Stier 2010)
s.S <- s.P
d <- 0.61 # annual litter decomposition (DeMeester and Richter 2010)


#### Kortessis model ####

k_sim_fun = function(A0, S0, P0, L0, simtime, disease, iter, d_in){
  
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
                exp(E_A_df$int_fun[iter])/(1 + exp(E_A_df$int_fun[iter])))
  e.P <- ifelse(disease == 1, 
                exp(E_S_df$int_wat[iter])/(1 + exp(E_S_df$int_wat[iter])), 
                exp(E_S_df$int_fun[iter])/(1 + exp(E_S_df$int_fun[iter])))
  
  # maximum seed production
  y.A <- ifelse(disease == 1, 
                Y_A_dens$int_wat[iter], 
                Y_A_dens$int_fun[iter])
  y.P <- u.P * ifelse(disease == 1, 
                      Y_P_dens$int_wat[iter], 
                      Y_P_dens$int_fun[iter])
  y.S <- ifelse(disease == 1, 
                Y_S_dens$int_wat[iter], 
                Y_S_dens$int_fun[iter])
  y.1 <- y.S / y.P
  
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
    E.A <- E_A_bh_fun(disease, L[t], iter)
    E.P <- E_S_bh_fun(disease, L[t], iter)
    
    # reduce seed production due to competition
    Y.A <- y.A / (1 + alpha.A * g.A * E.A * A[t] + alpha.P * P[t] + alpha.P * gamma * g.P * E.P * S[t])
    Y.P <- y.P / (1 + alpha.A * g.A * E.A * A[t] + alpha.P * P[t] + alpha.P * gamma * g.P * E.P * S[t])
    
    # population size
    A[t+1] = s.A * (1 - g.A) * A[t] + g.A * E.A * Y.A * A[t]  
    L[t+1] = g.A * E.A * b.A * A[t] + b.P * P[t] + (1 - d) * L[t]    
    S[t+1] = s.P * (1 - g.P) * S[t] + g.P * E.P * Y.P * y.1 * S[t] + Y.P * P[t]
    P[t+1] = s * P[t] + g.P * E.P * S[t]  
    
    # correct to prevent negative numbers
    A[t+1] = ifelse(A[t+1] < 1, 0, A[t+1])
    L[t+1] = ifelse(L[t+1] < 1, 0, L[t+1])
    S[t+1] = ifelse(S[t+1] < 1, 0, S[t+1])
    P[t+1] = ifelse(P[t+1] < 1, 0, P[t+1])
    
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


#### experimental model model ####

e_sim_fun = function(A0, S0, P0, L0, simtime, disease, iter){
  
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


#### simulation conditions ####

# simulation time
simtimeR <- 600
simtimeI <- 1000

# samples from parameter distributions
samps <- n_samps


#### resident E. virginicus ####

# initial conditions
A0E <- 0 # initial annual population size
S0E <- 0 # initial perennial seedling population size
P0E <- 1 # initial perennial adult population size
L0E <- 0 # initial annual litter amount

# initiate lists
k_abundE <- list()
k_abundE2 <- list()
e_abundE <- list()
e_abundE2 <- list()

# simulations
for(iter in 1:samps){
  
  # initial Elymus dynamics
  k_mod_estE <- k_sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 0, iter, d)
  e_mod_estE <- e_sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 0, iter)
  k_mod_estE2 <- k_sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 1, iter, d)
  e_mod_estE2 <- e_sim_fun(A0E, S0E, P0E, L0E, simtimeR, disease = 1, iter)
  
  # new initial conditions
  k_S0Ei <- k_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Perennial seedling") %>% 
    pull(N)
  k_P0Ei <- k_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Perennial adult") %>% 
    pull(N)
  k_L0Ei <- k_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Litter") %>% 
    pull(N)
  
  e_S0Ei <- e_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Perennial seedling") %>% 
    pull(N)
  e_P0Ei <- e_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Perennial adult") %>% 
    pull(N)
  e_L0Ei <- e_mod_estE[[1]] %>% 
    filter(time == simtimeR & species == "Litter") %>% 
    pull(N)
  
  k_S0Ei2 <- k_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Perennial seedling") %>% 
    pull(N)
  k_P0Ei2 <- k_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Perennial adult") %>% 
    pull(N)
  k_L0Ei2 <- k_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Litter") %>% 
    pull(N)
  
  e_S0Ei2 <- e_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Perennial seedling") %>% 
    pull(N)
  e_P0Ei2 <- e_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Perennial adult") %>% 
    pull(N)
  e_L0Ei2 <- e_mod_estE2[[1]] %>% 
    filter(time == simtimeR & species == "Litter") %>% 
    pull(N)
  
  # invasion simulation
  k_mod_invE <- k_sim_fun(A0 = 10, k_S0Ei, k_P0Ei, k_L0Ei, simtimeI, disease = 0, iter, d)
  e_mod_invE <- e_sim_fun(A0 = 10, e_S0Ei, e_P0Ei, e_L0Ei, simtimeI, disease = 0, iter)
  k_mod_invE2 <- k_sim_fun(A0 = 10, k_S0Ei2, k_P0Ei2, k_L0Ei2, simtimeI, disease = 1, iter, d)
  e_mod_invE2 <- e_sim_fun(A0 = 10, e_S0Ei2, e_P0Ei2, e_L0Ei2, simtimeI, disease = 1, iter)
  
  # save abundances
  k_abundE[[iter]] <- k_mod_invE[[1]]
  e_abundE[[iter]] <- e_mod_invE[[1]]
  k_abundE2[[iter]] <- k_mod_invE2[[1]]
  e_abundE2[[iter]] <- e_mod_invE2[[1]]
  
}

# convert to dataframes
k_abund_datE <- do.call(rbind, k_abundE)
k_abund_datE2 <- do.call(rbind, k_abundE2)
e_abund_datE <- do.call(rbind, e_abundE)
e_abund_datE2 <- do.call(rbind, e_abundE2)


#### process data ####

rel_abund_dat <- k_abund_datE %>%
  mutate(disease = "without disease",
         model = "Kortessis et al.") %>%
  full_join(k_abund_datE2 %>%
              mutate(disease = "with disease",
                     model = "Kortessis et al.")) %>%
  full_join(e_abund_datE %>%
              mutate(disease = "without disease",
                     model = "experiment")) %>%
  full_join(e_abund_datE2 %>%
              mutate(disease = "with disease",
                     model = "experiment")) %>%
  pivot_wider(names_from = species, values_from = N) %>%
  rename_with(~ gsub(" ", "_", .x, fixed = T)) %>%
  mutate(Perennial = Perennial_seedling + Perennial_adult,
         tot_abund = Annual + Perennial,
         ann_rel_abund = Annual / tot_abund,
         per_rel_abund = Perennial / tot_abund,
         iterationF = as.factor(iteration))

# relative abundance summary
rel_sum_dat <- rel_abund_dat %>%
  filter(time >= simtimeI - 100) %>%
  group_by(model, disease, iteration, iterationF) %>%
  summarise(ann_rel_abund = mean(ann_rel_abund),
            per_rel_abund = mean(per_rel_abund)) %>%
  ungroup() %>%
  mutate(disease = fct_relevel(disease, "without disease"))


#### figure ####

pdf("output/annual_perennial_model_comparison_2019_density_exp.pdf", width = 8, height = 3)
ggplot(rel_sum_dat, aes(x = ann_rel_abund, fill = disease)) +
  geom_histogram(binwidth = 0.05, position = position_dodge()) +
  facet_wrap(~ model) +
  scale_fill_viridis_d(direction = -1, end = 0.6) +
  xlab(expression(paste("Relative abundance (", italic("M. vimineum"), "/", italic("E. virginicus"), ")", sep = ""))) +
  ylab("Number of draws") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10),
        legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10),
        strip.background = element_blank())
dev.off()


#### output ####
write_csv(rel_abund_dat, "output/annual_perennial_model_comparison_2019_density_exp.csv")

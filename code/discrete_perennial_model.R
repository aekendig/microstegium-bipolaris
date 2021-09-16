disc_perennial_mod <- function(F0, P0, L0, C0, IF0, IP0, BP0, BF0, simtime, gs_days, disc_parms, con_parms){
  
  # initialize populations
  F1 <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  B_F <- rep(NA,simtime)
  B_P <- rep(NA,simtime)
  S_F <- rep(NA,simtime)
  S_P <- rep(NA,simtime)
  I_F <- rep(NA,simtime)
  I_P <- rep(NA,simtime)
  D <- rep(NA,simtime)
  C_0 <- rep(NA,simtime)
  I_F0 <- rep(NA,simtime)
  I_P0 <- rep(NA,simtime)
  B_P0 <- rep(NA,simtime)
  B_F0 <- rep(NA,simtime)
  
  F1[1] <- F0
  P[1] <- P0
  L[1] <- L0
  C_0[1] <- C0
  I_F0[1] <- IF0
  I_P0[1] <- IP0
  B_P0[1] <- BP0
  B_F0[1] <- BF0
  
  # parameters
  parms <- disc_parms
  
  s_P <- as.numeric(parms[parms$Parameter == "s_P", "Estimate"])
  g_P <- as.numeric(parms[parms$Parameter == "g_P", "Estimate"])
  c_F <- as.numeric(parms[parms$Parameter == "c_F", "Estimate"])
  c_P <- as.numeric(parms[parms$Parameter == "c_P", "Estimate"])
  b_F <- as.numeric(parms[parms$Parameter == "b_F", "Estimate"])
  b_P <- as.numeric(parms[parms$Parameter == "b_P", "Estimate"])
  e_P <- as.numeric(parms[parms$Parameter == "e_P", "Estimate"])
  gamma_P <- as.numeric(parms[parms$Parameter == "gamma_P", "Estimate"])
  l_P <- as.numeric(parms[parms$Parameter == "l_P", "Estimate"])
  l_B <- as.numeric(parms[parms$Parameter == "l_B", "Estimate"])
  
  d <- as.numeric(parms[parms$Parameter == "d", "Estimate"])
  h <- as.numeric(parms[parms$Parameter == "h", "Estimate"])
  
  grow_days <- seq(0, gs_days)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # seedling establishment
    E_P <- e_P/(1 + gamma_P * L[t])
    
    # new adult growth establishment
    L_B <- l_B/(1 + gamma_P * L[t])
    
    # initial conditions, prevent taking log of negative values or zero
    LogB_F0 <- log(g_P * E_P * b_F * F1[t])
    LogB_P0 <- log(L_B * (B_P0[t] + B_F0[t]))
    
    # biomass at the end of the growing season
    if(F1[t] > 0 | (B_P0[t] + B_F0[t]) > 0){
      out <- ode(func = cont_perennial_mod,
                 y = c(LogB_F = LogB_F0,
                       LogB_P = LogB_P0,
                       I_F = 0, I_P = 0, D = 0,
                       C = C_0[t] + h * (I_F0[t] + I_P0[t])),
                 times = grow_days,
                 parms = con_parms) %>%
        as.data.frame()
    } else{
      print("perennial gone")
    }
    
    # modify output
    out2 <- out %>%
      as_tibble() %>%
      filter(time == gs_days) %>%
      mutate(B_F = exp(LogB_F),
             B_P = exp(LogB_P),
             S_F = B_F - I_F,
             S_P = B_P - I_P) %>%
      select(-c(LogB_F, LogB_P))
    
    # save biomass
    B_F[t] <- as.numeric(out2$B_F)
    B_P[t] <- as.numeric(out2$B_P)
    
    S_F[t] <- as.numeric(out2$S_F)
    S_P[t] <- as.numeric(out2$S_P)
    
    I_F[t] <- as.numeric(out2$I_F)
    I_P[t] <- as.numeric(out2$I_P)
    
    D[t] <- as.numeric(out2$D)
    
    # correct to prevent negative/crazy small numbers
    B_F[t] = ifelse(B_F[t] < 1e-10, 0, B_F[t])
    B_P[t] = ifelse(B_P[t] < 1e-10, 0, B_P[t])
    
    S_F[t] = ifelse(S_F[t] < 1e-10, 0, S_F[t])
    S_P[t] = ifelse(S_P[t] < 1e-10, 0, S_P[t])
    
    I_F[t] = ifelse(I_F[t] < 1e-10, 0, I_F[t])
    I_P[t] = ifelse(I_P[t] < 1e-10, 0, I_P[t])
    
    D[t] = ifelse(D[t] < 1e-10, 0, D[t])
    
    # values for next year
    C_0[t+1] <- as.numeric(out2$C) # conidia
    C_0[t+1] = ifelse(C_0[t+1] < 1e-10, 0, C_0[t+1])
    I_F0[t+1] <- I_F[t] # conidia
    I_P0[t+1] <- I_P[t] # conidia
    B_P0[t+1] <- B_P[t] # perennial biomass
    B_F0[t+1] <- B_F[t] # perennial biomass
    
    # population size
    L[t+1] = B_F[t] + B_P[t] + D[t] + (1 - d) * L[t] 
    F1[t+1] = s_P * (1 - g_P) * F1[t] + c_F * B_F[t] + c_P * B_P[t]
    P[t+1] = l_P * P[t] + g_P * E_P * l_P * F1[t]
    
    # correct to prevent negative numbers
    L[t+1] = ifelse(L[t+1] < 0, 0, L[t+1])
    F1[t+1] = ifelse(F1[t+1] < 0, 0, F1[t+1])
    P[t+1] = ifelse(P[t+1] < 0, 0, P[t+1])
    
  } 
  
  # last growing season
  # seedling establishment
  E_P <- e_P/(1 + gamma_P * L[simtime])
  
  # new adult growth establishment
  L_B <- l_B/(1 + gamma_P * L[simtime])
  
  # initial conditions
  LogB_F0 <- log(g_P * E_P * b_F * F1[simtime])
  LogB_P0 <- log(L_B * (B_P0[simtime] + B_F0[simtime]))
  
  # biomass at the end of the growing season
  if(F1[simtime] > 0 | (B_P0[simtime] + B_F0[simtime]) > 0){
    out3 <- ode(func = cont_perennial_mod,
                y = c(LogB_F = LogB_F0,
                      LogB_P = LogB_P0,
                      I_F = 0, I_P = 0, D = 0,
                      C = C_0[simtime] + h * (I_F0[simtime] + I_P0[simtime])),
                times = grow_days,
                parms = con_parms) %>%
      as.data.frame()
  } else{
    print("perennial gone")
  }
  
  # modify output
  out4 <- out3 %>%
    as_tibble() %>%
    filter(time == gs_days) %>%
    mutate(B_F = exp(LogB_F),
           B_P = exp(LogB_P),
           S_F = B_F - I_F,
           S_P = B_P - I_P) %>%
    select(-c(LogB_F, LogB_P))
  
  # save biomass
  B_F[simtime] <- as.numeric(out4$B_F)
  B_P[simtime] <- as.numeric(out4$B_P)
  
  S_F[simtime] <- as.numeric(out4$S_F)
  S_P[simtime] <- as.numeric(out4$S_P)
  
  I_F[simtime] <- as.numeric(out4$I_F)
  I_P[simtime] <- as.numeric(out4$I_P)
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 10),
                species = rep(c("Perennial first-year", 
                                "Perennial adult", 
                                "Litter",
                                "Perennial first-year biomass",
                                "Perennial adult biomass",
                                "Perennial first-year susceptible biomass",
                                "Perennial adult susceptible biomass",
                                "Perennial first-year infected biomass",
                                "Perennial adult infected biomass",
                                "Conidia"), 
                              each = simtime),
                N = c(F1, P, L, B_F, B_P, S_F, S_P, I_F, I_P, C_0))
  
  # return
  return(dfN)
}
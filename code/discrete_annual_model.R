disc_annual_mod <- function(AS0, AI0, L0, C0, IA0, simtime, gs_days, disc_parms, con_parms){
  
  # initialize populations
  A_S <- rep(NA,simtime)
  A_I <- rep(NA,simtime)
  L <- rep(NA,simtime)
  B_A <- rep(NA,simtime)
  S_A <- rep(NA,simtime)
  I_A <- rep(NA,simtime)
  D <- rep(NA,simtime)
  C_0 <- rep(NA,simtime)
  I_A0 <- rep(NA,simtime)
  
  A_S[1] <- AS0
  A_I[1] <- AI0
  L[1] <- L0
  C_0[1] <- C0
  I_A0[1] <- IA0
  
  # parameters
  parms <- disc_parms
  
  s_A <- as.numeric(parms[parms$Parameter == "s_A", "Estimate"])
  g_S <- as.numeric(parms[parms$Parameter == "g_S", "Estimate"])
  g_I <- as.numeric(parms[parms$Parameter == "g_I", "Estimate"])
  p0 <- as.numeric(parms[parms$Parameter == "p0", "Estimate"])
  p1 <- as.numeric(parms[parms$Parameter == "p1", "Estimate"])
  c_A <- as.numeric(parms[parms$Parameter == "c_A", "Estimate"])
  b_A <- as.numeric(parms[parms$Parameter == "b_A", "Estimate"])
  e_A <- as.numeric(parms[parms$Parameter == "e_A", "Estimate"])
  gamma_A <- as.numeric(parms[parms$Parameter == "gamma_A", "Estimate"])
  
  d <- as.numeric(parms[parms$Parameter == "d", "Estimate"])
  h <- as.numeric(parms[parms$Parameter == "h", "Estimate"])
  
  grow_days <- seq(0, gs_days)
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # seedling establishment
    E_A <- e_A/(1 + gamma_A * L[t])
    
    # initial conditions
    LogB_A0 <- log(g_S * E_A * b_A * A_S[t] + g_I * E_A * b_A * A_I[t])
    
    # biomass at the end of the growing season
    if((A_S[t] + A_I[t]) > 0){
      out <- ode(func = cont_annual_mod,
                 y = c(LogB_A = LogB_A0,
                       I_A = 0, D = 0,
                       C = C_0[t] + h * I_A0[t]),
                 times = grow_days,
                 parms = con_parms) %>%
        as.data.frame()
    }else{
      print("annual gone")
    }
    
    
    # modify output
    out2 <- out %>%
      as_tibble() %>%
      filter(time == gs_days) %>%
      mutate(B_A = exp(LogB_A),
             S_A = B_A - I_A) %>%
      select(-LogB_A)
    
    # save biomass
    B_A[t] <- as.numeric(out2$B_A)
    S_A[t] <- as.numeric(out2$S_A)
    I_A[t] <- as.numeric(out2$I_A)
    D[t] <- as.numeric(out2$D)
    
    # correct to prevent negative/crazy small numbers
    B_A[t] = ifelse(B_A[t] < 1e-10, 0, B_A[t])
    S_A[t] = ifelse(S_A[t] < 1e-10, 0, S_A[t])
    I_A[t] = ifelse(I_A[t] < 1e-10, 0, I_A[t])
    D[t] = ifelse(D[t] < 1e-10, 0, D[t])
    
    # values for next year
    C_0[t+1] <- as.numeric(out2$C) # conidia
    C_0[t+1] = ifelse(C_0[t+1] < 1e-10, 0, C_0[t+1])
    I_A0[t+1] <- I_A[t] # conidia
    
    # proportion seeds infected
    p <- ifelse(B_A[t] > 0, exp(p0 + p1 * I_A[t]/B_A[t])/(1 + exp(p0 + p1 * I_A[t]/B_A[t])), 0)
    
    # population size
    A_S[t+1] = s_A * (1 - g_S) * A_S[t] +  B_A[t] * c_A * (1-p)
    A_I[t+1] = B_A[t] * c_A * p
    L[t+1] = B_A[t] + D[t] + (1 - d) * L[t] 
    
    # correct to prevent negative numbers
    A_S[t+1] = ifelse(A_S[t+1] < 0, 0, A_S[t+1])
    A_I[t+1] = ifelse(A_I[t+1] < 0, 0, A_I[t+1])
    L[t+1] = ifelse(L[t+1] < 0, 0, L[t+1])
    
  } 
  
  # last growing season
  # seedling establishment
  E_A <- e_A/(1 + gamma_A * L[simtime])
  
  # initial conditions
  LogB_A0 <- log(g_S * E_A * b_A * A_S[simtime] + g_I * E_A * b_A * A_I[simtime])
  
  # biomass at the end of the growing season
  if((A_S[simtime] + A_I[simtime]) > 0){
    out3 <- ode(func = cont_annual_mod,
                y = c(LogB_A = LogB_A0,
                      I_A = 0, D = 0,
                      C = C_0[simtime] + h * I_A0[simtime]),
                times = grow_days,
                parms = con_parms) %>%
      as.data.frame()
  }else{
    print("annual gone")
  }
  
  
  # modify output
  out4 <- out3 %>%
    as_tibble() %>%
    filter(time == gs_days) %>%
    mutate(B_A = exp(LogB_A),
           S_A = B_A - I_A) %>%
    select(-LogB_A)
  
  # save biomass
  B_A[simtime] <- as.numeric(out4$B_A)
  S_A[simtime] <- as.numeric(out4$S_A)
  I_A[simtime] <- as.numeric(out4$I_A)
  
  # total annual density
  A_T <- A_S + A_I
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 8),
                species = rep(c("Susceptible annual",
                                "Infected annual",
                                "Annual",
                                "Litter",
                                "Annual biomass",
                                "Annual susceptible biomass",
                                "Annual infected biomass",
                                "Conidia"), 
                              each = simtime),
                N = c(A_S, A_I, A_T, L, B_A, S_A, I_A, C_0))
  
  # return
  return(dfN)
}
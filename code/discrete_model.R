disc_AFP_mod <- function(A0, F0, P0, L0, C0, BA0, BP0, BF0, simtime, gs_time, disc_parms, cont_parms){
  
  # initialize populations
  A <- rep(NA,simtime)
  F1 <- rep(NA,simtime)
  P <- rep(NA,simtime)
  L <- rep(NA,simtime)
  B_A <- rep(NA,simtime)
  B_F <- rep(NA,simtime)
  B_P <- rep(NA,simtime)
  S_A <- rep(NA,simtime)
  S_F <- rep(NA,simtime)
  S_P <- rep(NA,simtime)
  I_A <- rep(NA,simtime)
  I_F <- rep(NA,simtime)
  I_P <- rep(NA,simtime)
  D <- rep(NA,simtime)
  C <- rep(NA,simtime)
  
  A[1] <- A0
  F1[1] <- F0
  P[1] <- P0
  L[1] <- L0
  I_A[1] <- 0
  I_F[1] <- 0
  I_P[1] <- 0
  B_A[1] <- BA0
  B_P[1] <- BP0
  B_F[1] <- BF0
  S_A[1] <- BA0
  S_P[1] <- BP0
  S_F[1] <- BF0
  D[1] <- 0
  C[1] <- C0
  
  # parameter sets
  ctrl_parms <- disc_parms %>%
    filter(Treatment == "control")
  
  fung_parms <- disc_parms %>%
    filter(Treatment == "fungicide")
  
  # parameters unaffected by disease
  grow_days <- seq(0, gs_time)

  s_A <- as.numeric(disc_parms[disc_parms$Parameter == "s_A", "Estimate"])
  b_A <- as.numeric(disc_parms[disc_parms$Parameter == "b_A", "Estimate"])
  gamma_A <- as.numeric(disc_parms[disc_parms$Parameter == "gamma_A", "Estimate"])
  
  s_P <- as.numeric(disc_parms[disc_parms$Parameter == "s_P", "Estimate"])
  b_F <- as.numeric(disc_parms[disc_parms$Parameter == "b_F", "Estimate"])
  gamma_P <- as.numeric(disc_parms[disc_parms$Parameter == "gamma_P", "Estimate"])
  
  d <- as.numeric(disc_parms[disc_parms$Parameter == "d", "Estimate"])
  h <- as.numeric(disc_parms[disc_parms$Parameter == "h", "Estimate"])
  
  dis_thresh <- as.numeric(disc_parms[disc_parms$Parameter == "dis_thresh", "Estimate"])
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # community disease severity
    comm_bio <- B_A[t] + B_F[t] + B_P[t]
    
    if(!is.na(comm_bio) & comm_bio > 0) {
      
      comm_dis <- (I_A[t] + I_F[t] + I_P[t])/comm_bio
    
    } else {
    
      comm_dis <- 0    
      
    }
    
    # change parameters with disease
    if(comm_dis > dis_thresh) {

      g_A <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "g_A", "Estimate"])
      c_A <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "c_A", "Estimate"])
      e_A <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "e_A", "Estimate"])

      g_P <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "g_P", "Estimate"])
      c_F <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "c_F", "Estimate"])
      c_P <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "c_P", "Estimate"])
      e_P <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "e_P", "Estimate"])
      l_P <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "l_P", "Estimate"])
      n_P <- as.numeric(ctrl_parms[ctrl_parms$Parameter == "n_P", "Estimate"])
      
      cont_parms_trt <- cont_parms %>%
        filter(Treatment == "control" | is.na(Treatment))
  
    } else {
      
      g_A <- as.numeric(fung_parms[fung_parms$Parameter == "g_A", "Estimate"])
      c_A <- as.numeric(fung_parms[fung_parms$Parameter == "c_A", "Estimate"])
      e_A <- as.numeric(fung_parms[fung_parms$Parameter == "e_A", "Estimate"])
      
      g_P <- as.numeric(fung_parms[fung_parms$Parameter == "g_P", "Estimate"])
      c_F <- as.numeric(fung_parms[fung_parms$Parameter == "c_F", "Estimate"])
      c_P <- as.numeric(fung_parms[fung_parms$Parameter == "c_P", "Estimate"])
      e_P <- as.numeric(fung_parms[fung_parms$Parameter == "e_P", "Estimate"])
      l_P <- as.numeric(fung_parms[fung_parms$Parameter == "l_P", "Estimate"])
      n_P <- as.numeric(fung_parms[fung_parms$Parameter == "n_P", "Estimate"])
      
      cont_parms_trt <- cont_parms %>%
        filter(Treatment == "fungicide" | is.na(Treatment))
      
    }
    
    # seedling establishment
    E_A <- e_A/(1 + gamma_A * L[t])
    E_P <- e_P/(1 + gamma_P * L[t])
    
    # new adult growth establishment
    N_P <- n_P/(1 + gamma_P * L[t])
    
    # initial conditions
    B_A_in <- log(g_A * E_A * b_A * A[t])
    B_F_in <- log(g_P * E_P * b_F * F1[t])
    B_P_in <- log(N_P * (B_P[t] + B_F[t]))
    C_in <- C[t] + h * (I_A[t] + I_F[t] + I_P[t])
    
    # print warning if population becomes NA
    if(is.na(A[t]) | is.na(F1[t]) | is.na(B_P[t]) | is.na(B_F[t])){
      print(cat("A:", A[t], "F:", F1[t], "BP:", B_P[t], "BF:", B_F[t], "time:", t))
    }
    
    # choose model function and starting conditions
    if(A[t] > 0 & F1[t] > 0 & (B_P[t] + B_F[t]) > 0) {
      
      cont_fun <- cont_AFP_mod
      ystart <- c(LogB_A = B_A_in, LogB_F = B_F_in, LogB_P = B_P_in, 
                  I_A = 0, I_F = 0, I_P = 0, D = 0, C = C_in)
      
    } else if(A[t] > 0 & F1[t] > 0 & (B_P[t] + B_F[t]) <= 0) {
      
      cont_fun <- cont_AF_mod
      ystart <- c(LogB_A = B_A_in, LogB_F = B_F_in, I_A = 0, I_F = 0, D = 0, C = C_in)
    
    } else if(A[t] > 0 & F1[t] <= 0 & (B_P[t] + B_F[t]) > 0) {
      
      cont_fun <- cont_AP_mod
      ystart <- c(LogB_A = B_A_in, LogB_P = B_P_in, I_A = 0, I_P = 0, D = 0, C = C_in)
    
    } else if(A[t] <= 0 & F1[t] > 0 & (B_P[t] + B_F[t]) > 0) {
      
      cont_fun <- cont_FP_mod
      ystart <- c(LogB_F = B_F_in, LogB_P = B_P_in, I_F = 0, I_P = 0, D = 0, C = C_in)
    
    } else if(A[t] > 0 & F1[t] <= 0 & (B_P[t] + B_F[t]) <= 0) {
      
      cont_fun <- cont_A_mod
      ystart <- c(LogB_A = B_A_in, I_A = 0, D = 0, C = C_in)
    
    } else if(A[t] <= 0 & F1[t] > 0 & (B_P[t] + B_F[t]) <= 0) {
      
      cont_fun <- cont_F_mod
      ystart <- c(LogB_F = B_F_in, I_F = 0, D = 0, C = C_in)
    
    } else if(A[t] <= 0 & F1[t] <= 0 & (B_P[t] + B_F[t]) > 0) {
      
      cont_fun <- cont_P_mod
      ystart <- c(LogB_P = B_P_in, I_P = 0, D = 0, C = C_in)
    
    } else {
      
      print("both species gone")
      
    }
    
    # biomass at the end of the growing season
    out <- ode(func = cont_fun,
               y = ystart,
               times = grow_days,
               parms = cont_parms_trt) %>%
      as_tibble() %>%
      filter(time == gs_time)
    
    # convert log values
    if("LogB_A" %in% colnames(out)){
      out$B_A <- exp(out$LogB_A)}
    if("LogB_F" %in% colnames(out)){
      out$B_F <- exp(out$LogB_F)}
    if("LogB_P" %in% colnames(out)){
      out$B_P <- exp(out$LogB_P)}
    
    # add missing values
    if(!("LogB_A" %in% colnames(out))){
      out$LogB_A <- 0
      out$B_A <- 0
      out$I_A <- 0 }
    if(!("LogB_F" %in% colnames(out))){
      out$LogB_F <- 0
      out$B_F <- 0
      out$I_F <- 0 }
    if(!("LogB_P" %in% colnames(out))){
      out$LogB_P <- 0
      out$B_P <- 0
      out$I_P <- 0 }
    
    # save biomass
    B_A[t+1] <- as.numeric(out$B_A)
    B_F[t+1] <- as.numeric(out$B_F)
    B_P[t+1] <- as.numeric(out$B_P)
    
    I_A[t+1] <- as.numeric(out$I_A)
    I_F[t+1] <- as.numeric(out$I_F)
    I_P[t+1] <- as.numeric(out$I_P)
    
    S_A[t+1] <- B_A[t] - I_A[t]
    S_F[t+1] <- B_F[t] - I_F[t]
    S_P[t+1] <- B_P[t] - I_P[t]
    
    D[t+1] <- as.numeric(out$D)
    C[t+1] <- as.numeric(out$C)
    
    # correct to prevent negative/very small numbers
    B_A[t+1] = ifelse(B_A[t+1] < 1e-10, 0, B_A[t+1])
    B_F[t+1] = ifelse(B_F[t+1] < 1e-10, 0, B_F[t+1])
    B_P[t+1] = ifelse(B_P[t+1] < 1e-10, 0, B_P[t+1])
    
    I_A[t+1] = ifelse(I_A[t+1] < 1e-10, 0, I_A[t+1])
    I_F[t+1] = ifelse(I_F[t+1] < 1e-10, 0, I_F[t+1])
    I_P[t+1] = ifelse(I_P[t+1] < 1e-10, 0, I_P[t+1])
    
    S_A[t+1] = ifelse(S_A[t+1] < 1e-10, 0, S_A[t+1])
    S_F[t+1] = ifelse(S_F[t+1] < 1e-10, 0, S_F[t+1])
    S_P[t+1] = ifelse(S_P[t+1] < 1e-10, 0, S_P[t+1])
    
    D[t+1] = ifelse(D[t+1] < 1e-10, 0, D[t+1])
    C[t+1] = ifelse(C[t+1] < 1e-10, 0, C[t+1])
    
    # population size
    A[t+1] = s_A * (1 - g_A) * A[t] +  c_A * B_A[t+1]
    L[t+1] = B_A[t+1] + B_F[t+1] + B_P[t+1] + D[t+1] + (1 - d) * L[t] 
    F1[t+1] = s_P * (1 - g_P) * F1[t] + c_F * B_F[t+1] + c_P * B_P[t+1]
    P[t+1] = l_P * P[t] + g_P * E_P * l_P * F1[t]
    
    # correct to prevent negative/very small numbers
    A[t+1] = ifelse(A[t+1] < 1e-10, 0, A[t+1])
    L[t+1] = ifelse(L[t+1] < 1e-10, 0, L[t+1])
    F1[t+1] = ifelse(F1[t+1] < 1e-10, 0, F1[t+1])
    P[t+1] = ifelse(P[t+1] < 1e-10, 0, P[t+1])
    
  } 
  
  # population data
  dfN <- tibble(time = rep(1:simtime, 14),
                species = rep(c("Annual",
                                "Perennial first-year", 
                                "Perennial adult", 
                                "Litter",
                                "Annual biomass",
                                "Perennial first-year biomass",
                                "Perennial adult biomass",
                                "Annual susceptible biomass",
                                "Perennial first-year susceptible biomass",
                                "Perennial adult susceptible biomass",
                                "Annual infected biomass",
                                "Perennial first-year infected biomass",
                                "Perennial adult infected biomass",
                                "Conidia"), 
                              each = simtime),
                N = c(A, F1, P, L, B_A, B_F, B_P, S_A, S_F, S_P, I_A, I_F, I_P, C))
  
  # return
  return(dfN)
}
cont_P_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_P <- x[1]
  I_P <- x[2]
  D <- x[3]
  C <- x[4]
  
  # parameters
  r_P <- as.numeric(params[params$Parameter == "r_P", "Estimate"])
  alpha_PP <- as.numeric(params[params$Parameter == "alpha_PP", "Estimate"])
  beta_PC <- as.numeric(params[params$Parameter == "beta_PC", "Estimate"])
  beta_PP <- as.numeric(params[params$Parameter == "beta_PP", "Estimate"])
  v_P <- as.numeric(params[params$Parameter == "v_P", "Estimate"])
  m_P <- as.numeric(params[params$Parameter == "m_P", "Estimate"])
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_P <- exp(LogB_P)
  S_P <- B_P - I_P
  
  # model with asymptotic transmission
  dLogBPdt <- r_P * (1 - alpha_PP * B_P) - m_P - v_P * I_P / B_P
  dIPdt <- beta_PC * S_P * C + beta_PP * S_P * I_P - (m_P + v_P) * I_P
  dDdt <-  m_P * B_P + v_P * I_P
  dCdt <- h * (m_P + v_P) * I_P - a * C
  
  # combine values
  dxdt <- c(dLogBPdt, dIPdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
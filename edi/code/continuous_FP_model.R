cont_FP_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_F <- x[1]
  LogB_P <- x[2]
  I_F <- x[3]
  I_P <- x[4]
  D <- x[5]
  C <- x[6]
  
  # parameters
  r_F <- as.numeric(params[params$Parameter == "r_F", "Estimate"])
  r_P <- as.numeric(params[params$Parameter == "r_P", "Estimate"])
  
  alpha_FF <- as.numeric(params[params$Parameter == "alpha_FF", "Estimate"])
  alpha_FP <- as.numeric(params[params$Parameter == "alpha_FP", "Estimate"])
  alpha_PF <- as.numeric(params[params$Parameter == "alpha_PF", "Estimate"])
  alpha_PP <- as.numeric(params[params$Parameter == "alpha_PP", "Estimate"])
  
  beta_FC <- as.numeric(params[params$Parameter == "beta_FC", "Estimate"])
  beta_PC <- as.numeric(params[params$Parameter == "beta_PC", "Estimate"])
  beta_FF <- as.numeric(params[params$Parameter == "beta_FF", "Estimate"])
  beta_FP <- as.numeric(params[params$Parameter == "beta_FP", "Estimate"])
  beta_PF <- as.numeric(params[params$Parameter == "beta_PF", "Estimate"])
  beta_PP <- as.numeric(params[params$Parameter == "beta_PP", "Estimate"])
  
  v_F <- as.numeric(params[params$Parameter == "v_F", "Estimate"])
  v_P <- as.numeric(params[params$Parameter == "v_P", "Estimate"])
  
  m_F <- as.numeric(params[params$Parameter == "m_F", "Estimate"])
  m_P <- as.numeric(params[params$Parameter == "m_P", "Estimate"])
  
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_F <- exp(LogB_F)
  B_P <- exp(LogB_P)
  S_F <- B_F - I_F
  S_P <- B_P - I_P
  
  # model with asymptotic transmission
  dLogBFdt <- r_F * (1 - alpha_FF * B_F - alpha_FP * B_P) - m_F - v_F * I_F / B_F
  dLogBPdt <- r_P * (1 - alpha_PF * B_F - alpha_PP * B_P) - m_P - v_P * I_P / B_P
  dIFdt <- beta_FC * S_F * C + beta_FF * S_F * I_F + beta_FP * S_F * I_P - (m_F + v_F) * I_F
  dIPdt <- beta_PC * S_P * C + beta_PF * S_P * I_F + beta_PP * S_P * I_P - (m_P + v_P) * I_P
  dDdt <- m_F * B_F + v_F * I_F + m_P * B_P + v_P * I_P
  dCdt <- h * ((m_F + v_F) * I_F + (m_P + v_P) * I_P) - a * C
  
  # combine values
  dxdt <- c(dLogBFdt, dLogBPdt, dIFdt, dIPdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
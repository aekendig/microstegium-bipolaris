cont_AFP_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_A <- x[1]
  LogB_F <- x[2]
  LogB_P <- x[3]
  I_A <- x[4]
  I_F <- x[5]
  I_P <- x[6]
  D <- x[7]
  C <- x[8]
  
  # parameters
  r_A <- as.numeric(params[params$Parameter == "r_A", "Estimate"])
  r_F <- as.numeric(params[params$Parameter == "r_F", "Estimate"])
  r_P <- as.numeric(params[params$Parameter == "r_P", "Estimate"])
  
  alpha_AA <- as.numeric(params[params$Parameter == "alpha_AA", "Estimate"])
  alpha_AF <- as.numeric(params[params$Parameter == "alpha_AF", "Estimate"])
  alpha_AP <- as.numeric(params[params$Parameter == "alpha_AP", "Estimate"])
  alpha_FA <- as.numeric(params[params$Parameter == "alpha_FA", "Estimate"])
  alpha_FF <- as.numeric(params[params$Parameter == "alpha_FF", "Estimate"])
  alpha_FP <- as.numeric(params[params$Parameter == "alpha_FP", "Estimate"])
  alpha_PA <- as.numeric(params[params$Parameter == "alpha_PA", "Estimate"])
  alpha_PF <- as.numeric(params[params$Parameter == "alpha_PF", "Estimate"])
  alpha_PP <- as.numeric(params[params$Parameter == "alpha_PP", "Estimate"])
  
  beta_AC <- as.numeric(params[params$Parameter == "beta_AC", "Estimate"])
  beta_FC <- as.numeric(params[params$Parameter == "beta_FC", "Estimate"])
  beta_PC <- as.numeric(params[params$Parameter == "beta_PC", "Estimate"])
  beta_AA <- as.numeric(params[params$Parameter == "beta_AA", "Estimate"])
  beta_AF <- as.numeric(params[params$Parameter == "beta_AF", "Estimate"])
  beta_AP <- as.numeric(params[params$Parameter == "beta_AP", "Estimate"])
  beta_FA <- as.numeric(params[params$Parameter == "beta_FA", "Estimate"])
  beta_FF <- as.numeric(params[params$Parameter == "beta_FF", "Estimate"])
  beta_FP <- as.numeric(params[params$Parameter == "beta_FP", "Estimate"])
  beta_PA <- as.numeric(params[params$Parameter == "beta_PA", "Estimate"])
  beta_PF <- as.numeric(params[params$Parameter == "beta_PF", "Estimate"])
  beta_PP <- as.numeric(params[params$Parameter == "beta_PP", "Estimate"])
  
  v_A <- as.numeric(params[params$Parameter == "v_A", "Estimate"])
  v_F <- as.numeric(params[params$Parameter == "v_F", "Estimate"])
  v_P <- as.numeric(params[params$Parameter == "v_P", "Estimate"])
  
  m_A <- as.numeric(params[params$Parameter == "m_A", "Estimate"])
  m_F <- as.numeric(params[params$Parameter == "m_F", "Estimate"])
  m_P <- as.numeric(params[params$Parameter == "m_P", "Estimate"])
  
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_A <- exp(LogB_A)
  B_F <- exp(LogB_F)
  B_P <- exp(LogB_P)
  S_A <- B_A - I_A
  S_F <- B_F - I_F
  S_P <- B_P - I_P
  
  # model with asymptotic transmission
  dLogBAdt <- r_A * (1 - alpha_AA * B_A - alpha_AF * B_F - alpha_AP * B_P) - m_A - v_A * I_A / B_A
  dLogBFdt <- r_F * (1 - alpha_FA * B_A - alpha_FF * B_F - alpha_FP * B_P) - m_F - v_F * I_F / B_F
  dLogBPdt <- r_P * (1 - alpha_PA * B_A - alpha_PF * B_F - alpha_PP * B_P) - m_P - v_P * I_P / B_P
  dIAdt <- beta_AC * S_A * C + beta_AA * S_A * I_A + beta_AF * S_A * I_F + beta_AP * S_A * I_P - (m_A + v_A) * I_A
  dIFdt <- beta_FC * S_F * C + beta_FA * S_F * I_A + beta_FF * S_F * I_F + beta_FP * S_F * I_P - (m_F + v_F) * I_F
  dIPdt <- beta_PC * S_P * C + beta_PA * S_P * I_A + beta_PF * S_P * I_F + beta_PP * S_P * I_P - (m_P + v_P) * I_P
  dDdt <- m_A * B_A + v_A * I_A + m_F * B_F + v_F * I_F + m_P * B_P + v_P * I_P
  dCdt <- h * ((m_A + v_A) * I_A + (m_F + v_F) * I_F + (m_P + v_P) * I_P) - a * C
  
  # combine values
  dxdt <- c(dLogBAdt, dLogBFdt, dLogBPdt, dIAdt, dIFdt, dIPdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
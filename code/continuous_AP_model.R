cont_AP_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_A <- x[1]
  LogB_P <- x[2]
  I_A <- x[3]
  I_P <- x[4]
  D <- x[5]
  C <- x[6]
  
  # parameters
  r_A <- as.numeric(params[params$Parameter == "r_A", "Estimate"])
  r_P <- as.numeric(params[params$Parameter == "r_P", "Estimate"])
  
  alpha_AA <- as.numeric(params[params$Parameter == "alpha_AA", "Estimate"])
  alpha_AP <- as.numeric(params[params$Parameter == "alpha_AP", "Estimate"])
  alpha_PA <- as.numeric(params[params$Parameter == "alpha_PA", "Estimate"])
  alpha_PP <- as.numeric(params[params$Parameter == "alpha_PP", "Estimate"])
  
  beta_AC <- as.numeric(params[params$Parameter == "beta_AC", "Estimate"])
  beta_PC <- as.numeric(params[params$Parameter == "beta_PC", "Estimate"])
  beta_AA <- as.numeric(params[params$Parameter == "beta_AA", "Estimate"])
  beta_AP <- as.numeric(params[params$Parameter == "beta_AP", "Estimate"])
  beta_PA <- as.numeric(params[params$Parameter == "beta_PA", "Estimate"])
  beta_PP <- as.numeric(params[params$Parameter == "beta_PP", "Estimate"])
  
  v_A <- as.numeric(params[params$Parameter == "v_A", "Estimate"])
  v_P <- as.numeric(params[params$Parameter == "v_P", "Estimate"])
  
  m_A <- as.numeric(params[params$Parameter == "m_A", "Estimate"])
  m_P <- as.numeric(params[params$Parameter == "m_P", "Estimate"])
  
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_A <- exp(LogB_A)
  B_P <- exp(LogB_P)
  S_A <- B_A - I_A
  S_P <- B_P - I_P
  
  # model with asymptotic transmission
  dLogBAdt <- r_A * (1 - alpha_AA * B_A - alpha_AP * B_P) - m_A - v_A * I_A / B_A
  dLogBPdt <- r_P * (1 - alpha_PA * B_A - alpha_PP * B_P) - m_P - v_P * I_P / B_P
  dIAdt <- beta_AC * S_A * C + beta_AA * S_A * I_A + beta_AP * S_A * I_P - (m_A + v_A) * I_A
  dIPdt <- beta_PC * S_P * C + beta_PA * S_P * I_A + beta_PP * S_P * I_P - (m_P + v_P) * I_P
  dDdt <- m_A * B_A + v_A * I_A + m_P * B_P + v_P * I_P
  dCdt <- h * ((m_A + v_A) * I_A + (m_P + v_P) * I_P) - a * C
  
  # combine values
  dxdt <- c(dLogBAdt, dLogBPdt, dIAdt, dIPdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
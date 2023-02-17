cont_AF_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_A <- x[1]
  LogB_F <- x[2]
  I_A <- x[3]
  I_F <- x[4]
  D <- x[5]
  C <- x[6]
  
  # parameters
  r_A <- as.numeric(params[params$Parameter == "r_A", "Estimate"])
  r_F <- as.numeric(params[params$Parameter == "r_F", "Estimate"])
  
  alpha_AA <- as.numeric(params[params$Parameter == "alpha_AA", "Estimate"])
  alpha_AF <- as.numeric(params[params$Parameter == "alpha_AF", "Estimate"])
  alpha_FA <- as.numeric(params[params$Parameter == "alpha_FA", "Estimate"])
  alpha_FF <- as.numeric(params[params$Parameter == "alpha_FF", "Estimate"])
  
  beta_AC <- as.numeric(params[params$Parameter == "beta_AC", "Estimate"])
  beta_FC <- as.numeric(params[params$Parameter == "beta_FC", "Estimate"])
  beta_AA <- as.numeric(params[params$Parameter == "beta_AA", "Estimate"])
  beta_AF <- as.numeric(params[params$Parameter == "beta_AF", "Estimate"])
  beta_FA <- as.numeric(params[params$Parameter == "beta_FA", "Estimate"])
  beta_FF <- as.numeric(params[params$Parameter == "beta_FF", "Estimate"])
  
  v_A <- as.numeric(params[params$Parameter == "v_A", "Estimate"])
  v_F <- as.numeric(params[params$Parameter == "v_F", "Estimate"])
  
  m_A <- as.numeric(params[params$Parameter == "m_A", "Estimate"])
  m_F <- as.numeric(params[params$Parameter == "m_F", "Estimate"])
  
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_A <- exp(LogB_A)
  B_F <- exp(LogB_F)
  S_A <- B_A - I_A
  S_F <- B_F - I_F
  
  # model with asymptotic transmission
  dLogBAdt <- r_A * (1 - alpha_AA * B_A - alpha_AF * B_F) - m_A - v_A * I_A / B_A
  dLogBFdt <- r_F * (1 - alpha_FA * B_A - alpha_FF * B_F) - m_F - v_F * I_F / B_F
  dIAdt <- beta_AC * S_A * C + beta_AA * S_A * I_A + beta_AF * S_A * I_F - (m_A + v_A) * I_A
  dIFdt <- beta_FC * S_F * C + beta_FA * S_F * I_A + beta_FF * S_F * I_F - (m_F + v_F) * I_F
  dDdt <- m_A * B_A + v_A * I_A + m_F * B_F + v_F * I_F
  dCdt <- h * ((m_A + v_A) * I_A + (m_F + v_F) * I_F) - a * C
  
  # combine values
  dxdt <- c(dLogBAdt, dLogBFdt, dIAdt, dIFdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
cont_F_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_F <- x[1]
  I_F <- x[2]
  D <- x[3]
  C <- x[4]
  
  # parameters
  r_F <- as.numeric(params[params$Parameter == "r_F", "Estimate"])
  alpha_FF <- as.numeric(params[params$Parameter == "alpha_FF", "Estimate"])
  beta_FC <- as.numeric(params[params$Parameter == "beta_FC", "Estimate"])
  beta_FF <- as.numeric(params[params$Parameter == "beta_FF", "Estimate"])
  v_F <- as.numeric(params[params$Parameter == "v_F", "Estimate"])
  m_F <- as.numeric(params[params$Parameter == "m_F", "Estimate"])
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_F <- exp(LogB_F)
  S_F <- B_F - I_F
  
  # model with asymptotic transmission
  dLogBFdt <- r_F * (1 - alpha_FF * B_F) - m_F - v_F * I_F / B_F
  dIFdt <- beta_FC * S_F * C + beta_FF * S_F * I_F - (m_F + v_F) * I_F
  dDdt <- m_F * B_F + v_F * I_F
  dCdt <- h * (m_F + v_F) * I_F - a * C
  
  # combine values
  dxdt <- c(dLogBFdt, dIFdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
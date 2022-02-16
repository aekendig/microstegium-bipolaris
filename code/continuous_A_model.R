cont_A_mod <- function(t, x, params) {
  
  # initial conditions
  LogB_A <- x[1]
  I_A <- x[2]
  D <- x[3]
  C <- x[4]
  
  # parameters
  r_A <- as.numeric(params[params$Parameter == "r_A", "Estimate"])
  alpha_AA <- as.numeric(params[params$Parameter == "alpha_AA", "Estimate"])
  beta_AC <- as.numeric(params[params$Parameter == "beta_AC", "Estimate"])
  beta_AA <- as.numeric(params[params$Parameter == "beta_AA", "Estimate"])
  v_A <- as.numeric(params[params$Parameter == "v_A", "Estimate"])
  m_A <- as.numeric(params[params$Parameter == "m_A", "Estimate"])
  h <- as.numeric(params[params$Parameter == "h", "Estimate"])
  a <- as.numeric(params[params$Parameter == "a", "Estimate"])
  
  # derived values
  B_A <- exp(LogB_A)
  S_A <- B_A - I_A
  
  # model with asymptotic transmission
  dLogBAdt <- r_A * (1 - alpha_AA * B_A) - m_A - v_A * I_A / B_A
  dIAdt <- beta_AC * S_A * C + beta_AA * S_A * I_A - (m_A + v_A) * I_A
  dDdt <- m_A * B_A + v_A * I_A
  dCdt <- h * (m_A + v_A) * I_A - a * C
  
  # combine values
  dxdt <- c(dLogBAdt, dIAdt, dDdt, dCdt)
  
  # output
  list(dxdt)
}
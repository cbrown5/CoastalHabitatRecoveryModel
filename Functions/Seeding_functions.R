
seeding_old_function <- function (x, phi = 400 , l_max = 1, nu = .01){
  # Deterministic Logistic repopulation probability based on x seeds
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift 
  # x_half: f(x) at half way through logistic curve 
  # ======================================
  l_max / (1 + exp(-nu*( x - phi)))  
}

seeding_check_par_function <- function (x, phi = 400, psi = 1 , l_max = 1, nu = .01){
  # Deterministic generalised logistic function (Richard's growth curve)
  # repopulation probability based on x seeds
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift 
  # x_half: f(x) at half way through logistic curve 
  # ======================================
  l_max / (1 + psi*exp(-nu*( x - phi)))^(1/psi)  
}

seeding_deter_function <- function (B, Q, d, dt, phi = 400 , l_max = 1, nu = .01){
  # Deterministic generalised logistic function (Richard's growth curve)
  # repopulation probability based on x seeds
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift - f(x) at half way through logistic curve 
  # ======================================
  l_max / (1 + exp(-nu*(t(Q) %*% B - phi))) * d * dt
}

seeding_stoch_function <- function (B, Q, d, dt, phi = 400, psi, l_max = 1, nu = .01){
  # Stochastic generalised logistic function (Richard's growth curve)
  # repopulation probability based on x seeds
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift - f(x) at half way through logistic curve  
  # ======================================
  
  p_rc <- l_max / (1 + psi*exp(-nu*(t(Q) %*% B - phi)))^(1/psi) # Calculate the probability of recolonisation over a year
  p_rc <-  1 - (1 - p_rc)^(1/(365*dt)) # Convert to probability of recolonisation over a time step
  # Probabilistically recolonise d amount if you are currently extinct
  (p_rc > runif(length(B))) * d * (B == 0)
  
}

seeding_simplified_function <- function (B, Q, d, dt, phi = 400, psi, l_max = 1, nu = .01){
  # Stochastic generalised logistic function (Richard's growth curve) SIMPLIFIED
  # repopulation probability based on x seeds
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift - f(x) at half way through logistic curve  
  # ======================================
  x <- t(Q) %*% B
  p_rc <- l_max / (1 + psi*exp(-nu*(t(Q) %*% B - phi)))^(1/psi) * (x > 0)# Calculate the probability of recolonisation over a year
  p_rc <-  1 - (1 - p_rc)^(1/(365*dt)) # Convert to probability of recolonisation over a time step
  # Probabilistically recolonise d amount if you are currently extinct
  (p_rc > runif(length(p_rc))) * d
  
}

rec_dependent_seeding_function <- function (B, Q, d, dt, phi = 400, psi, l_max = 1, nu = .01){
  # Stochastic generalised logistic function (Richard's growth curve) SIMPLIFIED
  # repopulation probability based on x seeds
  # This is to be used with for case 2 of the recruitment model
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift - f(x) at half way through logistic curve  
  # ======================================
  x <- t(Q) %*% B
  p_rc <- l_max / (1 + psi*exp(-nu*(x - phi)))^(1/psi) * (x > 0) # Calculate the probability of recolonisation over a year
  p_rc <-  1 - (1 - p_rc)^(1/(365*dt)) # Convert to probability of recolonisation over a time step
  
}

disp_dependent_seeding_function <- function (x, dt, phi = 400, psi, l_max = 1, nu = .01){
  # Stochastic generalised logistic function (Richard's growth curve) SIMPLIFIED
  # repopulation probability based on x seeds
  # This is to be used with for case 2 of the dispersal model
  # ======================================
  # lmax: maximum proportion of seeding
  # nu: the steepness parameter
  # phi: phase shift - f(x) at half way through logistic curve  
  # ======================================
  
  p_rc <- l_max / (1 + psi*exp(-nu*(x - phi)))^(1/psi) * (x > 0) # Calculate the probability of recolonisation over a year
  p_rc <-  1 - (1 - p_rc)^(1/(365*dt)) # Convert to probability of recolonisation over a time step
  p_rc
  }

check_seeding_fun <- function(B.max = 667, phi = 667/2, psi = 1, nu = 0.014, l_max = 1){  
  
  plot(y = seeding_check_par_function(0:B.max, phi = phi, psi = psi, nu = nu, l_max = l_max), x= 0:B.max,
       xlim = c(-100, B.max*1.5), ylim = c(0, 1),
       xlab = "Biomass_connected", ylab = "Probability of recolonisation (per yr)") # 667 is B.max
  
  summary(seeding_check_par_function(0:B.max, phi = phi, psi = psi, nu = nu))

}

0 ~ 0
Bmax ~ 1
# # Check that my probability interpretation is correct (which it is)
# ac <- 1 - 0.9241418
# rec <- 1 - (ac)^(1/(365*dt))  
# 
# sim_prob <- function(rec, t_periods){
#   x <- runif(t_periods)
#   any(x < rec)
# }
# 
# num <- 1000
# sum(replicate(num, expr = sim_prob(rec, 365*dt))) / num

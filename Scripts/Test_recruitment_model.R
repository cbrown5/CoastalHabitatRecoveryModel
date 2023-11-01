library(lme4)
library(tidyverse)

# source functions 
source("Functions/Conversion_functions.R")
source("Scripts/Stochastic_recruitment_model.R")
source("Functions/Seeding_functions.R")

# Get n
n <- 4

# Generate Q
Q <- matrix(data = rep(1, n^2), ncol = n, nrow = n)
diag(Q) <- 0
Q[1:2,] <-0 
Q[,1:2] <-0 
Q[2,1] <- 1
Q[1,2] <- 1


# Make the base parameter set
this_param_set <- list(
  
  # Patch parameters
  B.max = rep(667, n), 
  T. = rep(35,n),      # Temp                          degrees C
  I. = rep(1000, n),    # Irradiance                    mmol m-2 second-1
  M = rep(0.004,n), # Default 0.004 # NEED to have this specified
  B_init = c(0,0, rep(667, n-2)),
  
  # Inherent parameters
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),  # Max gross production (temp)   mg C g-1 dry wt h-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 319,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1), #,  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6,
  
  # Global parameters
  n = n,
  tmax = 365*100, # default 100 yrs
  dt = 1, # Step size
  Q = Q, # Seed transfer per biomass matrix
  E = 2, # Extinction threshold population - has to be less than 1% of B.max to repopulate
  
  # Seed transfer function parameters
  nu = 0.01, # steepness default 0.01 
  phi = 500, # phase shift of curve default 500
  psi = 1, # generalised logistic function parameter: logistic function = 1, Gompertz = 0+.
  d = 0.01 * 667, # Re-population biomass # needs d/dt > E otherwise can't repopulate
  l_max = 0.2 # Maximum probability of recovery in a year
)

rm(Q, n)

this_param_set[["B_init"]] <- Find_B_star(params = this_param_set, keep_perturbed = TRUE)

system.time({look <- Stochastic_recruitment_model(this_param_set, keep_timeseries = FALSE)})

attach(this_param_set)

look <- rec_dependent_seeding_function(B = c(0,0,0,0), Q, d, dt, phi, psi, l_max, nu)

detach(this_param_set)


# SEEDING FUNCTION WORKS
# Model never recovers when patches are not connected

# MC thinks everything is working fine




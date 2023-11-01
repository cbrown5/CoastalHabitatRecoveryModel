# This is the Halodule spp. parameter set 

halodule_param_set <- list(
  
  # Patch parameters
  B.max = 667, 
  M = 0.004, # Default 0.004 # NEED to have this specified
  
  # Inherent parameters
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),  # Max gross production (temp)   mg C g-1 dry wt h-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 223,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1), #,  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6,
  
  # Global parameters
  tmax = 365*100, # default 100 yrs
  dt = 1, # Step size
  E = 2, # Extinction threshold population - has to be less than 1% of B.max to repopulate
  
  # Seed transfer function parameters
  nu = 0.01, # steepness default 0.01 
  phi = 500, # phase shift of curve default 500
  psi = 1, # generalised logistic function parameter: logistic function = 1, Gompertz = 0+.
  d = 0.01 * 667, # Re-population biomass # needs d/dt > E otherwise can't repopulate
  l_max = 0.38 # Maximum probability of recovery in a year
  
)

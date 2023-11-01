# This is the zostera muelleri parameter set 

zostera_param_set <- list(
  
  # Patch parameters
  B.max = 558, 
  M = 0.017, # Default 0.004 # NEED to have this specified
  
  # Inherent parameters
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.3),  # Max gross production (temp)   mg C g-1 dry wt h-1
  T.opt = 30.9, # optimum temperature           degrees C
  T.max = 43.6, # Maximum temperature           degrees C
  Ik = 144,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(0.8), #,  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 36.7,  # Resp optimum temperature      degrees C
  RT.max = 44,
  
  # Global parameters
  tmax = 365*100, # default 100 yrs
  dt = 1, # Step size
  E = 2, # Extinction threshold population - has to be less than 1% of B.max to repopulate
  
  # Seed transfer function parameters
  nu = 0.01, # steepness default 0.01 
  phi = 500, # phase shift of curve default 500 #renamed tau in the ms
  psi = 1, # generalised logistic function parameter: logistic function = 1, Gompertz = 0+.
  d = 0.01 * 667, # Re-population biomass # needs d/dt > E otherwise can't repopulate
  l_max = 0.38 # Maximum probability of recovery in a year
  
)

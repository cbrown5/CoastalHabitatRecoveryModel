
Recovery_times_from_recolonisation <- 
  function(params, recovery_definition = c("B_equil","B_equil_half", "B.max_hlf")[[2]]){
  
  # ===============================
  # This Function computes the recovery times of a patches given their 
  # Stressor levels and Bmax
  #
  # Input: the parameter set of the experiment
  # Output: a vector of the same size as the perturbed patches
  # ===============================
  
  with(as.list(params), {
   
    # Create the recovery vector specified
    B_rec <- switch(recovery_definition,
           "B_equil" = Find_B_star(params),
           "B_equil_half" = Find_B_star(params)/2,
           "B.max_half" = B.max/2)
    
    # Compute fixed parameters
    ohm <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    g <- ohm - R - M
    alph <- -ohm/B.max
    F. <-  d/(alph*d + g) # where d is the starting biomass
    
    # Compute the recovery times
    Recovery_times <- log(B_rec/(F.*(alph*B_rec+g)))/g
    
    # Take only those from the perturbed patches
    return(Recovery_times[B_init == 0])
    
  })
}

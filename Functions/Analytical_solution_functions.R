
B_t <- function(params, t){
  # Compute the analytical solution at each time t based on the parameters
  # As per "Population_model_solution.docx" by Max Campbell
  
  with(params, {
    
    # Compute fixed parameters
    ohm <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    g <- ohm - R - M
    alph <- -ohm/B.max
    F. <-  B_init/(alph*B_init + g)
    
    # Compute B at time t
    return(g*F.*exp(g*t)/(1 - alph*F.*exp(g*t)))
  })
}


B_t_not_independent <- function(params, t, t_0){
  # Compute the analytical solution at each time t based on the parameters
  # As per "Population_model_solution.docx" by Max Campbell
  # This is for the case when patches are not independent
  # t_0 is the starting time of recolonisation
  
  with(params, {
    
    # Compute fixed parameters
    ohm <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    g <- ohm - R - M
    alph <- -ohm/B.max
    F. <-  d/(alph*d + g)
    
    
    # Compute B at time t
    return( g*F.*exp(g*(pmin(t-t_0, 1000)))/(1 - alph*F.*exp(g*(pmin(t-t_0, 1000)))) *  ( t >= t_0 )  )
  })
}


Find_B_star <- function(params, keep_perturbed = FALSE ){
  
  # =========================
  # Direct computation of the equilibrium biomass vector
  # =========================
  # keep_perturbed = TRUE you want to keep the perturbed patches at zero biomass 
  # =========================
  
  with(params, {
    
    # Compute fixed parameters
    ohm <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    
    output <- ifelse((ohm - R - M) <= 0, 0, (1-(R+M)/ohm)*B.max)
    
    # In case the none of the variables are vectors of the correct size
    if (length(output) != n){
      output<- rep(output, n)
    }
    
    # If specified replace the zeros
    if (keep_perturbed){
      output[B_init == 0] <- 0
    }
    
    return(output)
      
  })
}


source("Functions/resp_function.R")
source("Functions/photosynthesis_function.R")
source("Functions/Recovery_times_from_recolonisation.R")
source("Functions/Analytical_solution_functions.R")

Stochastic_dispersal_model <- function(params, keep_timeseries = TRUE){
  
  # ===============================
  # This is the implementation of the Stochastic dispersal model (one recolonisation probability
  # per perturbed edge connected to patch). The model was proposed by Chris Brown where each edge has an independent probability of recovering
  # the perturbed patch.
  #
  # This model can currently only handle degree 1 perturbations 
  #
  # Now if Xi was time at which the ith edge recovered the patch then we say that the system
  # is begins recovery at the time min(X_1, ..., X_N) ~ Exp(lambda_1 + ... + lambda_N)
  # ================================
  
params <- with(params, {

  # Check whether the system is in equilibrium
  if (any(Find_B_star(params = params, keep_perturbed = TRUE) != B_init)) {
    stop("ERROR ... Currently do not support source patches not at equilbrium for dispersal model :(")
  }
  
  B <- B_init
  
  # Case 1: The patches are independent and degree1 perturbation -------------------------------------
  
  if (sum(B_init==0) == 1){
  
    # Calculate the probability of recolonisation over a year based on generalised logistic function
    p_rc <- l_max / (1 + psi*exp(-nu*((B*Q[,B == 0]) - phi)))^(1/psi)
    #p_rc <- (l_max / (1 + psi*exp(-nu*(t(Q) %*% B - phi)))^(1/psi))[B == 0]
      
    lambda <- log(1-p_rc)/(-365) # Convert yearly probability into a daily rate for exponential rates
      
    recovery_rate <- ceiling(Q)[,B == 0] %*% lambda # Sum up the rates that are connected to perturbed patch
    
    # Calculate the recovery time of the patch in days
    recol_time <- Recovery_times_from_recolonisation(params)
    
    # Calculate quantities of interest
    expectation_R <- 1/recovery_rate + recol_time
    var_R <- 1/recovery_rate^2
    
    output <-  list(model="Edge model",
                    type = "Independent patches",
                    disturbance_type = "Degree 1",
                    
                    recovery_time = list( expectation = expectation_R, 
                                          variance = var_R,
                                          cdf = function (x, lambda, recol_time) {
                                            if (x > recol_time) {
                                            1 - exp(-lambda*(x - recol_time))
                                            } else { 0 }
                                            }, 
                                          lambda = recovery_rate,
                                          recol_time = recol_time,
                                          recolonisation_time = "Fixed analytical"),
                    parameters = params)
    
    return(output)

  # Case 2: The patches are NOT independent or many degree1 perturbations  -------------------------------------
  
  } else {
    
    which.perturbed <- which(B < 0.0000000001) # Which patches are perturbed
    
    
    # Find equilibrium biomass
    B_star <- Find_B_star(params = params)
    
    # time variables
    t_0 <- rep(NA, length(which.perturbed)) %>% setNames(which.perturbed)
    n_timesteps <-  tmax/dt
    times_vect <- seq(0, tmax, by = dt)
    
    
    # Initialise matrix to store results
    results <- cbind(B, matrix(B_star, n, n_timesteps))
    results[which.perturbed,-1] <- 0 # Choose zero here because simplifies computation
    
    # Simulate seeding probs to test
    seeding_prob_sims <- matrix(runif(n_timesteps*length(t_0)), nrow = length(t_0) )
    
    # Set perturbed patches at the beginning
    which.perturbed.still <- which.perturbed
    
    # Take only the connectivity we need 
    Q_short <- Q[,which.perturbed.still, drop = FALSE]
    
    
    # Compute seeding probabilities at t = 0
    seed_prob <- colSums(disp_dependent_seeding_function(x = (Q_short * B), dt = dt, 
                                                 phi = phi, psi = psi, l_max = l_max, nu = nu))
    
    
    # Find out when the next recolonisation event occurs
    when_realised <- unlist(lapply(
      lapply( split(seeding_prob_sims < as.numeric(seed_prob), seq(nrow(seeding_prob_sims))),FUN = which),
      FUN = min))
    
    # If no recolonisation is ever realised in the timeframe end experiment
    if (is.infinite(min(when_realised))){
      
      output <-  list(model="Stochastic dispersal",
                      type = "Dependent patches",
                      disturbance_type = "Degree 2 or higher",
                      recovery_time = list(recovery_time = NA, time_elasped = tmax),
                      parameters = params,
                      state_vars_time = NA)
      return(output) 
    }
    
    # Update t_0 with the first patch to recover
    which.realised <- when_realised == when_realised[which.min(when_realised)]
    t_0[as.character(which.perturbed.still[which.realised])] <- when_realised[which.realised]
    
    
    # Function to get the parameter set for a particular patch
    get_param <- function(x, patch_n){
      if (length(x) == 1){
        x
      } else {
        x[[patch_n]]
      }
    }
    
    
    # Add Biomass for patch/es that recolonise first
    for (j in seq_along(which.perturbed.still)){
      
      if (which.realised[[j]]) {
        patch_params <- lapply(params, get_param, patch_n = which.perturbed.still[[j]])
        # Compute the B for all time within the experiment
        results[(which.perturbed.still[[j]]),]<- B_t_not_independent(params = patch_params, 
                                                                     t = times_vect , t_0 = when_realised[[j]])
      }
    }
    
    # Now loop through updating the recolonisation probability, then finding the next recolonisation time
    # and the update the data, until we have all the recolonsation times
    while (any(is.na(t_0))) {
      
      # Update vars # Done first so we get the first recolonisation
      which.perturbed.still <- which.perturbed.still[!which.realised, drop = FALSE]
      seeding_prob_sims <- seeding_prob_sims[!which.realised, , drop = FALSE]
      Q_short <- Q_short[,!which.realised, drop = FALSE]
      
      # Compute new probabilities over time
      if (any((Q_short %% 1) != 0)){ # If graph is weighted
        
        # Find the new seeding probability for each timestep
        Q_list <- apply(Q_short, MARGIN = 2, FUN = function(x) list( "Q_vect"= x))
        p_rc <- t(sapply(Q_list, FUN = 
                          function(x) { colSums(disp_dependent_seeding_function(x = (results[,-1] * x[["Q_vect"]]) , dt = dt, 
                                                                        phi = phi, psi = psi, l_max = l_max, nu = nu)) }))
      
      } else { # Shortcut if graph is unweighted
      
        p_rc <- t(Q_short) %*% disp_dependent_seeding_function(x = results[,-1], dt = dt, 
                                        phi = phi, psi = psi, l_max = l_max, nu = nu)
      }
      
      # Figure out when the next recolonisation occurs
      when_realised <- unlist(lapply(
        lapply( split(seeding_prob_sims < p_rc, seq(nrow(seeding_prob_sims))),FUN = which),
        FUN = min))
      
      # Update the patch that is recolonised
      which.realised <- when_realised == when_realised[which.min(when_realised)]
      t_0[as.character(which.perturbed.still[which.realised])] <- when_realised[which.realised]
      
      # If no recolonisation is ever realised in the timeframe end experiment
      if (is.infinite(min(when_realised))){
        
        output <-  list(model="Stochastic dispersal",
                        type = "Dependent patches",
                        disturbance_type = "Degree 2 or higher",
                        recovery_time = list(recovery_time = NA, time_elasped = tmax),
                        parameters = params,
                        state_vars_time = NA)
        return(output) 
      }
      
      # Add Biomass for patch/es that recolonise first
      for (j in seq_along(which.perturbed.still)){
        
        if (which.realised[[j]]) {
          # Compute the B for all time within the experiment
          patch_params <- lapply(params, get_param, patch_n = which.perturbed.still[[j]])
          results[(which.perturbed.still[[j]]),]<- B_t_not_independent(params = patch_params, 
                                                                       t = times_vect , t_0 = when_realised[[j]])
        }
      }
      
    }
    
    # Create output -----------------------------------------------------------
    
    
    time_from_recol <- Recovery_times_from_recolonisation(patch_params)
    
    ## State variables over time data
    if (keep_timeseries){
      state_vars_time <- data.frame(state_var = rep(paste0(rep("B_", each = n), 
                                                           rep(1:n)), 
                                                    each = ncol(results)), 
                                    time = rep(seq(0, tmax, by = dt), nrow(results)), 
                                    value = as.numeric(matrix(t(results), ncol = 1)))
    } else {
      
      state_vars_time <- NA
      
    }
    
    output <-  list(model="Stochastic dispersal",
                    type = "Dependent patches",
                    disturbance_type = "Degree 2 or higher",
                    recovery_time = list(recovery_time = max(time_from_recol + t_0), time_elasped = tmax),
                    parameters = params,
                    state_vars_time = state_vars_time )
    return(output) 
    
  }
  })

}


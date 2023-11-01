
source("Functions/resp_function.R")
source("Functions/photosynthesis_function.R")
source("Functions/Recovery_times_from_recolonisation.R")
source("Functions/Analytical_solution_functions.R")

Stochastic_recruitment_model <- function(params, keep_timeseries = TRUE){ 
  
  # ==========================
  # This is the implementation of the Stochastic recruitment model (one recolonisation probability
  # per perturbed patch
  #
  # This model can handle all connectivity scenarios: 
  #
  # 1. Degree 1 (full analytical solution - returns variance and expectation),
  #
  # 2. Multiple degree 1 (analytical solution for recolonisation, recovery sometimes simulated 
  # - returns variance and expectation)
  #
  # 3. Dependent patches (no analytical available - returns a single replicate based on the parameters)
  #
  # keep_timeseries = FALSE option will speed up 3. by not returning the timeseries
  #
  # ==========================
  
  with(as.list(params), { 
    
    # Check whether the system is in equilibrium
    if (any(Find_B_star(params = params, keep_perturbed = TRUE) != B_init)) {
      stop("ERROR ... Currently do not support source patches not at equilbrium for recruitment model :(")
    }
    
    B <- B_init
    which.perturbed <- which(B < 0.0000000001) # Which patches are perturbed
    

# Case 1: The patches are independent -------------------------------------
    if(all(Q[which.perturbed, which.perturbed] == 0)){
      
      ## Find the recovery times based on stressors 
      time_from_recol <- Recovery_times_from_recolonisation(params)
      
      
      ## Find the distribution of time for last patch to recover (tlpr)
      
      # Calculate the probability of recolonisation for year
      p_rc <- l_max / (1 + psi*exp(-nu*(t(Q[,which.perturbed]) %*% B - phi)))^(1/psi)
      
      # input probability of recovery per year, output rate parameter in days
      find_rate_parameter <- function(p){ log(1-p)/(-365) } 
      
      ### When we only have 1 patch we can find everything directly from here
      if (length(p_rc) == 1){ # Recovery is exponentially distributed
        
        # Compute the rate of recovery per day
        lambda <- find_rate_parameter(p_rc)
        
        expectation_R <- 1/lambda + time_from_recol # Expected time till start of recolonisation + recolonisation time
        var_R <- 1/lambda^2 # Variance in the time till start of recolonisation 
        
        output <-  list(model="Stochastic recruitment",
                        type = "Independent patches",
                        disturbance_type = "Degree 1", 
                        recovery_time = list(expectation = expectation_R, variance = var_R,
                                             cdf = function (x, lambda, recol_time) {
                                               if (x > recol_time) {
                                                 1 - exp(-lambda*(x - recol_time))
                                               } else { 0 }
                                             }, 
                                             lambda = lambda,
                                             recol_time = time_from_recol,
                                              recolonisation_time = "Fixed analytical"),
                        parameters = params)
        
        return(output)
        
      }
      ###
      
    if (all(as.numeric(time_from_recol) == time_from_recol)) { # If the recolonisation times are the same
      
      # Assign the mean and variance (no randomness here)
      expectation_tlpr <- time_from_recol[1]
      var_tlpr <- 0
      
      recol_time_type <- "Fixed analytical"
      
      } else if (all(as.numeric(p_rc) == p_rc)) { # All probabilities are the same 
      
        # Directly compute the distribution of the recovery times
        distrib_final_patch <- rep(1/length(p_rc), length(p_rc)) 
        
        # Compute the expectation and variance estimation
        expectation_tlpr <- sum(distrib_final_patch * time_from_recol)
        var_tlpr <- sum((time_from_recol - expectation_tlpr)^2*distrib_final_patch)
        
        recol_time_type <- "Random analytical"
        
      } else { # Probabilities and recovery times differ - (here we use simulation)
      
      # Convert to per week for more efficiency (if we need speed can be gained here - tradeoff n sims and timesteps)
      p_rc_dt <-  1 - (1 - p_rc)^(1/(365/1)) 
      
      # Estimate the mean recovery time from recolonisation 
      # note this simulator will provide no result multiple final recolonisations start at the same time
      sims_vect <- unlist(replicate( { 
        recol_indicator <- rep(NA, length(p_rc_dt))
        
        # A sim (ends when we know which patch is last)
        while(sum(is.na(recol_indicator))>1){
          recol_indicator[runif(length(p_rc_dt)) < p_rc_dt] <- 1
        }
        which(is.na(recol_indicator))
        
        }, n = 100000) ) # How many sims?
      
      # Calculate the distribution
      sims_info <- table(sims_vect)/length(sims_vect)
      distrib_final_patch <- numeric(length(p_rc_dt))
      distrib_final_patch[as.numeric(names(sims_info))] <- as.numeric(sims_info)
      
      # Compute the expectation and variance estimation
      expectation_tlpr <- sum(distrib_final_patch * time_from_recol)
      var_tlpr <- sum((time_from_recol - expectation_tlpr)^2*distrib_final_patch)
      
      recol_time_type <- "Random simulated"
      
      }
      
      ##
      
      ### Calculate the expectation of recovery of all patches 
      # as per Info_analytical_computations_DE_system.RMD
      
      # Compute rates
      lambda_vec <- find_rate_parameter(p_rc)
      
      ## Compute the expectation of the time the final patch begins recolonisation (fsrt = "final starting recolonation time") 
      expectation_fsrt <- 0
      
      # Setup for the first sum
      sign_com <- 1
      m_temp <- 1
      
      while(m_temp <= length(lambda_vec)){
       
        comb <- combn(x = c(1:length(lambda_vec)), m = m_temp) # Find combinations
        expectation_fsrt <- expectation_fsrt + sign_com * sum(1 / colSums(matrix(as.numeric(lambda_vec)[comb], nrow = m_temp))) # Add sum to the expectation
        # update sign and sum
        sign_com <- sign_com * -1 
        m_temp <- m_temp + 1
      
      }
      ## 
      
      ## Compute the variance of the final starting recolonation time
      
      var_fsrt <- -(expectation_fsrt)^2
      
      # Setup for the first sum
      sign_com <- 1
      m_temp <- 1
      
      while(m_temp <= length(lambda_vec)){
        
        comb <- combn(x = c(1:length(lambda_vec)), m = m_temp) # Find combinations
        var_fsrt <- var_fsrt + 2 * sign_com * sum(1 / colSums(matrix(as.numeric(lambda_vec)[comb], nrow = m_temp))^2) # Add sum to the expectation
        # update sign and sum
        sign_com <- sign_com * -1 
        m_temp <- m_temp + 1
        
      }
      ##
      
      ### Compute the Variance and Expectation of the total recovery time of the system
      # This assumes independence which is not necessarily true (but should be close enough)
      expectation_R <- expectation_tlpr + expectation_fsrt
      var_R <- var_fsrt + var_tlpr
      
      # Create summary
      output <-  list(model="Stochastic recruitment",
                      type = "Independent patches",
                      disturbance_type = "Multiple Degree 1s",
                      recovery_time = list( expectation = expectation_R, variance = var_R, 
                                            recolonisation_time = recol_time_type),
                      parameters = params)
      
      return(output)
        
      
# Case 2: The patches are NOT independent -------------------------------------
      
      } else {
        
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
        
        # Start seeding probs via seeding function
        seed_prob <- rec_dependent_seeding_function( Q = Q_short, B = B, d = d , dt = dt, 
                                                   phi = phi, psi = psi, nu = nu, l_max = l_max)
        
        # Find out when the next recolonisation event occurs
        when_realised <- unlist(lapply(
          lapply( split(seeding_prob_sims < as.numeric(seed_prob), seq(nrow(seeding_prob_sims))),FUN = which),
          FUN = min))
        
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
        
        # If no recolonisation is ever realised in the timeframe end experiment
        if (is.infinite(min(when_realised))){
          
          output <-  list(model="Stochastic recruitment",
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
        
        # Find the new biomass for each timestep
        Biomass_sum <- t(Q_short) %*% results[,-1, drop = FALSE]
        
        # Compute the recolonisation probability for each timestep
        p_rc <- l_max / (1 + psi*exp(-nu*(Biomass_sum - phi)))^(1/psi) # Yearly
        p_rc <-  1 - (1 - p_rc)^(1/(365*dt)) # For each timestep 
        
        # Figure out when the next recolonisation occurs
        when_realised <- unlist(lapply(
          lapply( split(seeding_prob_sims < p_rc, seq(nrow(seeding_prob_sims))),FUN = which),
          FUN = min))
        
        # Update the patch that is recolonised
        which.realised <- when_realised == when_realised[which.min(when_realised)]
        t_0[as.character(which.perturbed.still[which.realised])] <- when_realised[which.realised]
        
        # If no recolonisation is ever realised in the timeframe end experiment
        if (is.infinite(min(when_realised))){
          
          output <-  list(model="Stochastic recruitment",
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
        
        output <-  list(model="Stochastic recruitment",
                        type = "Dependent patches",
                        disturbance_type = "Degree 2 or higher",
                        recovery_time = list(recovery_time = max(time_from_recol + t_0), time_elasped = tmax),
                        parameters = params,
                        state_vars_time = state_vars_time )
        return(output) 
      
      }
  })
}






#################################################################################################
#
# Sensitivity analysis - Biomass model with connectivity layer
#
# Zostera
#
# ====================
# Calculate expected recovery time and its quantiles for one perturbed patch where we 
# vary number of source patches between 1-20. Plot recovery time vs number of
# source patches for each 
# stressor combo and model (facets)
# 
# 1-20 source nodes, 1 perturbations
# 20 connectivity scenarios
# Full stressor combos 
# Analytical solutions available
#
################################################################################################

library(tidyverse)

# source functions 
source("Functions/Conversion_functions.R")
source("Functions/Stochastic_dispersal_model.R")
source("Functions/Stochastic_recruitment_model.R")
source("Functions/Seeding_functions.R")
source("Scripts/zostera_params.R")

# Specify experiment name HERE (and stress levels)
light_stress <- 2.4
Experiment_name <- gsub(".", "-", 
                        sprintf("Deg1_stochastic-rec-disp_L%s", light_stress),
                        fixed = TRUE)

######################## Set Parameters ##########################################


SD_model_list <- list()
SR_model_list <- list()
treatment_df <- list()

# Specify how many nodes (patches)
n_vect <- 2:10
 
# Setup the stressor treatments
# Choose light and temp levels plus name
lhour <- function(x) x*1000/24 # Function to convert mol per day to mmol per hour
light_lev <- lhour(c(1.3, 2.6, 5.2))#lhour(seq(1, 5.2, by = 0.1))
temp_lev <- c( "normal" = 30.3) 
templight <- expand.grid( T. = temp_lev, I. = light_lev)

# For each number of nodes
for (i in seq_along(n_vect)){
  
  # Get n
  n <- n_vect[[i]]
  
  # Generate Q
  Q <- matrix(data = numeric(n^2), ncol = n, nrow = n)
  Q[1,] <- 1
  Q[, 1] <- 1
  Q[1,1] <- 0

  this_param_set <- within(zostera_param_set, {
    
    n <- n
    B_init <-  c(0, rep(667, n-1))
    Q <- Q
  
    })
  
  rm(Q, n)
  
  # For each stressor treatment
  for (j in 1:nrow(templight)){
    
    temp_stressors <- templight[j,, drop = FALSE]
    
    # Update the parameters
    temp_params <- this_param_set
    temp_params[["I."]] <- as.vector(light_lev[j])
    temp_params[["T."]] <- temp_lev
    
    # Start source patches at equilibrium
    temp_params[["B_init"]] <- Find_B_star(params = temp_params, keep_perturbed = TRUE)

    SD_model_list <- c(SD_model_list, list(Stochastic_dispersal_model(temp_params)))
    SR_model_list <- cbind(SR_model_list, list(Stochastic_recruitment_model(temp_params)))
    
    treatment_df <- rbind(treatment_df, within(temp_stressors, n <-  n_vect[[i]]))
    
    
  }

}

# Get the variables we need
treatment_df <- within(treatment_df, {
  n_sources <- n-1
  # Stochastic Dispersal results
  disp_mod_lambda <- lapply(SD_model_list, FUN = 
                              function(x) x[["recovery_time"]][["lambda"]]) %>% unlist()
  disp_mod_recol_time <- lapply(SD_model_list, FUN = 
                                  function(x) x[["recovery_time"]][["recol_time"]]) %>% unlist()
  
  # Stochastic recruitment results
  rec_mod_lambda <- lapply(SR_model_list, FUN = 
                             function(x) x[["recovery_time"]][["lambda"]]) %>% unlist()
  rec_mod_recol_time <- lapply(SR_model_list, FUN = 
                                 function(x) x[["recovery_time"]][["recol_time"]]) %>% unlist()
  
  })

treatment_df <- treatment_df %>% pivot_longer(cols = disp_mod_lambda:rec_mod_recol_time,
                                              names_sep= "_mod_",
                                              names_to = c("model", "variable"),
                                              values_to = c("par_value")) %>% 
  pivot_wider(names_from = "variable", values_from = "par_value")  

# Both models have the same cdf function
mod_cdf <- SR_model_list[[1]][["recovery_time"]][["cdf"]]

# By letting the cdf = the percentile and rearranging to find x we get the following function
percentile_recovery <- function(percentile, lambda, recol_time) {
  -log(1-percentile/100)/lambda + recol_time 
  }


# Check cdf is working as expected # which it is
# v_mod_cdf <- Vectorize(mod_cdf, vectorize.args = "x")
# plot(y = v_mod_cdf(x = 0:30, lambda = treatment_df$lambda[[1]], recol_time = 16), x = 0:30, type = "l")

# Compute the percentiles
treatment_df <- within(treatment_df, {
  
  median_recovery <- percentile_recovery(50, lambda = lambda, recol_time = recol_time)
  per5 <- percentile_recovery(5, lambda = lambda, recol_time = recol_time)
  per95 <- percentile_recovery(95, lambda = lambda, recol_time = recol_time)
  per25 <- percentile_recovery(25, lambda = lambda, recol_time = recol_time)
  per75 <- percentile_recovery(75, lambda = lambda, recol_time = recol_time)
  model<- case_when(model == "rec" ~ "Stochastic recolonisation",
                    model == "disp" ~ "Stochastic dispersal")
  species <- "zostera"
  
  })

save(treatment_df, 
     file = "Outputs/zostera-stressors-scenarios.rda")


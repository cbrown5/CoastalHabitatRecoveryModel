
# Zostera analysis, mgmt simulations 
#CJ Brown 2024-02-16

library(tidyverse)


# source functions 
source("Functions/Conversion_functions.R")
source("Functions/Stochastic_dispersal_model.R")
source("Functions/Stochastic_recruitment_model.R")
source("Functions/Seeding_functions.R")
source("Scripts/zostera_params.R")

lhour <- function(x) x*1000/24 # Function to convert mol per day to mmol per hour
# Specify experiment name HERE (and stress levels)
temp_stress <- 34.5
light_stress <- seq(1, 5.1, by = 0.1)
  # c(2.6, 1.3)
######################## Set Parameters ##########################################

SD_model_list <- list()
SR_model_list <- list()
treatment_df <- list()

# Specify how many nodes (patches)
n_vect <- 2:11
 
# Setup the stressor treatments
# Choose light and temp levels plus name
temp_lev <- 30.3# c( "normal" = 30.3, "high" =  temp_stress) 
light_lev <-lhour(c(5.2, light_stress)) #c("normal" = lhour(5.2), "low" = lhour(light_stress))
templight <- expand.grid( T. = temp_lev, I. = light_lev)

#change parameters to saturate recruitment and/or dispersal
all_param_combos <- rbind(
  data.frame(templight, 
             phi = zostera_param_set$phi,
             l_max = zostera_param_set$l_max),
  data.frame(templight, 
             phi = zostera_param_set$phi*1.5,
             l_max = zostera_param_set$l_max),
  #This represents recruit/settlement saturation (ie restoration interventions, 
  # with type of intervention depending on which model we are looking at): 
  data.frame(templight, 
             #phi: biomass at which recruit probability is half Lmax (called tau in the paper)
            phi = 20,
             l_max = 0.999)
  )



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
  for (j in 1:nrow(all_param_combos)){
    
    this_combo <- all_param_combos[j,, drop = FALSE]
    
    # Update the parameters
    temp_params <- this_param_set
    temp_params[["phi"]] <- as.vector(this_combo$phi)
    temp_params[["l_max"]] <- as.vector(this_combo$l_max)
    temp_params[["T."]] <- as.vector(this_combo$T.)
    temp_params[["I."]] <- as.vector(this_combo$I.)
    
    # Start source patches at equilibrium
    temp_params[["B_init"]] <- Find_B_star(params = temp_params, 
                                           keep_perturbed = TRUE)

    SD_model_list <- c(SD_model_list, 
                       list(Stochastic_dispersal_model(temp_params)))
    SR_model_list <- cbind(SR_model_list, 
                           list(Stochastic_recruitment_model(temp_params)))
    
    treatment_df <- rbind(treatment_df,
                          within(this_combo, n <-  n_vect[[i]]))
    
    
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




save(treatment_df, all_param_combos, 
     file = "Outputs/zostera-mgmt-results.rda")


# ------------ 
# PLOTS 
# ------------ 


#
# Figure 3, just two models 
#

# ggplot(filter(treatment_df, I. == templight$I.[1])) + 
  ggplot(treatment_df) +
  aes(x = n_sources, y = median_recovery/365, colour = T.,
      shape = model, fill = I.) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = model), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1,2)) + 
  xlab("Number of source patches") + 
  ylab("Recovery time (years)") + 
  scale_x_continuous(breaks = 1:10) +
  facet_grid(phi ~ l_max, scales = "free") + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10)
  )


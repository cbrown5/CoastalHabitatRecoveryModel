
#################################################################################################
#
# Sensitivity analysis for STRESSNET MODEL - Biomass model with connectivity layer
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
temp_stress <- 34.5
light_stress <- 2.4
Experiment_name <- gsub(".", "-", 
                        sprintf("Deg1_stochastic-rec-disp_T%s_L%s", temp_stress, light_stress),
                        fixed = TRUE)

######################## Set Parameters ##########################################


SD_model_list <- list()
SR_model_list <- list()
treatment_df <- list()

# Specify how many nodes (patches)
n_vect <- 2:11
 
# Setup the stressor treatments
# Choose light and temp levels plus name
lhour <- function(x) x*1000/24 # Function to convert mol per day to mmol per hour
temp_lev <- c( "normal" = 30.3, "high" =  temp_stress) 
light_lev <- c("normal" = lhour(5.2), "low" = lhour(light_stress))
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
    temp_params[["T."]] <- as.vector(temp_stressors$T.)
    temp_params[["I."]] <- as.vector(temp_stressors$I.)
    
    # Start source patches at equilibrium
    temp_params[["B_init"]] <- Find_B_star(params = temp_params, keep_perturbed = TRUE)

    SD_model_list <- c(SD_model_list, list(Stochastic_dispersal_model(temp_params)))
    SR_model_list <- cbind(SR_model_list, list(Stochastic_recruitment_model(temp_params)))
    
    treatment_df <- rbind(treatment_df, within(temp_stressors, n <-  n_vect[[i]]))
    
    
  }

}

# Get the variables we need
treatment_df <- within(treatment_df, {
  light <- if_else( I.== light_lev[["normal"]], "normal", "low")
  temp <- if_else(T. == temp_lev[["normal"]], "normal", "high")
  Stressors <- case_when(light == "normal" & temp == "normal" ~ "Control", 
                        light == "low" & temp == "normal" ~ "Low light",
                        light == "normal" & temp == "high" ~ "High temp",
                        light == "low" & temp == "high" ~ "Low light + high temp")
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


# ------------ 
# PLOTS 
# ------------ 


#
# Figure 3, just two models 
#

mod2_res <- treatment_df %>% 
  filter(Stressors == "Control")

ggplot(mod2_res) + 
  aes(x = n_sources, y = median_recovery/365, colour = model,
      shape = model) +
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
  scale_color_manual(values = c("black", "#d41515")) + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = c(0.8,0.8),
        legend.title = element_blank()
  )

#
# Effect of stressors
#

ggplot(treatment_df) + aes(x = n_sources, y = median_recovery/365, colour = Stressors) +
  geom_hline(yintercept = 0) +
  geom_point(aes(shape = Stressors), position = position_dodge(width = 0.6)) + theme_classic()+
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(1,2,4,5)) + 
  facet_wrap(.~ model, nrow = 2, scales = "free_y") + xlab("Number of source patches") + 
  scale_x_continuous(breaks = 1:10) +
  ylab("Recovery time (years)") + 
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), legend.position = c(0.8,0.3),
          legend.title = element_blank()
        )
  
ggsave(sprintf("Shared/Figures/%s_recovery_time_vs_npatches.png", 
               Experiment_name), width = 18, height = 15, units = "cm", dpi = 300, 
       scale = 0.9)



#################################################################################################
#
# Sensitivity analysis for STRESSNET MODEL - Biomass model with connectivity layer
#
# ====================
# Experiment X:
# X nodes,
# X connectivity scenarios
# Full stressor combos 
# Ran every X for X years
#
################################################################################################

library(tidyverse)
library(data.table)
library(parallel)
library(SOAR)
library(igraph)

# source functions 
source("Functions/Conversion_functions.R")
source("Scripts/Stochastic_dispersal_model.R")
source("Scripts/Stochastic_recruitment_model.R")
source("Functions/Seeding_functions.R")
source("Functions/Generate_graphs_functions.R")

# Specify experiment name HERE
Experiment_name <- ""

######################## Set Parameters ##########################################

# Specify how many nodes (patches)
n <- 60

## Set seed transfer per biomass matrix

#Q <- matrix(0, n, n) # No transfer
#diag(Q) <- 1

Q <- matrix(1/(n-1), n, n) # Fully distributed
diag(Q) <- 0

##

## Help choosing seeding parameters # Try these ones
# png("Figures/recolonisation_functions.png", width = 2000, height = 500, pointsize = 20)
# par(mfrow = c(1,3)) # Set to however many functions we use
check_seeding_fun(B.max = 667, phi = 667/2, psi = 1, nu = 0.014)

# dev.off()

this_param_set <- list(
  
  # Patch parameters
  B.max = rep(667, n), 
  T. = rep(35, n),      # Temp                          degrees C
  I. = rep(1000, n),    # Irradiance                    mmol m-2 second-1
  M = rep(0.004, n), # Default 0.004
  B_init = c(rep(667, n-2), 0,0),
  
  # Inherent parameters
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),  # Max gross production (temp)   mg C g-1 dry wt h-1
  #PL.max = umol_to_mg(422), # Max gross production (light)  mmol O2 g dry wt-1 hr-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 319,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1), #,  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6,
  # theta = 24 * (1 / 0.33) / 1000, # Conversion factor for:   mg C h-1  ->  g dry wt day-1
  
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

# Start source patches at equilibrium
this_param_set[["B_init"]] <- Find_B_star(params = this_param_set, keep_perturbed = TRUE)

#look <- run_DE(this_param_set)

source("Scripts/Stochastic_recruitment_model_analytical_trial.R")
system.time({look3 <- replicate(n = 100,{
  Stochastic_recruitment_model(params = this_param_set, keep_timeseries = FALSE)
  })})

# Run connectivity experiment ---------------------------------------------

# Params experiment
reps <- 50 # default 50
set.seed(02022021) # default today's date
quantile_thinning <- 1 # default is dt - in days (to speed up make higher), only change resolution of the quantile plot
SAVE <- FALSE

# Specify connectivity scenarios

sample_models <- list()

# SPECIFY MODELS TO GENERATE GRAPHS FROM FOR EXPERIMENT
sample_models[[1]] <- expression(sample_pa(n = n, m = 2, power = 0, directed = FALSE))
sample_models[[2]] <- expression(sample_pa(n = n, m = 2, power = 1, directed = FALSE))
sample_models[[3]] <- expression(sample_pa(n = n, m = 2, power = 2, directed = FALSE))
sample_models[[4]] <- expression(sample_k_regular(no.of.nodes = n, k = 2, directed = FALSE, multiple = FALSE)) 
sample_models[[5]] <- expression(sample_smallworld(1, 60, nei = 1, 0.5))

sim_output <- simulate_graphs(sample_models, reps = 5)
graph_metric_info <- sim_output[[1]]
simulated_networks <- sim_output[[2]]

c1 <- list(nu = 0.032, 
           phi = 50,
           psi = 1e-12)

c2 <- list(nu = 0.014, 
           phi = 667/2,
           psi = 1)

c3 <- list(nu = 0.032, 
           phi = 600,
           psi = 1)

connect_sets <- list(c1, c2, c3)

# Choose light and temp levels plus name
temp_lev <- c( "normal" = 35, "high" =  42.58) 
light_lev <- c("normal" = 1000, "low" = 254.55)
templight <- expand.grid( T. = temp_lev, I. = light_lev)

# Initialise storage
param_lists <- NULL 
param_list_id <- NULL


# Generate parameter sets
for (i in 1:nrow(templight)){
  
  for (j in 1:length(connect_sets)){
    
    # Specify the parameter set for the treatment
    connect_set <- connect_sets[[j]]
    stress_connect_params <- within(this_param_set,{
      phi <-  connect_set[["phi"]]
      psi <-  connect_set[["psi"]]
      nu <-  connect_set[["nu"]]
      T. <- as.vector(templight[i,"T."])
      I. <-  as.vector(templight[i,"I."])
      
    })
    
    # Set B_init at equilibrium given the stressor levels
    stress_connect_params[["B_init"]] <- Find_Binit_burnin(stress_connect_params, abs.tol = 0.001, 
                                                           timestep = 0.5)
    
    # Duplicate parameter set for replicates
    param_lists <- c(param_lists, rep(list(stress_connect_params), reps))
    # Keep rep metadata
    param_list_id <- rbind(param_list_id, data.frame(connectivity = LETTERS[j], 
                                                     light = names(light_lev)[light_lev == templight[i,"I."]],
                                                     temp = names(temp_lev)[temp_lev == templight[i,"T."]]))
  }
}

# Run the model parallel based on the specified parameter sets  
DE_results <- mclapply(param_lists, FUN = run_DE, 
                       mc.cores = detectCores() - 1)

#state_vars_time <- NULL
DE_time_summary <- NULL
quantile_data <- NULL

# Combine data and add treatment and connectivity info
system.time({
  for (i in 1:(nrow(templight) * length(connect_sets))){ # For the number of treatment:connectivities
    # Add treatment and connectivity info to time series data
    system.time({ DE_treatment <- within(do.call("rbind", DE_results[((i-1) * reps+1):(i * reps)]), {
      connectivity <- param_list_id[i,]$connectivity
      temp <- param_list_id[i,]$temp
      light <- param_list_id[i,]$light
      # Add reps based on data structure: tmax/dt * n * number of vars per node + number of vars per node for initial conditions
      rep_no <- rep(1:reps, each = (this_param_set[["tmax"]] / this_param_set[["dt"]]*this_param_set[["n"]]*2 +
                                      this_param_set[["n"]]*2)) # Initial conditions
      
    })})
    
    # Compute the recovery times for each rep
    system.time({ DE_time_summary_temp <- DE_treatment %>% 
      filter(state_var %in% paste0("B_", which(param_lists[[reps*i]]$B_init == 0)), # only the empty patches
             ( value >= max(param_lists[[reps*i]]$B_init)/2 | # only the one that are over half equilibrium biomass
                 time == this_param_set[["tmax"]])) %>%  
      group_by(rep_no, temp, light, connectivity) %>% 
      summarise(Recovery_time = min(time)) %>% ungroup() %>% 
      mutate(B_max_stressors = max(param_lists[[reps*i]]$B_init)) })
    
    # Compute quantile data - this is around 50 times quicker than tidyverse way
    quantile_treat <-  data.table(DE_treatment)[grepl("B_", DE_treatment$state_var) & 
                                                  ((DE_treatment$time %% quantile_thinning) == 0), ]
    quantile_treat <- setDT(quantile_treat)[,.(med_value = median(value), quant25  = quantile(value, probs = 0.25), 
                                               quant75  = quantile(value, probs = 0.75)), 
                                            by = 'state_var,time,connectivity,temp,light']
    
    # Bind data
    # state_vars_time <- rbind(state_vars_time, DE_treatment) # This is too large to practically store
    DE_time_summary <- rbind(DE_time_summary, DE_time_summary_temp)
    quantile_data <- rbind(quantile_data, quantile_treat)
    
  }
  
  # Add treatment variable
  DE_time_summary <- within(DE_time_summary, {
    treatment <- case_when(light == "normal" & temp == "normal" ~ "control", 
                           light == "low" & temp == "normal" ~ "low light",
                           light == "normal" & temp == "high" ~ "high temp",
                           light == "low" & temp == "high" ~ "low light + high temp")
    Recovery_time[Recovery_time == this_param_set[["tmax"]]] <- NA
  })
  
  # Add treatment variable
  quantile_data <- within(quantile_data, {
    treatment <- case_when(light == "normal" & temp == "normal" ~ "control", 
                           light == "low" & temp == "normal" ~ "low light",
                           light == "normal" & temp == "high" ~ "high temp",
                           light == "low" & temp == "high" ~ "low light + high temp")
  })
})

# Store locally to free up RAM - takes a while
# Store(DE_results)
# Attach()

# Clear unwanted objects
rm(DE_treatment, quantile_treat, DE_time_summary_temp)

gc()


# Save OBJECTS if specified
if (SAVE) {
  # Generate filename  
  filename_details <- paste0(Experiment_name, "_",  Sys.Date(), "_", 
                             reps, "reps_", this_param_set[["dt"]], "daysteps_", this_param_set[["tmax"]]/365, "yrs", 
                             ifelse(quantile_thinning == this_param_set[["dt"]],"" ,paste0("_", quantile_thinning , "qthin")))
  save(DE_time_summary, this_param_set, quantile_data, param_lists, param_list_id, 
       file = paste0("Data/", filename_details, "_objects.RDA"))
}



# Recovery calculations ---------------------------------------------------

# Need to extend the number of timesteps fairly far way
av_recovery_times <- DE_time_summary %>%
  group_by(treatment, temp, light, connectivity) %>%
  summarise(mean_recovery_time = mean(c(Recovery_time, rep(this_param_set[["tmax"]], times = sum(is.na(Recovery_time)))), na.rm = TRUE), 
            median_recovery_time = median(c(Recovery_time, rep(this_param_set[["tmax"]], times = sum(is.na(Recovery_time)))), na.rm = TRUE), 
            median_recovery_time_yrs = median(c(Recovery_time, rep(this_param_set[["tmax"]], times = sum(is.na(Recovery_time)))), na.rm = TRUE)/365,
            quant25  = quantile(c(Recovery_time, rep(this_param_set[["tmax"]], times = sum(is.na(Recovery_time)))), probs = 0.25, na.rm = TRUE), 
            quant75  = quantile(c(Recovery_time, rep(this_param_set[["tmax"]], times = sum(is.na(Recovery_time)))), probs = 0.75, na.rm = TRUE),
            Prop_recovered = sum(!is.na(Recovery_time))/reps)

# # Proportion of simulations where B_2 recovered to B.max
# prop_year_recovered_Bmax <- state_vars_time %>% filter(time == this_param_set[["tmax"]], state_var == "B_2") %>% 
#   group_by(treatment, temp, light, connectivity) %>% mutate(recovered = (value >= (667/2))) %>% 
#   summarise( recovery_prop = mean(recovered))

# # Proportion of simulations where B_2 recovered to half of B_1
# prop_year_recovered <- state_vars_time %>% filter(time == this_param_set[["tmax"]]) %>%
#   pivot_wider(names_from = "state_var", values_from = "value") %>%
#   group_by(treatment, temp, light, connectivity) %>%
#   mutate(recovered = ((B_1/2) <= B_2)) %>%
#   summarise( recovery_prop = mean(recovered))

if (SAVE) { # Save RECOVERY SUMMARY if specified
  write_csv(av_recovery_times, path = paste0("Outputs/", filename_details,"_recovery_summary.csv"))
}

# Save some outputs
# write_csv(prop_year_recovered_Bmax, path = "Outputs/28-01-Experiment1_recovered1yr_halfBmax_100sims.csv")
# write_csv(prop_year_recovered, path = "Outputs/28-01-Experiment1_recovered1yr_halfB1_100sims.csv")

# Make plots --------------------------------------------------------------

# Look at state variables quantiles together over time 
p1 <- ggplot(quantile_data) + 
  aes(x = time/365, y = med_value, colour = state_var) + 
  geom_line() + 
  geom_line(aes(x = time/365, y = quant25), linetype="dashed") +
  geom_line(aes(x = time/365, y = quant75), linetype="dashed") + 
  theme_bw() + scale_x_continuous(breaks = seq(0, this_param_set[["tmax"]]/365, 10)) +
  facet_wrap(.~ treatment + connectivity, scales = "free", nrow = 4) + 
  geom_hline(yintercept = 0, linetype="dashed" ) + ylab("Biomass (g dry wt m-2)") +
  xlab("Time (Years)") 

p1 

if (SAVE) { # Save QUANTILES OVER TIME if specified
  ggsave(paste0("Figures/", filename_details,"_overtime_25-50-75quantile.png"), width = 10, height = 10)
}

# Change NAs to tmax so they show up in the histogram
DE_time_summary_na_rep <- DE_time_summary
DE_time_summary_na_rep$Recovery_time[is.na(DE_time_summary_na_rep$Recovery_time)] <- this_param_set[["tmax"]]

p2 <-  ggplot(DE_time_summary_na_rep) + 
  aes(x = Recovery_time/365) + 
  geom_histogram(stat = "bin", bins = 100) + theme_bw() +
  facet_wrap(.~ treatment + connectivity, nrow = 4) + 
  geom_hline(yintercept = 0) + ylab("Frequency of patches") +
  xlab("Recovery time (Years)") 

p2

if (SAVE) { # Save RECOVERY HISTOGRAM if specified
  ggsave(paste0("Figures/", filename_details,"_recovery_hists.png"), width = 10, height = 10)
}

# Cut off parts of the data where the empty patches sit at the stable equilibrium
quantile_data_zoomed <- quantile_data %>% group_by(treatment, temp, light, connectivity, state_var) %>% 
  summarise( location = which.max(quant25)) %>% 
  group_by(treatment, temp, light, connectivity) %>% 
  summarise( time_zoom = max(location)) %>% 
  ungroup() %>% mutate(time_zoom = ifelse(time_zoom == 1, this_param_set[["tmax"]], 
                                          ceiling((time_zoom*this_param_set[["dt"]]*quantile_thinning)*1.1))) %>% 
  right_join(quantile_data, by = c("treatment", "temp", "light", "connectivity")) %>% 
  filter(time_zoom >= time)

# Look at state variables quantiles together over time 
p3 <- ggplot(quantile_data_zoomed) + 
  aes(x = time/365, y = med_value, colour = state_var) + 
  geom_line() + 
  geom_line(aes(x = time/365, y = quant25), linetype="dashed") +
  geom_line(aes(x = time/365, y = quant75), linetype="dashed") + 
  theme_bw() +
  facet_wrap(.~ treatment + connectivity, scales = "free", nrow = 4) + 
  geom_hline(yintercept = 0, linetype="dashed" ) + ylab("Biomass (g dry wt m-2)") +
  xlab("Time (Years)") 

p3 

if (SAVE) { # Save QUANTILES OVER TIME ZOOMED if specified
  ggsave(paste0("Figures/", filename_details,"_overtime_zoom_25-50-75quantile.png"), width = 10, height = 10)
}

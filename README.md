# Coastal habitat recovery models

Code to support the manuscript: "Speeding up the recovery of coastal habitats through management interventions that address constraints on dispersal and recruitment"

Christopher J. Brown1,2, Max D. Campbell2, Catherine J. Collier3, Mischa P. Turschwell2, Megan I. Saunders4, Rod M. Connolly2

Code Developers: Max Campbell, Dr Chris Brown

# Title: Connectivity facilitates restoration of foundation habitats in the presence of multiple stressors

## Summary

Habitat restoration plans would benefit from predictions of timescales for recovery. Theoretical models have been a powerful tool for informing practical guidelines in planning marine protected areas, suggesting restoration planning could also benefit from a theoretical framework. We developed a model that can predict recovery times following restoration action, under dispersal, connectivity and recruitment constraints. We apply the model to a case-study of seagrass restoration and find recovery times following restoration action can vary greatly, from <1 to >20 years. The model also shows how recovery can be accelerated when restoration actions are matched to the constraints on recovery. For example, spreading propagules can be used when connectivity is restricted. The recovery constraints we articulated mathematically also apply to the restoration of coral reefs, mangroves, saltmarsh, shellfish reefs and macroalgal forests, so our model provides a general framework for choosing restoration actions that accelerate coastal habitat recovery. 


## How to navigate the code

The code is is two locations "Scripts" and "Functions", where the user should always work from the "Scripts" folder and these scripts source the functions they need to run the desired analysis.  The number prefixes indicate the logical order you should run the self contained scripts in.

### Scripts/

**halodule-simulations.R** and **zostera-simulations.R** - evaluate the models and creates a file with all model output for plotting

**halodule-simulations-sensitivity.R** and **zostera-simulations-sensitivity.R** - as above, but also runs models for different parameter settings to create sensitivity analysis in the supplemental files. 

**zostera-params.R** and **halodule-params.R** Zostera and halodule parameter settings.

Other scripts create the figures in the manuscript as per their labels. 


### Functions/

**Analytical_solution_functions.R** - analytical solutions and computations of the deterministic part of the models 

**Conversion_functions.R** - Conversions of units for the parameters

**photosynthesis_function.R** - Function for photosynthesis rate given temperature and light

**Recovery_times_from_recolonisation.R** - Analytical solution for the recovery times of patches from the beginning of recolonisation

**resp_function.R** - Function for respiration rate given temperature

**Run_DE_system_base.R** - Version to solve the dynamical system numerically (OLD)

**Seeding_functions.R** - Functions for probability of recolonation

**Stochastic_dispersal_model.R** - The dispersal model used for both experiments

**Stochastic_recruitment_model.R** - The recruitment model used for both experiments

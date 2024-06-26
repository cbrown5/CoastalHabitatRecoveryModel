# Coastal habitat recovery models

R code to support the manuscript: "Speeding up the recovery of coastal habitats through management interventions that address constraints on dispersal and recruitment"

Christopher J. Brown, Max D. Campbell, Catherine J. Collier, Mischa P. Turschwell, Megan I. Saunders, Rod M. Connolly

Code Developers: Dr Chris Brown, Max Campbell

2024-06-26

### Recommended citation: 

It is recommended you look up the peer-reviewed study by the same title and cite that. 

Brown, Campbell, 2024. R code to support the manuscript: "Speeding up the recovery of coastal habitats through management interventions that address constraints on dispersal and recruitment" Zenodo. DOI: 10.5281/zenodo.12541372

[![DOI](https://zenodo.org/badge/712717461.svg)](https://zenodo.org/doi/10.5281/zenodo.12541371)


## How to navigate the code

The code is is two locations "Scripts" and "Functions", where the user should always work from the "Scripts" folder and these scripts source the functions they need to run the desired analysis.  The number prefixes indicate the logical order you should run the self contained scripts in.


## Overview and Aims

Create a model that can predict recovery time of coastal habitats under ecological constraints on dispersal, management and connectivity to propagule sources. Investigate the influence of management efforts on reducing recovery time. 

## Data sources

No data was used in this study. Model parameters are taken from the literature and are descriped in the peer-reviewed study. 

### Scripts/

All code was written in R with version 4.4.0. The R packages `tidyverse`, `patchwork` and `plotly` are used for plotting and data wrangling. 

**zostera-params.R**  **halodule-params.R** 

Zostera and Halodule parameter settings.

**halodule-simulations.R**  **zostera-simulations.R** 

Evaluate the models and creates a file with all model output for plotting

**zostera-simulations-sensitivity.R**
 
As above, but also runs models for different parameter settings to create sensitivity analysis in the supplemental files. 

**zostera-simulations-mgmt-v2.R** **halodule-simulations-mgmt-v2**

Runs Zostera or Halodule model with management scenarios, input to create fig 5. 

**zostera-simulations-mgmt-phi-sens.R** 

Tests sensitivity of management benefits (time saved) to the parameter phi. For supplementary figures. 

**zostera-simulations-light-sensitivity.R** **halodule-simulations-light-sensitivity.R**

Sensitivity analysis of Zostera and Halodule to different stressor scenarios. Not used in final revised ms, but kept here for future use. 

**figure2-colonisation-prob.R** 

Figure of colonisation probability for the dispersal and recruitment models. 

**figure3-connections.R**

Figure of recovery time under different levels of connectivity. 

**figure4-connections-stressors.R**

Figure of recovery time under different levels of connectivity and stressors.

**figure5-mgmt-comparisons-v2.R**

Figure of recovery time under different management scenarios and light levels. 

**figure-param-sensitivity.R**

Creates figure for supp with the sensitivity of recovery time to key parameters. 
**figure-mgmt-comparisons-phisens-v2.R** 

Creates figure for supp with the sensitivity of management benefits to phi.


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

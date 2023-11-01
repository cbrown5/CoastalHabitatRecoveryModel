#Plots of simulations 
rm(list = ls())
library(tidyverse)
library(patchwork)
# source functions 
source("Functions/Conversion_functions.R")
source("Functions/Stochastic_dispersal_model.R")
source("Functions/Stochastic_recruitment_model.R")
source("Functions/Seeding_functions.R")
source("Scripts/zostera_params.R")
load(file = "Outputs/zostera-recruit-function-sensitivity.rda")

treatdf <-  treatment_df %>% 
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Dispersal limited",
    grepl("recolon", model) ~ "Recruitment limited"),
    `# sources` = factor(n_sources))



# ------------ 
# PLOTS 
# ------------ 
plotout <- NULL
for (imodels in c("disp", "rec")){
  g1_l_max <- filter(treatdf, 
                     nu == zostera_param_set$nu &
                       phi == zostera_param_set$phi) %>%
    filter(n_sources %in% c(1,5)) %>%
    filter(grepl(imodels, model)) %>%
    ggplot() + 
    aes(x = l_max, y = median_recovery/365, colour = `# sources`,
        shape = `# sources`) +
    geom_hline(yintercept = 0) +
    theme_classic()+
    geom_line(position = position_dodge(width = 0.1), alpha = 0.5) + 
    geom_point(aes(shape = `# sources`),
               position = position_dodge(width = 0.1))  +
    geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                   position = position_dodge(width = 0.1), size = 1, alpha = 0.5) +
    geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                   position = position_dodge(width = 0.1), alpha = 0.5) +
    facet_wrap(~ model) + 
    scale_color_manual(values = c("darkblue", "tomato")) + 
    scale_y_log10() +
    # scale_shape_manual(values = c(1,2,3)) + 
    xlab(expression("Max probability "(L[m][a][x]))) + 
    ylab("") + 
    # scale_x_continuous(breaks = 1:10) +
    # ylim(0, ymax) +
    # scale_color_manual(values = c("black", "#d41515")) + 
    theme(axis.text = element_text(size = 11), 
          axis.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()
          # legend.title = element_blank()
    )

  pd <- position_dodge(width = 0.001)
  g1_nu <- filter(treatdf, 
                  l_max == zostera_param_set$l_max &
                    phi == zostera_param_set$phi) %>%
    filter(n_sources %in% c(1,5)) %>%
    filter(grepl(imodels, model)) %>%
    ggplot() + 
    aes(x = nu, y = median_recovery/365, 
        colour = `# sources`,
        shape = `# sources`) +
    geom_point(position = pd)  +
    geom_hline(yintercept = 0) +
    geom_line(position = pd, alpha = 0.5) + 
    geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                   position =pd, size = 1, alpha = 0.5) +
    geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                   position = pd, alpha = 0.5) +
    facet_wrap(~ model) + 
    scale_color_manual(values = c("darkblue", "tomato")) +
    scale_y_log10() +
    # scale_shape_manual(values = c(1,2,3)) + 
    xlab(expression("Steepness "(nu))) + 
    ylab("Recovery time (years)") + 
    theme_classic()+
    # scale_x_continuous(breaks = 1:10) +
    # ylim(0, ymax) +
    # scale_color_manual(values = c("black", "#d41515")) + 
    theme(axis.text = element_text(size = 11), 
          axis.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()
          # legend.title = element_blank()
    )

  
  pd <- position_dodge(width = 0.001)
  g1_phi <- filter(treatdf, 
                   l_max == zostera_param_set$l_max &
                     nu == zostera_param_set$nu) %>%
    filter(n_sources %in% c(1,5)) %>%
    filter(grepl(imodels, model)) %>%
    ggplot() + 
    aes(x = phi, y = median_recovery/365, 
        colour = `# sources`,
        shape = `# sources`) +
    geom_point(position = pd)  +
    geom_hline(yintercept = 0) +
    geom_line(position = pd, alpha = 0.5) + 
    geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                   position =pd, size = 1, alpha = 0.5) +
    geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                   position = pd, alpha = 0.5) +
    facet_wrap(~ model) + 
    scale_color_manual(values = c("darkblue", "tomato")) + 
    scale_y_log10() +
    # scale_shape_manual(values = c(1,2,3)) + 
    xlab(expression("Biomass at 50% "(tau))) + 
    ylab("") + 
    theme_classic()+
    # scale_x_continuous(breaks = 1:10) +
    # ylim(0, ymax) +
    # scale_color_manual(values = c("black", "#d41515")) + 
    theme(axis.text = element_text(size = 11), 
          axis.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()
          # legend.title = element_blank()
    ) 
plotout <- c(plotout, list(g1_l_max),
             list(g1_nu),
             list(g1_phi))
  
  
}
plotout[[6]] <- plotout[[6]] +
  theme(legend.position = c(0.2, 0.75))

gall <- (plotout[[1]] + plotout[[4]]) / 
 ( plotout[[2]] + plotout[[5]]) /
 ( plotout[[3]] + plotout[[6]]) +
  plot_annotation(tag_levels = "A")
gall

ggsave(gall, 
       file ="Outputs/recovery-time-sensitivity.png",
       width = 6, height = 8)


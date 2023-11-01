#Conceptual figure

library(ggplot2)
library(patchwork)
source("Functions/Seeding_functions.R")

nu = 0.01 # steepness default 0.01 
phi = 500 # phase shift of curve default 500
l_max = 0.2 

x <- seq(0, 2000, by = 1)

# Dispersal limitation model
y <- seeding_old_function(x = x, phi = 500, 
                            l_max = 0.2, nu = 0.01) 
y2 <- seeding_old_function(x = x, phi = 500, 
                             l_max = 0.4, nu = 0.02) 

# Recruitment limitation model
recprob <- seeding_old_function(x = x, phi = 500, 
                                l_max = 0.38, nu = 0.01) 
recprob2 <- seeding_old_function(x = x, phi = 500, 
                                 l_max = 0.38*2, nu = 0.02) 


#
# GGPLOT2 version 
#

dat <- data.frame(x= c(x,x,x,x), 
                  y =c(y, y2,
                       recprob, recprob2),
                  Scenario = rep(rep(c("No stressors",
                                       "Stressors"), 
                                     each = length(x)),2),
                  Model = rep(c("Dispersal limited",
                                "Recruitment limited"),
                              each = 2*length(x)))

pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2")

g1 <- ggplot(subset(dat, Model == "Dispersal limited")) +
  aes(x = x, y = y, color = Scenario) +
  geom_line() +
  ylab("Probability of dispersal") +
  xlab("Patch biomass \n (proportional to propagule supply)") +
  theme_classic() +
  scale_color_manual(values = c( "black", "hotpink3"))+
  theme(legend.position = "none")
g1

g2 <- ggplot(subset(dat, Model == "Recruitment limited")) +
  aes(x = x, y = y, color = Scenario) +
  # geom_vline(xintercept = 1.3, linetype = 2,
  #            color = "darkblue") +
  # geom_vline(xintercept = 4, linetype = 2,
  #            color = "red") +
  geom_line() +
  ylab("Probability of recruitment") +
  xlab("Total connected biomass \n (proportional to propagule availability)") +
  theme_classic() +
  scale_color_manual(values =  c( "black", "hotpink3"))  +
  theme(legend.position = "none")
g2

#
# Save figures
#

ggsave("Outputs/conceptual-fig1-dispmodel.png", g1,
       width = 4, height = 3)
ggsave("Outputs/conceptual-fig1-recmodel.png", g2,
       width = 4, height = 3)

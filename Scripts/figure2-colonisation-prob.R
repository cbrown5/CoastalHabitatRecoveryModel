#Conceptual figure

library(ggplot2)
library(patchwork)
source("Functions/Seeding_functions.R")

nu = 0.01 # steepness default 0.01 
phi = 500 # phase shift of curve default 500
l_max = 0.2 

x <- seq(0, 10, by = 0.1)

# Dispersal limitation model
y <- 1-seeding_old_function(x = 750, phi = 500, 
                            l_max = 0.2, nu = 0.01) 
dispprob <- 1-(y^x)
y2 <- 1-seeding_old_function(x = 250, phi = 500, 
                             l_max = 0.2, nu = 0.01) 
dispprob2 <- 1-(y2^x)

# Recruitment limitation model
recprob <- seeding_old_function(x = 750*x, phi = 500, 
                                l_max = 0.2, nu = 0.01) 
recprob2 <- seeding_old_function(x = 250*x, phi = 500, 
                                 l_max = 0.2, nu = 0.01) 


#
# GGPLOT2 version 
#

dat <- data.frame(x= c(x,x,x,x), 
                  y =c(dispprob, dispprob2,
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
  ylim(0, 1)+
  ylab("Probability of \n colonisation") +
  xlab("Number of connected patches") +
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  theme_classic() +
  scale_color_manual(values = c( "Darkblue", "Red")) +
  theme(legend.position = "none")
g1

g2 <- ggplot(subset(dat, Model == "Recruitment limited")) +
  aes(x = x, y = y, color = Scenario) +
  # geom_vline(xintercept = 1.3, linetype = 2,
  #            color = "darkblue") +
  # geom_vline(xintercept = 4, linetype = 2,
  #            color = "red") +
  geom_line() +
  ylim(0, 1)+
  ylab("") +
  xlab("Number of connected patches") +
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  theme_classic() +
  scale_color_manual(values = c( "Darkblue", "Red"))  +
  theme(legend.position = c(0.6, 0.9))
g2

gboth <- g1 + g2 + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")
gboth
ggsave("Outputs/colonisation prob.png", 
       gboth,
       width = 7, height = 3)



dat2 <- data.frame(x = 0:10, y = 1:11)

g3 <- ggplot(dat2) +
  aes(x = x, y = y) +
  geom_line() +
  ylab("Number of \n patches disturbed") +
  xlab("Connectivity") +
  theme_classic() #+
theme(axis.text.x = element_blank(),
      axis.text.y  = element_blank())

#
# Recovery times 
#

x3 <- seq(0, 10, by = 0.1)
dat3 <- data.frame(x = c(x3,x3), 
                   y = c(exp(x3*-1),
                         exp(x3*-0.2)),
                   Scenario = rep(c("No stressors",
                                    "Stressors"), 
                                  each = length(x3)))

g4 <- ggplot(dat3) +
  aes(x = x, y = y, color = Scenario) +
  geom_line() +
  ylab("Recovery time") +
  xlab("Connectivity") +
  theme_classic() +
  scale_color_manual(values = c( "Darkblue", "Red"))  +
  theme(axis.text.y  = element_blank(),
      legend.position = "none")
g4

dat4 <- data.frame(x = c(x3,x3), 
                   y = c(exp((x3+1)*-2) +0.08*x3,
                         exp(x3*-1) + 0.08*x3),
                   Scenario = rep(c("No stressors",
                                    "Stressors"), 
                                  each = length(x3)))

g5 <- ggplot(dat4) +
  aes(x = x, y = y, color = Scenario) +
  geom_line() +
  ylab("Recovery time") +
  xlab("Connectivity") +
  theme_classic() +
  scale_color_manual(values = c( "Darkblue", "Red"))  +
   theme(axis.text.y  = element_blank(),
        legend.position = "none")
g5

#
# Save figures
#

ggsave("Outputs/conceptual-fig-dispmodel.png", g1,
       width = 3, height = 2)
ggsave("Outputs/conceptual-fig-recmodel.png", g2,
       width = 3, height = 2)
ggsave("Outputs/conceptual-fig-disturbance.png", g3,
       width = 3, height = 2)
ggsave("Outputs/conceptual-fig-disp-recovery.png", g4,
       width = 3, height = 2)
ggsave("Outputs/conceptual-fig-rec-recovery.png", g5,
       width = 3, height = 2)


#
# Probability vs biomass
#

x <- seq(0, 600, by = 1)

# Dispersal limitation model
y <- seeding_old_function(x = x, phi = 500, 
                            l_max = 0.38, nu = 0.01) 

plot(x,y, ylim = c(0, 0.4))

seeding_old_function(
                     x = 558*0.5, 
                     # x = 667*0.2,
                     phi = 500, 
                     l_max = 0.38, nu = 0.01)

# Recruitment limitation model

phivals <- seq(1, 1000, by = 10)
nuvals <- seq(0.005, 0.05, length.out = 500)
vals <- expand.grid(phivals = phivals, nuvals = nuvals)
vals$prob_bmax <- seeding_old_function(
  x = 558, 
  phi = vals$phivals, 
  l_max = 0.38, nu = vals$nuvals)
vals$prob_20 <- seeding_old_function(
  x = 558*0.2, 
  phi = vals$phivals, 
  l_max = 0.38, nu = vals$nuvals)


ggplot(vals) + 
  aes(x = phivals, y = nuvals, fill = prob_bmax) + 
  geom_raster()

ggplot(vals) + 
  aes(x = phivals, y = nuvals, fill = prob_20) + 
  geom_raster()

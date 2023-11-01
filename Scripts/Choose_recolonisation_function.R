source("Functions/Seeding_functions.R")
nu = 0.01 # steepness default 0.01 
phi = 500 # phase shift of curve default 500
l_max = 0.2 

x <- 1:2000
y <- seeding_old_function(x = x, phi = 500, l_max = 0.2, nu = 0.01) 
y1 <- seeding_old_function(x = x, phi = 500, l_max = 0.2, nu = 0.02) 
y2 <- seeding_old_function(x = x, phi = 500, l_max = 0.2, nu = 0.005) 

y1 <- seeding_old_function(x = x, phi = 500, l_max = 0.4, nu = 0.01) 
y2 <- seeding_old_function(x = x, phi = 500, l_max = 0.1, nu = 0.01) 

y1 <- seeding_old_function(x = x, phi = 1000, l_max = 0.2, nu = 0.01) 
y2 <- seeding_old_function(x = x, phi = 2000, l_max = 0.2, nu = 0.01) 

y1 <- seeding_old_function(x = x, phi = 1000, l_max = 0.2, nu = 0.005) 
y2 <- seeding_old_function(x = x, phi = 250, l_max = 0.2, nu = 0.02) 
y3 <- seeding_old_function(x = x, phi = 2000, l_max = 0.2, nu = 0.0025) 


# Zostera cases
y <- seeding_old_function(x = x, phi = 500, l_max = 0.2, nu = 0.01) 
y1 <- seeding_old_function(x = x, phi = 750, l_max = 0.2, nu = 0.007) 
y2 <- seeding_old_function(x = x, phi = 250, l_max = 0.2, nu = 0.02) 
y3 <- seeding_old_function(x = x, phi = 125, l_max = 0.2, nu = 0.035) 

png(filename = "Outputs/seeding_functions_for_sensitivity.png", 
    width =32, height = 15, units = "cm", res = 300)
par(mfrow = c(1,2))
plot(y = y, x = x, type = "l", ylab = "Probability of recruitment", xlab = "Total connected biomass (indicative of seeds in seedbank)")
lines(y = y1, x = x, col = "blue")
lines(y = y2, x = x, col = "red")
lines(y = y3, x = x, col = "purple")

plot(y = y, x = x, type = "l", ylab = "Probability of dispersal (single patch)", xlab = "Patch Biomass (indicative of the number of seeds dispersed)")
lines(y = y1, x = x, col = "blue")
lines(y = y2, x = x, col = "red")
lines(y = y3, x = x, col = "purple")
legend(x = 1700, y = 0.05, legend = c("phi = 500, nu = 0.01 (base)","phi = 1000, nu = 0.005", 
                  "phi = 250, nu = 0.02","phi = 2000, nu = 0.0025"),
       fill = c("black", "blue", "red", "purple"), cex = 0.75, bty = "n")

dev.off()

#
# GGPLOT2 version 
#
library(ggplot2)
library(patchwork)

dat <- data.frame(x = rep(x, 4), 
                  scnr = rep(c("phi = 500, nu = 0.01 (base)","phi = 1000, nu = 0.005", 
                               "phi = 250, nu = 0.02","phi = 2000, nu = 0.0025"),
                             each = length(x)), 
                  y = c(y, y1, y2, y3))

pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2")

g1 <- ggplot(dat) +
  aes(x = x, y = y, color = scnr, linetype = scnr) +
  geom_hline(yintercept = 0.2) +
  geom_line() +
  ylab("Probability of dispersal") +
  xlab("Patch Biomass \n (indicative of propagule supply)") +
  theme_classic() +
  scale_color_manual("Parameters", 
                     values = pal[2:5]) +
  scale_linetype("Parameters")

g2 <- ggplot(dat) +
  aes(x = x, y = y, color = scnr, linetype = scnr) +
  geom_hline(yintercept = 0.2) +
  geom_line() +
  ylab("Probability of recruitment") +
  xlab("Total connected biomass \n (indicative of propagule availability)") +
  theme_classic() +
  scale_color_manual("Parameters", 
                     values = pal[2:5]) +
  scale_linetype("Parameters")

gall <- g1 + g2 + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides='collect')
ggsave("Outputs/ggplot_seeding_functions_for_sensitivity.png",
       width = 8, height = 3)

#
# Lmax sens
#

chance_rec <- function(l_max, years){
  1 - (1 - l_max)^years
}

png(filename = "Outputs/l_max_sensitivity.png", width =17, height = 15, units = "cm", res = 300)
par(mfrow = c(1,1))
years_vect <- seq(0, 20, by = 0.25)
plot(y = chance_rec(0.38, years_vect), x = years_vect, type = "l", 
     ylab = "Max probability recovered", xlab = "Years")
lines(y = chance_rec(0.19, years_vect), x = years_vect, col = "purple")
lines(y = chance_rec(0.76, years_vect), x = years_vect, col = "blue")
legend(x = 10, y = 0.25, legend = c("lmax = 0.38 (base)","lmax = 0.19", "lmax = 0.76"),
       fill = c("black", "purple", "blue"), cex = 01, bty = "n")

dev.off()

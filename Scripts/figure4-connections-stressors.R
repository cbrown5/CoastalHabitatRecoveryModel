#Plots of simulations 

library(tidyverse)
library(patchwork)

load(file = "Outputs/zostera-stressors-scenarios.rda")
load(file = "Outputs/halodule-stressors-scenarios.rda")

treatdf <- rbind(treatment_df, halo_treatment_df)

# ------------ 
# PLOTS 
# ------------ 

species_to_plot <- "zost"
# species_to_plot <- "halo"


#
# Figure 3, just two models, noe for each species 
#

mod2_res <- treatdf %>% 
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Uncertain dispersal",
    grepl("recolon", model) ~ "Uncertain recruitment"),
    Light = case_when(I. < 60 ~ "Extreme stress - 25% light",
                      I. > 60  & I. < 200 ~ "Stressed - 50% light", 
                      I. > 200 ~ "Normal"))
mod2_res$Light <- factor(mod2_res$Light, levels = c("Extreme stress - 25% light",
                                                    "Stressed - 50% light",
                                                    "Normal"))
ymax <- 50
g1_Disp <- 
  filter(mod2_res, grepl(species_to_plot, species)) %>%
  filter(., grepl("dispersal", model)) %>%
  ggplot() + 
  aes(x = n_sources, y = median_recovery/365, colour = Light,
      shape = Light) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Light), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1,2, 3)) +
  xlab("Number of connected patches") + 
  ylab("Recovery time (years)") + 
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c( "#d41515",  "darkblue","black")) + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = "none",
        legend.title = element_blank()
  )
g1_Disp

g1_rec <- 
  filter(mod2_res, grepl(species_to_plot, species)) %>%
  filter(., grepl("rec", model)) %>%
  ggplot() + 
  aes(x = n_sources, y = median_recovery/365, colour = Light,
      shape = Light) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Light), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1,2, 3)) +
  xlab("Number of connected patches") + 
  ylab("") +
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c( "#d41515",  "darkblue","black")) + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.title = element_blank(),
        legend.position = c(.65, .80)
  )

g1_rec


gboth <- g1_Disp + g1_rec +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")
gboth


ggsave(gboth, 
       file =paste0("Outputs/fig-4-recovery-time-",species_to_plot,"-light-stress.png"),
       width = 7, height = 3)



#Results for 1 source. 
filter(mod2_res, n_sources == 1) %>% mutate(med = median_recovery/365) %>% View()




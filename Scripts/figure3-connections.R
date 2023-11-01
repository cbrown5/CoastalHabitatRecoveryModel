#Plots of simulations 

library(tidyverse)
library(patchwork)

load(file = "Outputs/zostera-stressors-scenarios.rda")
load(file = "Outputs/halodule-stressors-scenarios.rda")

treatdf <- rbind(treatment_df, halo_treatment_df)

# ------------ 
# PLOTS 
# ------------ 

#
# Figure 3, just two models, noe for each species 
#

mod2_res <- treatdf %>% 
  filter(Stressors == "Control") %>%
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Uncertain dispersal",
    grepl("recolon", model) ~ "Uncertain recruitment"))

ymax <- 12
g1_zost <- ggplot(filter(mod2_res, grepl("zost", species))) + 
  aes(x = n_sources, y = median_recovery/365, colour = Model,
      shape = Model) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Model), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1,2)) + 
  xlab("Number of connected patches") + 
  ylab("Recovery time (years)") + 
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c("black", "#d41515")) + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = "none",
        legend.title = element_blank()
  )


g1_halo <- ggplot(filter(mod2_res, grepl("halo", species))) + 
  aes(x = n_sources, y = median_recovery/365, colour = Model,
      shape = Model) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Model), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1,2)) + 
  xlab("Number of connected patches") + 
  ylab("") + 
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c("black", "#d41515")) + 
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = c(0.65,0.8),
        legend.title = element_blank()
  )

gboth <- g1_zost + g1_halo +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")
gboth


ggsave(gboth, 
       file ="Outputs/recovery-time-both-species.png",
       width = 7, height = 3)



#
# Effect of stressors
#

ggplot(treatment_df) + aes(x = n_sources, y = median_recovery/365, colour = Stressors) +
  geom_hline(yintercept = 0) +
  geom_point(aes(shape = Stressors), position = position_dodge(width = 0.6)) + theme_classic()+
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(1,2,4,5)) + 
  facet_wrap(.~ model, nrow = 2, scales = "free_y") + 
  xlab("Number of source patches") + 
  scale_x_continuous(breaks = 1:10) +
  ylab("Recovery time (years)") + 
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), legend.position = c(0.8,0.3),
          legend.title = element_blank()
        )
  
ggsave(sprintf("Shared/Figures/%s_recovery_time_vs_npatches.png", 
               Experiment_name), width = 18, height = 15, units = "cm", dpi = 300, 
       scale = 0.9)


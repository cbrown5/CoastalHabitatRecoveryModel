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
  filter(grepl("zost", species)) %>%
  # filter(n_sources %in% c(1, 5)) %>%
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Dispersal limited",
    grepl("recolon", model) ~ "Recruitment limited"),
    `# sources` = factor(n_sources))

ymax <- 20
g1_disp <-
  ggplot(filter(mod2_res, grepl("dispersal", model))) +
# ggplot(mod2_res) +
  aes(x = n_sources, y = median_recovery/365, 
      colour = Stressors,
      shape = Stressors) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Stressors), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1:5)) +
  xlab("Number of connected meadows") + 
  ylab("Recovery time (years)") + 
  scale_x_continuous(breaks = seq(0,10,by=2)) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c("blue", "black", "grey", "#d41515")) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = "none",
        legend.title = element_blank()
  )
g1_disp

g1_rec <-
  ggplot(filter(mod2_res, grepl("recolon", model))) +
  # ggplot(mod2_res) +
  aes(x = n_sources, y = median_recovery/365, 
      colour = Stressors,
      shape = Stressors) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_line(position = position_dodge(width = 0.6), alpha = 0.5) + 
  geom_point(aes(shape = Stressors), position = position_dodge(width = 0.6))  +
  geom_linerange(mapping = aes(ymin = per25/365, ymax = per75/365), 
                 position = position_dodge(width = 0.6), size = 1, alpha = 0.5) +
  geom_linerange(mapping = aes(ymin = per5/365, ymax = per95/365), 
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  scale_shape_manual(values = c(1:5)) +
  xlab("Number of connected meadows") + 
  ylab("") + 
  scale_x_continuous(breaks = seq(0,10,by=2)) +
  coord_cartesian(ylim=c(0, ymax))+
  scale_color_manual(values = c("blue", "black", "grey", "#d41515")) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        legend.position = c(0.6, 0.8),
        legend.title = element_blank()
  )

gboth <- g1_disp + g1_rec +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")
gboth
ggsave(gboth, 
       file ="Outputs/recovery-time-stressors-zostera.png",
       width = 7, height = 3)




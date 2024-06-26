# Plot of mgmt effects on recovery times
#
# CJ Brown 2024-02-19

# Figures for supplemental material

library(tidyverse)
library(patchwork)
rm(list=ls())
load(file = "Outputs/zostera-mgmt-results-phi-sens.rda")
# load(file = "Outputs/halodule-mgmt-results.rda")

#all light levels, two phi parameters
baseline <- filter(all_param_combos, 
                   phi > 20 & l_max<0.4) %>%
  left_join(treatment_df) 

#
#MGMT SCENARIOS
#

#increase light (just high light levels)
light <- filter(all_param_combos, 
                  I. < 200 & phi > 20 & l_max<0.4) %>%
  left_join(treatment_df) %>%
  mutate(Scenario = "+ light") %>%
  select(-I.)

unique(treatment_df$model)

#Sediment stabilisation
#Filter for: 
# Rec model with perfect parameters
sed_stabil1 <- filter(all_param_combos, 
                phi < 21 & l_max>0.4) %>%
  mutate(model = "Stochastic recolonisation") %>%
  inner_join(treatment_df)%>%
  mutate(Scenario = "+ sed. stabilisation")

# Dispersal limited model, basic parameters
sed_stabil2 <- filter(all_param_combos, 
                     phi > 20 & l_max<0.4) %>%
  mutate(model = "Stochastic dispersal") %>%
  inner_join(treatment_df) %>%
  mutate(Scenario = "+ sed. stabilisation")

#Dispersal mgmt
#Filter for: 
# Disp model with perfect parameters, max connectivity
disp1 <- filter(all_param_combos, 
                      phi < 21 & l_max>0.4) %>%
  mutate(model = "Stochastic dispersal") %>%
  #remove these, so it joins to all of them
  inner_join(treatment_df)%>%
  filter(n_sources == 10) %>%
  mutate(Scenario = "+ seeds/plant")

#
# Recruitment limited model, basic parameters, max connectivity
#
disp2 <- filter(all_param_combos, 
                phi > 20 & l_max<0.4) %>%
  mutate(model = "Stochastic recolonisation") %>%
  #remove these, so it joins to all of them
  inner_join(treatment_df)%>%
  filter(n_sources == 10) %>%
  mutate(Scenario = "+ seeds/plant")

#
# MGMT COMPARISONS
#

# Need to join baseline with each scenario. The columns to join by 
# vary, depending on the scenario and model. 

light_mgmt <- left_join(baseline, light, by = c("phi", "n_sources", "n", "model")) %>%
  mutate(`Gain (years)` = (median_recovery.x - median_recovery.y)/365,
         `Gain (times faster)` = median_recovery.x / median_recovery.y,
         median_recovery = median_recovery.y) %>%
  select(model, Scenario, n_sources, I., phi, `Gain (years)`, `Gain (times faster)`,median_recovery)

sed_mgmt1 <- inner_join(baseline, select(sed_stabil1, -phi), 
                      by = c("n_sources", "n", "model", "I.")) %>%
  mutate(`Gain (years)` = (median_recovery.x - median_recovery.y)/365,
         `Gain (times faster)` = median_recovery.x / median_recovery.y,
         median_recovery = median_recovery.y) %>%
  select(model, Scenario, n_sources, I., phi, `Gain (years)`, `Gain (times faster)`,median_recovery)

sed_mgmt2 <- inner_join(baseline, sed_stabil2, 
                        by = c("phi", "n_sources", "n", "model", "I.")) %>%
  mutate(`Gain (years)` = (median_recovery.x - median_recovery.y)/365,
         `Gain (times faster)` = median_recovery.x / median_recovery.y,
         median_recovery = median_recovery.y) %>%
  select(model, Scenario, n_sources, I., phi, `Gain (years)`, `Gain (times faster)`,median_recovery)

disp_mgmt1 <- inner_join(baseline, select(disp1, -phi, -n_sources, -n), 
                        by = c("model", "I.")) %>%
  mutate(`Gain (years)` = (median_recovery.x - median_recovery.y)/365,
         `Gain (times faster)` = median_recovery.x / median_recovery.y,
         median_recovery = median_recovery.y) %>%
  select(model, Scenario, n_sources, I., phi, `Gain (years)`, `Gain (times faster)`,median_recovery)

disp_mgmt2 <- inner_join(baseline, select(disp2, -n_sources, -n), 
                        by = c("phi", "model", "I.")) %>%
  mutate(`Gain (years)` = (median_recovery.x - median_recovery.y)/365,
         `Gain (times faster)` = median_recovery.x / median_recovery.y,
         median_recovery = median_recovery.y) %>%
  select(model, Scenario, n_sources, I., phi, `Gain (years)`, `Gain (times faster)`,median_recovery)

mgmt_scnrs <- bind_rows(light_mgmt, sed_mgmt1, sed_mgmt2, disp_mgmt1, disp_mgmt2) %>%
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Uncertain dispersal",
    grepl("recolon", model) ~ "Uncertain recruitment"),
    conn = ifelse(n_sources == 1, "Low connectivity", "High connectivity")
  ) %>%
  mutate(light_lev = I./max(I.))

# ------------ 
# PLOTS 
# ------------ 

theme_set(theme_classic())

pd <- position_dodge(width = 0.75)

g1 <- ggplot(mgmt_scnrs) +
  aes(x = n_sources, y = median_recovery/365, color = Scenario) +
  geom_hline(yintercept = 0) +
  geom_point(position = pd) +
  facet_grid(phi~Model, scales = "free") +
  # ylim(0, 10)+
  ylab("Recovery time (years)") +
  xlab("Number of source patches") + 
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c("blue", "darkgreen",
                                "purple", "black"))

g1

  
# ggsave(g1, file = "Outputs/recovery-time-mgmt-zostera.png",
       # width = 6, height = 4)
# ggsave(g1, file = "Outputs/recovery-time-mgmt-halodule.png",
       # width = 6)

#Note if we wanted to plot difference between control and 
# the scenarios we need to recalculate variacnes from lambda
# then sum across those. 

#
#Above but as a barplot of mgmt gain
#


g3 <- 
  filter(mgmt_scnrs, n_sources %in% c(1,8) & 
          light_lev<0.6) %>%
  ggplot() +
  aes(x = phi, 
      # y = `Gain (years)`,
      y = `Gain (times faster)`,
      color = conn,
      linetype = conn) +
  geom_hline(yintercept = 1) +
  geom_line(size = 1, alpha = 0.8) +
  facet_grid(Scenario~Model) +
  xlab("Biomass productivity parameter") + 
  # xlab(expression("Biomass at 50% "(tau))) +
  ylab("Benefit from management\n action (relative increase in recovery rate)")+
  scale_color_manual(values = c("lightblue",
                               "red")) + 
  scale_linetype_manual(values = c(1,2))
g3

plotly::ggplotly(g3)

ggsave(g3, file = "Outputs/mgmt-gain-zostera-curves-tau-sens.png",
       width = 6, height = 3.5)
# ggsave(g3, file = "Outputs/mgmt-gain-halodule-bars.png",
       # width = 6, height = 3.5)




#Plots of simulations 

library(tidyverse)
library(patchwork)
rm(list=ls())
load(file = "Outputs/zostera-mgmt-results.rda")
# load(file = "Outputs/halodule-mgmt-results.rda")

# treatdf <- rbind(treatment_df, halo_treatment_df)

#
# Setup mgmt comparisons 
#
#All scenarios have low light, unless its specifically managed
# all scenarios repeated for high temp and normal temp
LT_LL <- filter(treatment_df, T. == 30.3 & 
                       I. < 200 & 
                       l_max == 0.38 &
                       phi == 500) %>%
  select(n_sources, model, T.,  per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "No mgmt")

HT_LL <- filter(treatment_df, T. > 31 & 
                       I. <  200 & 
                       l_max == 0.38 &
                       phi == 500)%>%
  select(n_sources, model, T.,  per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "No mgmt")

HT_manage_light <- filter(treatment_df, T. > 31 & 
                             I. >  200 & 
                             l_max == 0.38 &
                             phi == 500)%>%
  select(n_sources, model, T.,  per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ light")

manage_light <- filter(treatment_df, T. < 31 & 
                            I. >  200 & 
                            l_max == 0.38 &
                            phi == 500) %>%
  select(n_sources, model, T.,  per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ light")

#Rec mgmt
#Combine: 
# Rec model with perfect parameters
# Dispersal limited model, basic parameters
HT_manage_recruitment <- filter(treatment_df, T. > 31 & 
                                I. <  200 & 
                                l_max > 0.5 &
                                phi < 400 & 
                                model == "Stochastic recolonisation") %>%
  bind_rows(filter(treatment_df, T. > 31 & 
                     I. <  200 & 
                     l_max < 0.5 &
                     phi > 400 & 
                     model == "Stochastic dispersal")) %>% 
  select(n_sources,  model, T.,  per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ sed. stabilisation")

manage_recruitment <- filter(treatment_df, T. < 31 & 
                                  I. <  200 & 
                                  l_max > 0.5 &
                                  phi < 400 & 
                                  model == "Stochastic recolonisation") %>%
  bind_rows(filter(treatment_df, T. < 31 & 
                     I. <  200 & 
                     l_max < 0.5 &
                     phi > 400 & 
                     model == "Stochastic dispersal")) %>% 
  select(n_sources,  model, T., per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ sed. stabilisation")

#Dispersal mgmt
#Combine: 
# Disp model with perfect parameters, max connectivity
# Recruitment limited model, basic parameters, max connectivity
HT_manage_dispersal <- filter(treatment_df, T. > 31 & 
                                  I. <  200 & 
                                  l_max > 0.5 &
                                  phi < 400 & 
                                  n_sources == 10 & 
                                  model == "Stochastic dispersal") %>%
  bind_rows(filter(treatment_df, T. > 31 & 
                     I. <  200 & 
                     l_max < 0.5 &
                     phi > 400 & 
                     n_sources == 10 & 
                     model == "Stochastic recolonisation")) %>% 
  select(n_sources,  model, T., per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ seeds/plant")

manage_dispersal <- filter(treatment_df, T. < 31 & 
                                I. <  200 & 
                                l_max > 0.5 &
                                phi < 400 & 
                                n_sources == 10 & 
                                model == "Stochastic dispersal") %>%
  bind_rows(filter(treatment_df, T. < 31 & 
                     I. <  200 & 
                     l_max < 0.5 &
                     phi > 400 & 
                     n_sources == 10 & 
                     model == "Stochastic recolonisation")) %>% 
  select(n_sources,  model, T., per5, per25, per75, per95, median_recovery) %>%
  mutate(Scenario = "+ seeds/plant")

#add different connectivity scenarios back in
#(these are all the same as low connectivity with seeding is 
# same as high connectivity with seeding)
n_sources <- unique(treatment_df$n_sources)
manage_dispersal <- manage_dispersal[rep(seq(nrow(manage_dispersal)),length(n_sources)),]
manage_dispersal$n_sources <- rep(n_sources, each = 2)

HT_manage_dispersal <- 
  HT_manage_dispersal[rep(seq(nrow(HT_manage_dispersal)),length(n_sources)),]
HT_manage_dispersal$n_sources <- rep(n_sources, each = 2)


# ------------ 
# PLOTS 
# ------------ 

# At arbitrary high and low connectivity vals plot:
#4 panels, model x HT/LT
theme_set(theme_classic())

x1 <- bind_rows(LT_LL) %>%
  bind_rows(manage_dispersal) %>%
  bind_rows(manage_recruitment) %>% 
  bind_rows(manage_light) %>% 
  bind_rows(HT_LL) %>%
  bind_rows(HT_manage_dispersal) %>%
  bind_rows(HT_manage_recruitment) %>% 
  bind_rows(HT_manage_light)%>%
  mutate(Temp = ifelse(T.<31, "Normal", "High"))%>%
  mutate(Model = case_when(
    grepl("dispersal", model) ~ "Uncertain dispersal",
    grepl("recolon", model) ~ "Uncertain recruitment"))

pd <- position_dodge(width = 0.75)

g1 <- ggplot(x1) +
  aes(x = n_sources, y = median_recovery/365, color = Scenario) +
  geom_hline(yintercept = 0) +
  geom_point(position = pd) +
  geom_linerange(aes(ymin = per5/365, ymax = per95/365),
                 position = pd) +
  facet_grid(Temp~Model) +
  # ylim(0, 10)+
  ylab("Recovery time (years)") +
  xlab("Number of source patches") + 
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c("blue", "darkgreen",
                                "purple", "black"))
  
ggsave(g1, file = "Outputs/recovery-time-mgmt-zostera.png",
       width = 6, height = 4)
# ggsave(g1, file = "Outputs/recovery-time-mgmt-halodule.png",
       # width = 6)

#Note if we wanted to plot difference between control and 
# the scenarios we need to recalculate variacnes from lambda
# then sum across those. 

#
#Above as a barplot to simplify 
#

x2 <- filter(x1, n_sources %in% c(1,8))

g2 <- ggplot(x2) +
  aes(x = factor(n_sources), y = median_recovery/365, 
       fill = Scenario) +
  geom_hline(yintercept = 0) +
  geom_col(position = "dodge2") +
  facet_grid(Temp~Model) +
  ylab("Recovery time (years)") +
  xlab("Number of source patches") +
  scale_fill_manual(values = c("lightblue", "seagreen",
                               "steelblue3", "black"))
g2  

ggsave(g2, file = "Outputs/recovery-time-mgmt-zostera-bars.png",
       width = 6, height = 3.5)
# ggsave(g2, file = "Outputs/recovery-time-mgmt-halodule-bars.png",
       # width = 6, height = 3.5)


#
#Above but as a barplot of mgmt gain
#
x_control <- filter(x2, Scenario == "No mgmt") %>%
  select(-c("model", "T.","Scenario"))

x3 <- filter(x2, Scenario != "No mgmt") %>%
  left_join(x_control, by = c("Model", "Temp",
                              "n_sources"),
            suffix = c("", "_no_mgmt")) %>%
  mutate(`Gain (years)` = (median_recovery_no_mgmt - median_recovery)/365,
         `Gain (times faster)` = (median_recovery_no_mgmt / median_recovery),
         conn = ifelse(n_sources == 1, "Low connectivity", "High connectivity")
  )

#now join two above and calculate difference

x3 %>%
  filter(Model == "Uncertain recruitment" & n_sources ==1) %>%
  select(Scenario, Temp, `Gain (years)`) %>%
  data.frame()


g3 <- ggplot(x3) +
  aes(x = Temp, y = `Gain (years)`, 
      fill = Scenario) +
  geom_hline(yintercept = 0) +
  geom_col(position = "dodge2") +
  facet_grid(conn~Model) +
  xlab("Temperature") +
  ylab("Benefit from management\n action (reduction in years)")+
  scale_fill_manual(values = c("lightblue", "navy",
                               "seagreen"))
plotly::ggplotly(g3)

ggsave(g3, file = "Outputs/mgmt-gain-zostera-bars.png",
       width = 6, height = 3.5)
# ggsave(g3, file = "Outputs/mgmt-gain-halodule-bars.png",
       # width = 6, height = 3.5)


g4 <- ggplot(x3) +
  aes(x = Temp, y = `Gain (times faster)`, 
      fill = Scenario) +
  geom_hline(yintercept = 0) +
  geom_col(position = "dodge2") +
  facet_grid(conn~Model) +
  xlab("Temperature") +
  scale_fill_manual(values = c("lightblue", "navy",
                               "seagreen"))
plotly::ggplotly(g4)

ggsave(g4, file = "Outputs/mgmt-gain-zostera-bars-multiples.png",
width = 6, height = 3.5)
# ggsave(g4, file = "Outputs/mgmt-gain-halodule-bars-multiples.png",
       # width = 6, height = 3.5)

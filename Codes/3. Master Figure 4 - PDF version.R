#Forest plot Figure
library(ggplot2)
library(ggpubr)
library(brms)
library(tidyverse)
library(tidyr)
library(ggdist)
library(tidybayes)
library(cowplot)
library(ggrepel)
library(gridExtra)

## Read in the posterior samples from the expanded analysis
post.samples_m1 <- readRDS("Output - filled/full_medium_model_full_posterior_samples_all.rds")

## data wrangle
post.samples_m1_concat <- post.samples_m1 %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                           names_to = "Predictor")

colnames(post.samples_m1_concat)[3] <- "Estimate"

## Rename Predictors
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"


## Pull draws for specific variables
post.samples_m1_concat_HWI <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Dispersal limitation")
post.samples_m1_concat_disturbance <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Historical disturbance")
post.samples_m1_concat_latitude <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Absolute latitude")
post.samples_m1_concat_BM <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Body mass")
post.samples_m1_concat_BMI <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Body mass:Dispersal limitation")
post.samples_m1_concat_Seasonality <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Seasonality")


## relevel
post.samples_m1_concat$Predictor <- factor(post.samples_m1_concat$Predictor, levels = c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality"))

## Plot
BM_medium <- ggplot(post.samples_m1_concat) +
  geom_density(aes(x = Estimate, y=..density.., fill = Predictor), lwd = 0.3, trim = TRUE, alpha = 0.3) +
  scale_fill_manual(values = c("blue", "red3", "gold", "grey", "mediumorchid", "deepskyblue")) +
  theme_classic2()  +
  ylab("Coefficient probability density") +
  xlab("Coefficient estimate") +
  ggtitle("(a)") +
  scale_y_continuous(limits = c(-1.55, 2.2), breaks = c( 0, 0.5, 1.0, 1.5), labels = scales::number_format(accuracy = 0.1), expand = c(0,-0.3)) +
  xlim(-7,7) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1.55, size = 1.7/3) +
  geom_vline(xintercept=0, lty=2, size = 0.8/3)  +
  stat_pointinterval(data=post.samples_m1_concat_HWI, 
                     aes(x=Estimate), 
                     fill = "blue", 
                     color = "blue",
                     point_fill = "white",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.15)) +
  stat_pointinterval(data=post.samples_m1_concat_disturbance, 
                     aes(x=Estimate), 
                     fill = "red", 
                     color = "red",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.3)) +
  stat_pointinterval(data=post.samples_m1_concat_latitude, 
                     aes(x=Estimate), 
                     fill = "gold", 
                     color = "gold",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.45)) +
  stat_pointinterval(data=post.samples_m1_concat_BM, 
                     aes(x=Estimate), 
                     fill = "grey", 
                     color = "grey",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.6)) +
  stat_pointinterval(data=post.samples_m1_concat_BMI, 
                     aes(x=Estimate), 
                     fill = "mediumorchid", 
                     color = "mediumorchid",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.75)) +
  stat_pointinterval(data=post.samples_m1_concat_Seasonality, 
                     aes(x=Estimate), 
                     fill = "deepskyblue", 
                     color = "deepskyblue",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.9)) +
  scale_size_continuous(range = c(5/3, 12/3)) +
  theme(legend.position = c(0.775,0.8),
        legend.key.size = unit(rel(0.9/3), "cm"),
        legend.text = element_text(size = 35/3),
        legend.title = element_blank(),
        plot.caption = element_text("Medium"),
        axis.title = element_text(size = 40/3),
        axis.text = element_text(size = 35/3),
        axis.ticks = element_line(size = 1.5/3),
        axis.ticks.length.x = unit(rel(0.5/3), "cm"),
        axis.ticks.length.y = unit(rel(0.5/3), "cm"),
        axis.line = element_line(size = 1.5/3),
        plot.title = element_text(size = 40/3)) +
  labs(title = "Expanded analysis",
       subtitle = "b") +
  theme(plot.title = element_text(size = 50/3, vjust = -4),
        plot.subtitle = element_text(size = 50/3, face = "bold", vjust = -8, hjust = 0.05))


#------------------------------------------------------------------------------------------------------------------
# Repeat with the restricted analysis
#-----------------------------------------------------------------------------------------------------------------


post.samples_h1 <- readRDS("Output - filled/full_high_model_full_posterior_samples_all.rds")

post.samples_h1_concat <- post.samples_h1 %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                           names_to = "Predictor")

colnames(post.samples_h1_concat)[3] <- "Estimate"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_h1_concat$Predictor[post.samples_h1_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"

post.samples_h1_concat_HWI <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Dispersal limitation")
post.samples_h1_concat_disturbance <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Historical disturbance")
post.samples_h1_concat_latitude <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Absolute latitude")
post.samples_h1_concat_BM <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Body mass")
post.samples_h1_concat_BMI <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Body mass:Dispersal limitation")
post.samples_h1_concat_Seasonality <- subset(post.samples_h1_concat, post.samples_h1_concat$Predictor == "Seasonality")


post.samples_h1_concat$Predictor <- factor(post.samples_h1_concat$Predictor, levels = c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality"))



BM_high <- ggplot(post.samples_h1_concat) +
  geom_density(aes(x = Estimate, y=..density.., fill = Predictor), lwd = 0.3, trim = TRUE, alpha = 0.3) +
  scale_fill_manual(values = c("blue", "red3", "gold", "grey", "mediumorchid", "deepskyblue")) +
  theme_classic2()  +
  ylab("Coefficient probability density") +
  xlab("") +
  ggtitle("(a)") +
  scale_y_continuous(limits = c(-1.55, 2.2), breaks = c( 0, 0.5, 1.0, 1.5), labels = scales::number_format(accuracy = 0.1), expand = c(0,-0.3)) +
  xlim(-7,7) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1.55, size = 1.7/3) +
  geom_vline(xintercept=0, lty=2, size = 0.8/3)  +
  stat_pointinterval(data=post.samples_h1_concat_HWI, 
                     aes(x=Estimate), 
                     fill = "blue", 
                     color = "blue",
                     point_fill = "white",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.15)) +
  stat_pointinterval(data=post.samples_h1_concat_disturbance, 
                     aes(x=Estimate), 
                     fill = "red", 
                     color = "red",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.3)) +
  stat_pointinterval(data=post.samples_h1_concat_latitude, 
                     aes(x=Estimate), 
                     fill = "gold", 
                     color = "gold",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.45)) +
  stat_pointinterval(data=post.samples_h1_concat_BM, 
                     aes(x=Estimate), 
                     fill = "grey", 
                     color = "grey",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.6)) +
  stat_pointinterval(data=post.samples_h1_concat_BMI, 
                     aes(x=Estimate), 
                     fill = "mediumorchid", 
                     color = "mediumorchid",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.75)) +
  stat_pointinterval(data=post.samples_h1_concat_Seasonality, 
                     aes(x=Estimate), 
                     fill = "deepskyblue", 
                     color = "deepskyblue",
                     .width = c(0.68, 0.95),
                     position = position_nudge(y = -0.9)) +
  scale_size_continuous(range = c(5/3, 12/3)) +
  theme(legend.position = "none", 
        legend.text = element_text(size = 35/3),
        legend.title = element_blank(),
        plot.caption = element_text("Medium"),
        axis.title = element_text(size = 40/3),
        axis.text = element_text(size = 35/3),
        axis.ticks = element_line(size = 1.5/3),
        axis.ticks.length.x = unit(rel(0.5/3), "cm"),
        axis.ticks.length.y = unit(rel(0.5/3), "cm"),
        axis.line = element_line(size = 1.5/3),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 40/3)) +
  labs(title = "Restricted analysis",
       subtitle = "a") +
  theme(plot.title = element_text(size = 50/3, vjust = -4),
        plot.subtitle = element_text(size = 50/3, face = "bold", vjust = -8, hjust = 0.05))

BM_high

png("Figures/Figure_4cut.png", width = 1572, height = 1675)

grid.arrange(BM_high, BM_medium, ncol = 1)

dev.off()



pdf("Figures/Figure_4_cut4.pdf", width = (1572/3)/72, height = (1675/3)/72)

grid.arrange(BM_high, BM_medium, ncol = 1)

dev.off()

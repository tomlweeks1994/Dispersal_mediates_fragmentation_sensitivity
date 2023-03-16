#Forest plot Figure
library(ggplot2)
library(ggpubr)
library(brms)
library(tidyverse)
library(cowplot)
library(gridExtra)

### pull the posterior samples from the sensitivity analysis
### pull restricted sample
post.samples_m1 <- readRDS("NEE main files/Second review/Data/results_data/Migration/full_high_model_full_posterior_samples_all_mig2.rds")

## conver to the longer format
post.samples_m1_concat <- post.samples_m1 %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                           names_to = "Predictor")

## rename predictors
colnames(post.samples_m1_concat)[3] <- "Estimate"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_m1_concat$Predictor[post.samples_m1_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"

## pull specific posterior distributions
post.samples_m1_concat_HWI <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Dispersal limitation")
post.samples_m1_concat_disturbance <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Historical disturbance")
post.samples_m1_concat_latitude <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Absolute latitude")
post.samples_m1_concat_BM <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Body mass")
post.samples_m1_concat_BMI <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Body mass:Dispersal limitation")
post.samples_m1_concat_Seasonality <- subset(post.samples_m1_concat, post.samples_m1_concat$Predictor == "Seasonality")


## relevel factors
post.samples_m1_concat$Predictor <- factor(post.samples_m1_concat$Predictor, levels = c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality"))
## plot 
BM_medium <- ggplot(post.samples_m1_concat) +
  geom_density(aes(x = Estimate, y=..density.., fill = Predictor), trim = TRUE, alpha = 0.3) +
  scale_fill_manual(values = c("blue", "red3", "gold", "grey", "mediumorchid", "deepskyblue")) +
  theme_classic2()  +
  ylab("Coefficient probability density") +
  xlab("Coefficient estimate") +
  #ggtitle("(a)") +
  scale_y_continuous(limits = c(-1.55, 2.2), breaks = c( 0, 0.5, 1.0, 1.5), labels = scales::number_format(accuracy = 0.1), expand = c(0,-0.3)) +
  xlim(-7,7) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1.55, size = 1.7) +
  geom_vline(xintercept=0, lty=2, size = 0.8)  +
  #annotation_custom(frigillidae_raster,xmin = -6, xmax = -4.5, ymin = 1, ymax = 1.5) +
  #annotation_custom(galliformes_raster,xmin = -5, xmax = -2.5, ymin = 1, ymax = 1.7) +
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
                     #size = 8,
                     position = position_nudge(y = -0.3)) +
  stat_pointinterval(data=post.samples_m1_concat_latitude, 
                     aes(x=Estimate), 
                     fill = "gold", 
                     color = "gold",
                     .width = c(0.68, 0.95),
                     #size = 8,
                     position = position_nudge(y = -0.45)) +
  stat_pointinterval(data=post.samples_m1_concat_BM, 
                     aes(x=Estimate), 
                     fill = "grey", 
                     color = "grey",
                     .width = c(0.68, 0.95),
                     #size = 8,
                     position = position_nudge(y = -0.6)) +
  stat_pointinterval(data=post.samples_m1_concat_BMI, 
                     aes(x=Estimate), 
                     fill = "mediumorchid", 
                     color = "mediumorchid",
                     .width = c(0.68, 0.95),
                     #size = 8,
                     position = position_nudge(y = -0.75)) +
  stat_pointinterval(data=post.samples_m1_concat_Seasonality, 
                     aes(x=Estimate), 
                     fill = "deepskyblue", 
                     color = "deepskyblue",
                     .width = c(0.68, 0.95),
                     #size = 8,
                     position = position_nudge(y = -0.9)) +
  scale_size_continuous(range = c(5, 12)) +
  theme(legend.position = c(0.775,0.8),
        legend.key.size = unit(rel(0.9), "cm"),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        plot.caption = element_text("Medium"),
        #axis.line.x = element_blank(),
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.line = element_line(size = 1.5),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text(size = 40)) +
  labs(title = "No long-distance migrants") +
  theme(plot.title = element_text(size = 50, vjust = -6, hjust = 0.05))


## plot PNG
png("Figures/Publication/Main/Figure_4_no_mig.png", width = 1500, height = 850)

grid.arrange(BM_medium, ncol = 1)

dev.off()




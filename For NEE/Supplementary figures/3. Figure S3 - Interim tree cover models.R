#Forest plot Figure
library(ggplot2)
library(ggpubr)
library(brms)
library(dplyr)
library(tidyr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(png)
library(tidyr)
library(dplyr)

high_large_change <- readRDS("NEE main files/Second Review/Data/results_data/Largest_change/OG/full_high_model_full_posterior_samples_largest_change.rds")
high_large_change_sum <- readRDS("NEE main files/Second Review/Data/results_data/Largest_change/OG/full_high_model_full_summary_largest_change.rds")


post.samples_hl_concat <- high_large_change %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                           names_to = "Predictor") %>%
  group_by(Predictor) %>%
  dplyr::summarise(Estimate = mean(value),
                   lower = quantile(value, 0.025),
                   upper = quantile(value, 0.975)
                  )

post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_hl_concat$Predictor[post.samples_hl_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"


post.samples_hl_concat$Predictor <- factor(post.samples_hl_concat$Predictor, levels = rev(c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality")))

post.samples_hl_concat$change <- "Large"


high_small_change <- readRDS("NEE main files/Second Review/Data/results_data/Smallest_change/OG/full_high_model_full_posterior_samples_smallest_change.rds")
post.samples_hs_concat <- high_small_change %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                               names_to = "Predictor") %>%
  group_by(Predictor) %>%
  dplyr::summarise(Estimate = mean(value),
                   lower = quantile(value, 0.025),
                   upper = quantile(value, 0.975)
  )

post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_hs_concat$Predictor[post.samples_hs_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"


post.samples_hs_concat$Predictor <- factor(post.samples_hs_concat$Predictor, levels = rev(c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality")))

post.samples_hs_concat$change <- "Small"


full_high <- rbind(post.samples_hs_concat, post.samples_hl_concat)

BM_high <- ggplot(full_high) +
  geom_point(aes(y = Estimate, x = Predictor, col = change), position = position_dodge(0.5), size = 5) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, col = change), size = 1.1, width = 0.4, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Red", "deepskyblue"), name = "Interim tree cover change:") +
  scale_y_continuous(breaks=c(-10, -5, 0, 5, 10)) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  #ylim(-13.85,13.85) +
  ylab("Coefficient Estimate") +
  coord_flip() +
  theme_classic() +
  ggtitle("a") +
  theme(legend.position = "top",
        #legend.key.size = unit(rel(0.9), "cm"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        plot.title = element_text(size = 40, face = "bold", vjust = -12, hjust = 0.075),
        axis.title = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 30),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.line = element_line(size = 1.5))


################################################################################

medium_large_change <- readRDS("NEE main files/Second Review/Data/results_data/Largest_change/OG/full_medium_model_full_posterior_samples_largest_change.rds")


post.samples_ml_concat <- medium_large_change %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                               names_to = "Predictor") %>%
  group_by(Predictor) %>%
  dplyr::summarise(Estimate = mean(value),
                   lower = quantile(value, 0.025),
                   upper = quantile(value, 0.975)
  )

post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_ml_concat$Predictor[post.samples_ml_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"


post.samples_ml_concat$Predictor <- factor(post.samples_ml_concat$Predictor, levels = rev(c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality")))

post.samples_ml_concat$change <- "Large"


medium_small_change <- readRDS("NEE main files/Second Review/Data/results_data/Smallest_change/OG/full_medium_model_full_posterior_samples_smallest_change.rds")
post.samples_ms_concat <- medium_small_change %>% pivot_longer(cols = c("b_nHWI_z", "b_abs_lat_z", "b_DisturbanceHigh", "b_Mass_log_z", "b_nHWI_z:Mass_log_z", "b_Seasonality_z"),
                                                               names_to = "Predictor") %>%
  group_by(Predictor) %>%
  dplyr::summarise(Estimate = mean(value),
                   lower = quantile(value, 0.025),
                   upper = quantile(value, 0.975)
  )

post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_nHWI_z"] <- "Dispersal limitation"
post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_DisturbanceHigh"] <- "Historical disturbance"
post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_abs_lat_z"] <- "Absolute latitude"
post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_Mass_log_z"] <- "Body mass"
post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_nHWI_z:Mass_log_z"] <- "Body mass:Dispersal limitation"
post.samples_ms_concat$Predictor[post.samples_ms_concat$Predictor == "b_Seasonality_z"] <- "Seasonality"


post.samples_ms_concat$Predictor <- factor(post.samples_ms_concat$Predictor, levels = rev(c("Dispersal limitation", 
                                                                                        "Historical disturbance",
                                                                                        "Absolute latitude",
                                                                                        "Body mass",
                                                                                        "Body mass:Dispersal limitation",
                                                                                        "Seasonality")))

post.samples_ms_concat$change <- "Small"


full <- rbind(post.samples_ms_concat, post.samples_ml_concat)

BM_medium <- ggplot(full) +
  geom_point(aes(y = Estimate, x = Predictor, col = change), position = position_dodge(0.5), size = 5) +
  geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper, col = change), size = 1.1, width = 0.4, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Red", "deepskyblue"), name = "Interim tree cover change:") +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  ylab("Coefficient Estimate") +
  coord_flip() +
  ggtitle("b") +
  theme_classic() +
  theme(legend.position = "none",
        #legend.key.size = unit(rel(0.9), "cm"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        plot.title = element_text(size = 40, face = "bold", vjust = -12, hjust = 0.075),
        axis.title = element_text(size = 30),
        #axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.line = element_line(size = 1.5))
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),

library(patchwork)

plot <- (BM_high | BM_medium) + plot_layout(guides = "collect") & theme(legend.position = "top")

png("Figures/Publication/Extended Data/Figure S3 interim tree cover change.png",
    width = 1450, height = 750)
plot
dev.off()

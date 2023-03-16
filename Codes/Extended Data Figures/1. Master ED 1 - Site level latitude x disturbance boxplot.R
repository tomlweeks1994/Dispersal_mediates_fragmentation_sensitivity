library(tidyverse)
library(boot)
library(broom)
library(cowplot)
library(dplyr)
library(ggplot2)
library(readxl)

## Read in Data
site_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Seasonality")

species_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Order", "nHWI", "Hand.Wing.Index", "Mass (g)", "Forest_dependency")

site_sp_list <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 4) %>%
  filter(Species_name %in% species_data$Species_name)

##Combine
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## Add Suffic
data$Disturbance <- paste(data$Disturbance, "disturbance")

## Convert Latitude to absolute latitude for modelling
data$abs_lat <- abs(data$Latitude)

## Drop to just the study_landscapes
data <- data %>%
  distinct(Dataset_ID, .keep_all = TRUE)


## Subset the different disturbance groups
dist_1 <- data %>% dplyr::filter(Disturbance == "High disturbance")
dist_0 <- data %>% dplyr::filter(Disturbance == "Low disturbance")

## run the wilcoxon test
test <- wilcox.test(dist_1$abs_lat, dist_0$abs_lat, alternative = "two.sided")

## save plot
dev.size("px")
jpeg("Figures/Extended data figure 1.jpg",
    width = 825/1.5,
    height = 691/1.5)

ggplot(data) + 
  xlab("") +
  ylab("Absolute latitude") +
  geom_boxplot(aes(x = Disturbance, y = abs_lat, fill = Disturbance), width = 0.4, draw_quantiles = 0.5, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("red", "deepskyblue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15))  +
  annotate("text", x= 2.33, y= 0, label= paste("W-statistic:", test$statistic), size = 5) +
  annotate("text", x= 2.34, y= - 3.3, label= paste("P-value:", round(test$p.value, 3)), size = 5)

dev.off()

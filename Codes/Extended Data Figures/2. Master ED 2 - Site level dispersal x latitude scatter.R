library(tidyverse)
library(boot)
library(broom)
library(cowplot)
library(dplyr)
library(ggplot2)
library(readxl)


## Read data
site_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Seasonality")

species_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Order", "nHWI", "Hand.Wing.Index", "Mass (g)", "Forest_dependency")

site_sp_list <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 4) %>%
  filter(Species_name %in% species_data$Species_name)

## Join datasets
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## Create response variables
## Restricted
data$dep_sense <- 0
data$dep_sense[data$Forest_dependency %in% c("High") & data$Classification == "Forest core"] <- 1

## Expanded
data$dep_sense_medium <- 0
data$dep_sense_medium[data$Forest_dependency %in% c("High", "Medium") & data$Classification == "Forest core"] <- 1

## take negative logs
data$log_nHWI <- -(log(data$Hand.Wing.Index))

## Summarise dataset to the tings we want
data1 <- data %>%
  dplyr::group_by(Dataset_ID) %>%
  dplyr::summarise(mean_HWI = mean(log_nHWI),
            abs_lat = abs(Latitude),
            n = n(),
            cores = length(which(Classification == "Forest core")),
            sensitive = length(which(dep_sense == 1))) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## Calculate proportions
data1$prop_core <- data1$cores/data1$n
data1$prop_sensitive <- data1$sensitive/data1$n

## reformat
data1$Disturbed <- as.factor(data$Disturbance[match(data1$Dataset_ID, data$Dataset_ID)])

## run the analysis mode
mod <- lm(mean_HWI ~ abs_lat, data = data1)

## create newdata for predictions and predict
newdata <- data.frame(abs_lat = seq(min(data1$abs_lat), max(data1$abs_lat), length.out = 100))
predicted <- predict(mod, newdata = newdata, se.fit = TRUE)
predicted$abs_lat <- newdata$abs_lat
predicted$lower <- predicted$fit - predicted$se.fit
predicted$upper <- predicted$fit + predicted$se.fit

## make into a df
predicted <- as.data.frame(predicted)

##ggplot
plot <- ggplot() +
  geom_line(data=predicted, aes(x = abs_lat, y= fit), color="purple4", size=1.5) +
  geom_ribbon(data=predicted, aes(x = abs_lat, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE, alpha=0.2) +
  geom_point(data = data1, aes(x = abs_lat, y = mean_HWI, colour = Disturbed), size = 5) +
  ylab("Mean assemblage dispersal limitation (nHWI)") +
  xlab("Absolute latitude") +
  #ggtitle("(a)") +
  scale_colour_manual(name = "Disturbance:",
                        values = c("red", "deepskyblue")) +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 27),
        legend.text = element_text(size = 27),
        axis.title.y = element_text(size = 25, vjust = 2.5),
        axis.title.x = element_text(size = 25, vjust = -0.5),
        axis.text = element_text(size = 20))

## add summary stats
pmain <-  plot  +
  annotate("text", x= 13.4, y= -3.39, label= paste0("Coefficient: ", round(as.numeric(mod$coefficients[2]), 3)), size = 7.5) +
  annotate("text", x= 15, y= - 3.422, label= paste0("P-value < 0.001"), size = 7.5)

## save plot
jpeg("Figures/Extended data figure 2.jpg",
    width = 823/1.25, height = 779/1.25)

plot(pmain)

dev.off()








library(tidyverse)
library(boot)
library(broom)
library(dplyr)
library(readxl)

#Read in the data
site_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "nHWI", "Mass (g)",
                "Forest_dependency", "Disturbance_range_breeding",
                "Seasonality_range_breeding", "TrophicNiche")

site_sp_list <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)

##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## take averages of the study_landscapes
data1 <- data %>%
  group_by(Dataset_ID) %>%
  summarise(mean_dist = mean(Disturbance_range_breeding),
            mean_local_seasonality = mean(Seasonality_range_breeding),
            mean_abs_lat = abs(mean(Latitude)),
            n = n()) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## scale
data1$mean_local_seasonality <- scale(data1$mean_local_seasonality)

## model
mod <- lm(mean_local_seasonality ~ mean_dist, data = data1)

##predict model
newdata <- data.frame(mean_dist = seq(min(data1$mean_dist), max(data1$mean_dist), length.out = 100))
predicted <- predict(mod, newdata = newdata, se.fit = TRUE)
predicted$mean_dist <- newdata$mean_dist
predicted$lower <- predicted$fit - predicted$se.fit
predicted$upper <- predicted$fit + predicted$se.fit

predicted <- as.data.frame(predicted)

## ggplot
plot <- ggplot() +
  geom_line(data=predicted, aes(x = mean_dist, y= fit), color="purple4", size=1.5) +
  geom_ribbon(data=predicted, aes(x = mean_dist, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE, alpha=0.2) +
  geom_point(data = data1, aes(x = mean_dist, y = mean_local_seasonality), size = 5, alpha = 0.5) +
  ylab("Mean seasonality (z-score)") +
  xlab("Mean disturbance") +
  #ggtitle("(a)") +
  scale_colour_discrete(name = "Disturbance:") +
  theme_bw() +
  ggtitle("b") +
  theme(legend.position = "top",
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.title.y = element_text(size = 25),# vjust = 2.5),
        axis.title.x = element_text(size = 25, vjust = -0.5),
        axis.text = element_text(size = 20),        
        plot.title = element_text(size = 40, face = "bold", hjust = 0.08, vjust = -9)) 

pmain <-  plot  +
  annotate("text", x= 0.7, y= -1.15, label= paste0("Coefficient: ", round(as.numeric(mod$coefficients[2]), 3)), size =7) +
  annotate("text", x= 0.703, y= -1.3, label= paste0("P-value = 0.0004"), size = 7)

grid.arrange(whisker,pmain, nrow = 1, ncol = 2)


##############t#################################################################
## reread the data as
################################################################################
#Read in the data
site_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
dplyr::select("Species_name", "Family", "Order", "nHWI", "Mass (g)",
              "Forest_dependency", "Disturbance_range_breeding",
              "Seasonality_range_breeding", "TrophicNiche")

site_sp_list <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)

##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## distinct data for each study_landscape
data <- data %>%
  distinct(Dataset_ID, .keep_all = TRUE)

## scale
data$Seasonality <- scale(data$Seasonality)

## subset data by disturbances
dist_1 <- data %>% dplyr::filter(Disturbance == "High")
dist_0 <- data %>% dplyr::filter(Disturbance == "Low")

## run a wilcox test
test <- wilcox.test(dist_1$Seasonality, dist_0$Seasonality, alternative = "two.sided")

## add a suffix
data$Disturbance <- paste(data$Disturbance, "disturbance")

## relevel
data$Disturbance <- factor(data$Disturbance, levels = c("High disturbance", 
                                                           "Low disturbance"))
## ggplot
whisker <- ggplot(data) + 
  xlab(" ") +
  ylab("Local seasonality (z-score)") +
  geom_boxplot(aes(x = Disturbance, y = Seasonality, fill = Disturbance), width = 0.4, draw_quantiles = 0.5, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("red", "deepskyblue")) +
  theme_bw() +
  ggtitle("a") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 25),# vjust = 2.5),
        axis.title.x = element_text(size = 25, vjust = -0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.margin = margin(0.2,1,0.2,0, unit = "cm"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.08, vjust = -9))  +
  annotate("text", x= 2.28, y= -1.65, label= paste("W-statistic:", test$statistic), size = 7) +
  annotate("text", x= 2.3, y= -1.9, label= paste("P-value:", round(test$p.value, 3)), size = 7)


#[1] create the dataframe
jpeg(filename = "NEE main files/Figures/Extended Data/ED_figure_9.jpg", width = 1198, height = 552)
grid.arrange(whisker,pmain, nrow = 1, ncol = 2)
dev.off()



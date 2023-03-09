library(EnvStats)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

#data <- read.csv("Analysis/HPC folder/data2009_medium.csv")
#Read in the data
site_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "nHWI", "Mass (g)", "Forest_dependency",
                "TrophicNiche")

site_sp_list <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)



##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## create HWI & abs_lat
data$abs_lat <- abs(data$Latitude)
data$HWI <- -data$nHWI

########## Climate zone split based on landscape Lat #######
data_temp <- data %>% filter(abs_lat > 20)
data_tropics <- data %>% filter(abs_lat <= 20)

## reduce tp required columns
nHWI_data_temp <- data_temp %>% dplyr::select(Species_name, HWI, nHWI, abs_lat, Dataset_ID, TrophicNiche) %>% distinct(Species_name, .keep_all = TRUE)

## only pull diet niches > 5 species in
toPlot <- names(which(table(nHWI_data_temp$TrophicNiche) >= 5 ))

# reduce to create df with only the correct Trop niches
nHWI_data_temp <- nHWI_data_temp %>% dplyr::filter(TrophicNiche %in% toPlot)

## get median nHWIs  for coloring plor
nHWI_data_temp <- nHWI_data_temp %>% group_by(TrophicNiche) %>%
  mutate(nHWI_median = median(nHWI))

## 
nHWI_data_temp$climatic_zone <- "a)"

nHWI_data_temp$TrophicNiche <- factor(nHWI_data_temp$TrophicNiche, levels = c("Invertivore",
                                                                              "Frugivore",
                                                                              "Omnivore",
                                                                              "Granivore",
                                                                              "Vertivore",
                                                                              "Nectarivore"))


## reduce tp required columns
nHWI_data_tropics <- data_tropics %>% dplyr::select(Species_name, HWI, nHWI, abs_lat, Dataset_ID, TrophicNiche) %>% distinct(Species_name, .keep_all = TRUE)

## only pull diet niches > 5 species in
toPlot <- names(which(table(nHWI_data_tropics$TrophicNiche) > 5 ))

# reduce to create df with only the correct Trop niches
nHWI_data_tropics <- nHWI_data_tropics %>% dplyr::filter(TrophicNiche %in% toPlot)

## get median nHWIs  for coloring plor
nHWI_data_tropics <- nHWI_data_tropics %>% group_by(TrophicNiche) %>%
  mutate(nHWI_median = median(nHWI))

##
nHWI_data_tropics$climatic_zone <- "b)"

nHWI_data_tropics$TrophicNiche <- factor(nHWI_data_tropics$TrophicNiche, levels = c("Invertivore",
                                                                                    "Frugivore",
                                                                                    "Omnivore",
                                                                                    "Granivore",
                                                                                    "Vertivore",
                                                                                    "Nectarivore"))

## bind the data
nHWI_data_all <- rbind(nHWI_data_temp, nHWI_data_tropics)

## run anova models to pull out F-stat
trop_mod <- lm(nHWI ~ TrophicNiche, data = nHWI_data_tropics)
temp_mod <- lm(nHWI ~ TrophicNiche, data = nHWI_data_temp)
f_stat_trop <- anova(trop_mod)$`F value`[1]
f_stat_temp <- anova(temp_mod)$`F value`[1]


## plot
temp_a <- nHWI_data_temp %>%
  ggplot(aes(x = reorder(TrophicNiche, HWI, na.rm = TRUE), y= nHWI, fill = nHWI_median)) +
  geom_boxplot() +
  scale_fill_gradientn(colors = c("white","#ffbcb5", "red3"),
                       values = c(0, 0.8, 1)) +
  stat_n_text(size = 10, vjust = 1.5) +
  labs(y = "Dispersal limitation (nHWI)", x = "Trophic Niche") +
  theme_classic() +
  ggtitle("a Temperate forest birds") +
  ylim(-80,0) +
  #facet_wrap(~climatic_zone, nrow = 2) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 25, face = "bold", vjust = -6, hjust = 0.05),
        axis.text.x = element_blank(),#text(size = 20, angle = 45, hjust = 1.1),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        # strip.background = element_blank(),
        # strip.text.x = element_text(size = 25, angle = 0, hjust = 0, vjust = 1),
        # panel.spacing = unit(2, "lines"),
        plot.margin = margin(0,0,1,0, "cm"),
        legend.position = "none") +
  annotate("text", x=5.55, y=-0.038, label= paste0("F-statistic: ", round(as.numeric(f_stat_temp), 2)), size = 9) +
  annotate("text", x=5.575, y=-6, label= paste0("P-value < 0.001"), size = 9)


#plot
trop_b <- nHWI_data_tropics %>%
  ggplot(aes(x = reorder(TrophicNiche, HWI, na.rm = TRUE), y= nHWI, fill = nHWI_median)) +
  geom_boxplot() +
  scale_fill_gradientn(colors = c("white","#ffbcb5", "red3"),
                       values = c(0, 0.8, 1)) +
  stat_n_text(size = 10) + #, vjust = -16) +
  labs(y = "Dispersal limitation (nHWI)", x = "Trophic Niche") +
  theme_classic() +
  ggtitle("b Tropical forest birds") +
  ylim(-80,0) +
  #facet_wrap(~climatic_zone, nrow = 2) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 25, face = "bold", vjust = -6, hjust = 0.05),
        axis.text.x = element_text(size = 25, angle = 45, hjust = 1.1),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25),
        # strip.background = element_blank(),
        # strip.text.x = element_text(size = 25, angle = 0, hjust = 0, vjust = 1),
        # panel.spacing = unit(2, "lines"),
        plot.margin = margin(1,0,0,0, "cm"),
        legend.position = "none") +
  annotate("text", x=5.55, y=-0, label= paste0("F-statistic: ", round(as.numeric(f_stat_trop), 2)), size = 9) +
  annotate("text", x=5.625, y=-6, label= paste0("P-value < 0.001"), size = 9)


## PNG ##
png("Figures/Publication/Extended Data/ED 8 - Diet x nHWI_.png",
    width = 843, height = 1124)

gridExtra::grid.arrange(temp_a, trop_b, ncol = 1, heights = c(4,5))

dev.off()
dev.size("px")
# [1]  843 1124




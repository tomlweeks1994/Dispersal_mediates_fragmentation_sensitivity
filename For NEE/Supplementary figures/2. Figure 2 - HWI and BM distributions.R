library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

#phylo_data <- read.csv("Analysis/HPC folder/phylo_data2009_medium.csv")
#Read in the data
site_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "nHWI", "Mass (g)", "Forest_dependency",
                "Disturbance_range_total", "Disturbance_range_breeding",
                "Centroid_latitude_range_total (Behrmann's projection)",
                "Centroid_latitude_range_breeding (Behrmann's projection)",
                "Seasonality_range_total", "Seasonality_range_breeding")

site_sp_list <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)


##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## create HWI
data$HWI <- -data$nHWI

## read in global sample
global_sample <- read.csv("../../Useful materials/AVONET1 - Jetz.csv") %>% dplyr::select(ï..Species3, Hand.Wing.Index, Mass)
colnames(global_sample) <- c("Species_name", "HWI", "Mass")
global_sample$Mass <- log(global_sample$Mass)


data$Passerine <- "Nonpasserine"
data$Passerine[data$Order == "Passeriformes"] <- "Passerine"

data$dep_sense <- "Insensitive"
data$dep_sense[data$Forest_dependency %in% c("High") & data$Classification == "Forest Core"] <- "Sensitive"
data$dep_sense_medium <- "Insensitive"
data$dep_sense_medium[data$Forest_dependency %in% c("High", "Medium") & data$Classification == "Forest Core"] <- "Sensitive"

data_sample <- data %>% dplyr::select(Species_name, HWI, `Mass (g)`, Passerine, Forest_dependency, dep_sense, dep_sense_medium)
colnames(data_sample) <- c("Species_name", "HWI", "Mass", "Passerine", "Forest_dep", "dep_sense", "dep_sense_medium")

data_sample$Mass <- log(data_sample$Mass)
data_sample$Forest_dep <- factor(data_sample$Forest_dep, levels = c("Low", "Medium", "High"))

global_sample$sample <- "Global"
data_sample$sample <- "Study"
global_vs_sample <- rbind(global_sample, data_sample %>% 
                            dplyr::select(Species_name, HWI, Mass, sample)) %>% 
  pivot_longer(cols = c("HWI", "Mass"), names_to = "Variable")


sample_plot <- ggplot(global_vs_sample, aes(x = sample, y = value, fill = sample)) +
    geom_violin(draw_quantiles = c(0.5),
                alpha = 0.5,
                lwd = 1) +
    ggtitle("a") +
    scale_fill_manual(name = "Sample:         ", values = c(Global = "Red", Study = "#ff928a")) +
    facet_wrap(~Variable, scales = "free") +
    theme_bw()+
    theme(legend.position = "right",# c(1.385, 0.5),
          legend.title = element_text(size = 30),
          strip.text = element_text(size = 25, face = "bold"),
          legend.text = element_text(size = 25),
          plot.title = element_text(size = 35, face = "bold",  hjust = -0.1, vjust = -3.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          plot.margin = unit(c(0,0,0.01,0.05), "npc"),
          axis.title.y = element_blank(),
          axis.text.x = element_blank())


data_sample_piv <- data_sample %>% 
  pivot_longer(cols = c("HWI", "Mass"), names_to = "Variable")

data_sample_piv$Forest_dep <- factor(data_sample_piv$Forest_dep, levels = c("High", "Medium", "Low"))

dep_plot <-ggplot(data_sample_piv, aes(x = Forest_dep, y = value, fill = Forest_dep)) +
  geom_violin(draw_quantiles = c(0.5),
              alpha = 0.5,
              lwd = 1) +
    scale_fill_manual(name = "Forest \ndependency:   ", values = c( High ="blue",  Medium = "#645ff5", Low = "#b5b3f5")) +
    ggtitle("b") +
  theme_bw()+
  facet_wrap(~Variable, scales = "free") +
  theme(legend.position = "right",# c(1.385, 0.5),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        strip.text = element_text(size = 25, face = "bold"),
        plot.title = element_text(size = 35, face = "bold",  hjust = -0.1, vjust = -3.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25),
        plot.margin = unit(c(0,0,0.01,0.05), "npc"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank())


pas_plot <- ggplot(data_sample_piv, aes(x = Passerine, y = value, fill = Passerine)) +
  geom_violin(draw_quantiles = c(0.5),
              alpha = 0.5,
              lwd = 1) +
  scale_fill_manual(name = "Clade:", values = c(Nonpasserine = "gold", Passerine = "#fcf992")) +
  ggtitle("c") +
  theme_bw()+
  facet_wrap(~Variable, scales = "free") +
  theme(legend.position = "right",# c(1.385, 0.5),
        legend.title = element_text(size = 30),
        strip.text = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 35, face = "bold", hjust = -0.1, vjust = -3.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25),
        plot.margin = unit(c(0,0,0.01,0.05), "npc"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank())


restricted <- ggplot(data_sample_piv, aes(x = dep_sense_medium, y = value, fill = dep_sense_medium)) +
  geom_violin(draw_quantiles = c(0.5),
              alpha = 0.5,
              lwd = 1) +
    scale_fill_manual(name = paste("Fragmentation\nsensitivity:"), values = c(Insensitive = "gold", Sensitive = "#fcf992")) +
    ggtitle("c") +
    theme_bw()+
  facet_wrap(~Variable, scales = "free") +
  theme(legend.position = "right",# c(1.385, 0.5),
        legend.title = element_text(size = 30),
        strip.text = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 35, face = "bold", hjust = -0.1, vjust = -3.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25),
        plot.margin = unit(c(0,0,0.01,0.05), "npc"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank())


#ggarrange(sample_plot, dep_plot, pas_plot, restricted, nrow = 4, ncol = 1)#, heights = c(4,5))

# > dev.size("px")
# [1] 1114 1358

library(patchwork)

plots <- sample_plot / dep_plot / restricted

png("Figures/Publication/Extended Data/ED 1 - HWI + BM densities.png", width = 820, height = 800)

cowplot::plot_grid(plotlist = list(sample_plot, dep_plot, restricted), ncol = 1,
                   align = 'v')
dev.off()


tiff("Figures/Publication/Extended Data/ED 1 - HWI + BM densities_tf.tif", width = 820, height = 800)

cowplot::plot_grid(plotlist = list(sample_plot, dep_plot, restricted), ncol = 1,
                   align = 'v')
dev.off()










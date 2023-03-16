library(raster)
library(ape)
library(geiger)
library(phytools)
library(ggtree)
library(ggnewscale)
library(ggplot2)
library(gridExtra)
library(readxl)
library(dplyr)
library(tidyr)

#phylo_data <- read.csv("Analysis/HPC folder/phylo_data2009_medium.csv")
#Read in the data
site_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Longitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 3) %>%
  dplyr::select("Species_name","Family", "Order", "nHWI", "Mass (g)", "Forest_dependency", "Migration",
                "Disturbance_range_total", "Disturbance_range_breeding",
                "Centroid_latitude_range_total (Behrmann's projection)",
                "Centroid_latitude_range_breeding (Behrmann's projection)",
                "Seasonality_range_total", "Seasonality_range_breeding")

site_sp_list <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 4)

##Join them up
phylo_data <- left_join(site_sp_list, site_data)
phylo_data <- left_join(phylo_data, species_data)

## Restricted analysis ##
phylo_data$dep_sense <- 0
phylo_data$dep_sense[phylo_data$Forest_dependency == "High" & phylo_data$Classification == "Forest Core"] <- 1

## Expanded analysis ##
phylo_data$dep_sense_medium <- 0
phylo_data$dep_sense_medium[phylo_data$Forest_dependency %in% c("High", "Medium") & phylo_data$Classification == "Forest Core"] <- 1


#-----------------------------------
#panel a) disturbance world map
#-----------------------------------
## read in forest loss raster
forest_loss <- raster("Rasters/forest_loss.tif")
forest_loss[forest_loss > 0] <- 5

## read in natural disturbance raster
natural_raster <- raster("Rasters/Natural_raster.grd")
natural_raster[which(getValues(natural_raster) == 0)] <- NA
natural_raster[which(getValues(natural_raster) > 0)] <- 1

## load environmental variables (will be used as a mask)
envar <- raster::getData("worldclim", var = "bio", res=2.5)

## pull first raster
worldmap_blank <- envar[[1]]

## set all NAs to 10
worldmap_blank[which(getValues(worldmap_blank) != "NA")] <- 10

## stack rasters
raster_stack <- natural_raster
raster_stack <- addLayer(raster_stack, forest_loss)
raster_stack <- addLayer(raster_stack, worldmap_blank)

## create disturbance raster by summing the rastr values
disturbance_raster <- stackApply(raster_stack, indices = rep(1, nlayers(raster_stack)), fun = sum, na.rm = T)

## create df
rasdf <- as.data.frame(disturbance_raster, xy=TRUE)

## remove 10 to pull all old NAs back to 0
rasdf$index_1 <- rasdf$index_1 - 10

## anything bellow 0 (old 0s/No data) turn to NA
rasdf$index_1[rasdf$index_1 < 0] <- NA

## remove these values and format df
rasdf <- rasdf %>% drop_na()
colnames(rasdf)[3] <- "Disturbance_score"
rasdf$Disturbance <- NA

## create character data for the values
rasdf$Disturbance[which(rasdf$Disturbance_score == 1)] <- "Natural"
rasdf$Disturbance[which(rasdf$Disturbance_score == 5)]  <- "Anthropogenic (forest loss)"
rasdf$Disturbance[which(rasdf$Disturbance_score == 6)] <- "Natural & Anthropogenic"
rasdf$Disturbance[which(rasdf$Disturbance_score == 0)] <- "Minimal disturbance"

## relevel
rasdf$Disturbance <- factor(rasdf$Disturbance, levels = c("Natural & Anthropogenic",
                                                          "Natural",
                                                          "Anthropogenic (forest loss)",
                                                          "Minimal disturbance"))

## set color scheme
myColors <- c("firebrick", "red", "tomato", "skyblue")
names(myColors) <- levels(rasdf$Disturbance)
colScale <- scale_fill_manual(name = paste0("Historical",'\n',"disturbance"), 
                              values = myColors,
                              na.value = "grey80")

##ggplot
a <- ggplot() +
  geom_tile(aes(x=x, y=y, fill= factor(Disturbance)), data=rasdf) +
  colScale +
  geom_point(aes(x= Longitude, y=Latitude), colour = "black", size = 8, data = phylo_data) +
  geom_point(aes(x= Longitude, y=Latitude), colour = "yellow", size = 6, data = phylo_data) +
  theme_void() +
  ggtitle("a") +
  theme(legend.position= c(0.1,0.4),
        plot.margin = unit(c(0,0,0,6.5), "cm"),
        legend.direction='vertical',
        legend.box = "vertical",
        plot.title = element_text(size = 80, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 50)) +
  ylab("Latitude") + 
  xlab("Longitude")


#--------------------------------------------------
#Panel c) Disturbance vs fragmentation sensitivity
#--------------------------------------------------
## create fake data for my correlation concept graph
dat <- data.frame(x = 50:5, y = 5:50)

## plot simple plot of disturbance vs frag sensitiivity
plot_disturbance <- ggplot(dat, aes(x, y)) + 
  geom_line(col="skyblue", size = 3) + 
  xlab(paste("Fragmentation",'\n',"sensitivity")) + 
  ylab(paste("Historical disturbance")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ggtitle("c") +
  coord_flip() +
  theme_classic() + 
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=45),
        axis.title.y = element_text(vjust = 3.5,  face = "bold"),
        axis.title.x = element_text(vjust = -1.5,  face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 80, face = "bold", hjust = -0.3, vjust = 6),
        plot.margin = unit(c(3,3,3,3), "cm")) +
  geom_segment(aes(x=1, xend = 55 , y=1, yend = 1), size= 4,
               arrow = arrow(length = unit(1,"cm"))) +
  geom_segment(aes(x=1, xend = 1, y=1, yend = 55), size= 4,
               arrow = arrow(length = unit(1,"cm")))

## have to remove boarders to allow plotting correctly
c <- ggplot_gtable(ggplot_build(plot_disturbance))
c$layout$clip[c$layout$name=="panel"] <- "off"



#-----------------------------------------------------
#Panel b) HWI distribution map
#-----------------------------------------------------
## read in my HWI raster
all_HWI <- raster("Rasters/All_species_HWI_SR5.grd")

## create df and drop NAs
all_rasdf <- as.data.frame(all_HWI, xy = TRUE) %>% drop_na()

## set all 0s(no data areas) to NA - this will be given a color in the plot
all_rasdf$layer[all_rasdf$layer == 0] <- NA

## create a df and the dispersal limitation values
colnames(all_rasdf)[3] <- "HWI"
all_rasdf$inv_HWI <- 1/all_rasdf$HWI 
all_rasdf$neg_log_HWI <- -log(all_rasdf$HWI) 
all_rasdf$neg_HWI <- -(all_rasdf$HWI)

## ggplot
b <- ggplot() +
  geom_tile(aes(x=x, y=y, fill= neg_log_HWI), data=all_rasdf) +
  scale_fill_gradientn(name = paste0("Dispersal", '\n', "limitation", '\n', "(nHWI)"),
                       colors = c("dodgerblue4", "dodgerblue3", "dodgerblue2","dodgerblue1", "skyblue", "tomato", "red1","firebrick"),
                       na.value = "grey80",
                       breaks = c(-3.2, -3.6, -4.0)) +
  
  ylab("Latitude") + 
  xlab("Longitude") +
  theme_void() +
  ggtitle("b") +
  theme(legend.position = c(0.045,0.40),
        plot.margin = unit(c(0,0,0,6.5), "cm"),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.key.height = unit(0.035, "npc"),
        legend.text = element_text(size = 50),
        plot.title = element_text(size = 80, face = "bold"),
        legend.title = element_text(size = 60, face = "bold")) +
  geom_point(aes(x= Longitude, y=Latitude), colour = "black", size = 8, data = phylo_data) +
  geom_point(aes(x= Longitude, y=Latitude), colour = "yellow", size = 6, data = phylo_data)


#--------------------------------------------------
#Panel d) Dispersal vs fragmentation sensitivity
#--------------------------------------------------
## create fake data for correlation concept graph
dat <- data.frame(x = 5:50, y = 5:50)

## ggplot (simple plot for dispersal vs frag sensitivity)
plot_HWI <- ggplot(dat, aes(x, y)) + 
  geom_line(col="skyblue", size = 3) + 
  xlab(paste("Fragmentation",'\n', "sensitivity")) + 
  ylab(paste("Dispersal limitation (nHWI)")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() + 
  coord_flip() +
  ggtitle("d") +
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=45),
        axis.title.y = element_text(vjust = 3.5,  face = "bold"),
        axis.title.x = element_text(vjust = -1.5,  face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size = 80, face = "bold", hjust = -0.3, vjust = 6),
        plot.margin = unit(c(3,3,3,3), "cm")) +
  geom_segment(aes(x=1, xend = 55 , y=1, yend = 1), size=4,
               arrow = arrow(length = unit(1,"cm"))) +
  geom_segment(aes(x=1, xend = 1, y=1, yend = 55), size=4,
               arrow = arrow(length = unit(1,"cm")))

## have to remove boarders to allow plotting correctly
d <- ggplot_gtable(ggplot_build(plot_HWI))
d$layout$clip[d$layout$name=="panel"] <- "off"

### SAVE TIFFS AND PNGS ###################
## TIFF
tiff("Figures/Figure_2.tif",
     width = 2820, height = 1580)

cowplot::ggdraw( grid.arrange(
  grobs = list(a, c, b, d),
  widths = c(5, 2),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)) +
  theme(plot.background = element_rect(fill="white", color = "white"))

dev.off()

## png
png("Figures/Figure_2.png",
     width = 2820, height = 1580)

cowplot::ggdraw( grid.arrange(
  grobs = list(a, c, b, d),
  widths = c(5, 2),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)) +
  theme(plot.background = element_rect(fill="white", color = "white"))

dev.off()
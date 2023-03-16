library(brms)
library(tidyverse)
library(ape)
library(geiger)
library(phytools)
library(future)
library(caper)
library(foreach)
library(doParallel)
library(glmmTMB)
library(raster)
library(lme4)
library(lmerTest)
library(gridExtra)
library(png)
library(RCurl)
library(grid)
library(tidyverse)
library(readxl)
library(plyr)
library(wesanderson)
library(hier.part)

#Read in the data
site_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "nHWI", "Mass (g)", "Forest_dependency", "Migration")

site_sp_list <- read_excel("NEE main files/First Review/Supplementatry Dataset 1.xlsx", sheet = 4)

##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## Edit Jetz name to same as tree
data$Jetz_name <- gsub(" ", "_", data$Species_name)

## Change disturbance score to numerica;
data$Disturbance_number <- ifelse(data$Disturbance == "High", 1, 0)

## Take averages at Species level
data_reduced <- data %>% group_by(Jetz_name) %>%
  dplyr::summarise(Migration = Migration[1],
            nHWI = nHWI[1], 
            Disturbance_mean = mean(Disturbance_number),
            Seasonality_mean = mean(Seasonality),
            Abs_latitude_mean = abs(mean(Latitude)))

## Scale
data_reduced$nHWI_z <- scale(data_reduced$nHWI)
data_reduced$Disturbance_mean <- scale(data_reduced$Disturbance_mean)
data_reduced$Seasonality_mean <- scale(data_reduced$Seasonality_mean)
data_reduced$Abs_latitude_mean <- scale(data_reduced$Abs_latitude_mean)


#### Split the model into different subsets according to migration
data_no_full_migs <- data_reduced %>% filter(Migration != 3)
data_sedentary_only <- data_reduced %>% filter(!Migration %in% c(2,3))
subset <- list(data_reduced, data_no_full_migs, data_sedentary_only)

## Take first subset
s <- subset[[1]]
rand_data <- as.data.frame(s)

## read tree and create a comparitive data for pgls
consensus_tree <- read.tree("Useful info/Consensus trees/consensus_tree_2021_current.nex")
comp_data <- comparative.data( consensus_tree, rand_data, "Jetz_name")

## Run PGLS
mod_i_full <- pgls(nHWI_z ~ Seasonality_mean  + Disturbance_mean + Abs_latitude_mean,
                   lambda = "ML",
                data = comp_data)


## Take statistics
summary_all_full <- summary(mod_i_full)
R2_all_full <- summary_all_full$r.squared
R2_all_full
AIC_all_full <- AIC(mod_i_full)

##############################################
## Repeat for next subset
##############################################
## pull non migrants
s <- subset[[2]]

## create comparitive data
rand_data_2 <- as.data.frame(s)
comp_data_2 <- comparative.data( consensus_tree, rand_data_2, "Jetz_name")

## run pgls
mod_no_mig_full <- pgls(nHWI_z ~ Seasonality_mean  + Disturbance_mean + Abs_latitude_mean,
                        lambda = "ML",
                              data = comp_data_2)

##pull summary stats
summary_no_mig_full <- summary(mod_no_mig_full)
R2_no_mig_full <- summary_no_mig_full$r.squared
AIC_no_mig_full <- AIC(mod_no_mig_full)

##############################################
## Repeat for next subset
##############################################
## pull residents only subset
s <- subset[[3]]

## create comparitive data
rand_data_3 <- as.data.frame(s) 
comp_data_3 <- comparative.data( consensus_tree, rand_data_3, "Jetz_name")

## run model
mod_resident_full <- pgls(nHWI_z ~  Seasonality_mean  + Disturbance_mean + Abs_latitude_mean ,
                          lambda = "ML",
                      data = comp_data_3)

## pull summary statistics
summary_resident_full <- summary(mod_resident_full)
R2_resident_full <- summary_resident_full$r.squared
AIC_resident_full <- AIC(mod_resident_full)


##############################################
#### Pull information from the model outpurs and create a df
mod_all_df <- as.data.frame(summary_all_full$coefficients)
mod_all_df <- mod_all_df[-1,]

## create CIs
colnames(mod_all_df)[2] <- "se"
mod_all_df$CI_interval <- as.numeric(mod_all_df$se) * 1.96
mod_all_df$Coefficient <- c("Seasonality", "Disturbance", "Latitude")
mod_all_df$Sample <- "All birds"

## pull coefs except intercept
mod_no_mig_df <- as.data.frame(summary_no_mig_full$coefficients)
mod_no_mig_df <- mod_no_mig_df[-1,]

## create CIs
colnames(mod_no_mig_df)[2] <- "se"
mod_no_mig_df$CI_interval <- as.numeric(mod_no_mig_df$se) * 1.96
mod_no_mig_df$Coefficient <- c("Seasonality", "Disturbance", "Latitude")
mod_no_mig_df$Sample <- "Residents and partial migrants"

## pull coefs (except intercept)
mod_resident_df <- as.data.frame(summary_resident_full$coefficients)
mod_resident_df <- mod_resident_df[-1,]

## create CIs
colnames(mod_resident_df)[2] <- "se"
mod_resident_df$CI_interval <- as.numeric(mod_resident_df$se) * 1.96
mod_resident_df$Coefficient <- c("Seasonality", "Disturbance", "Latitude")
mod_resident_df$Sample <- "Residents only"


## Combine results from models into a df for whole thing
df_all <- rbind(mod_all_df, mod_no_mig_df, mod_resident_df)


## Relevel a few things
df_all$Sample <- factor(df_all$Sample, levels = c("Residents only",
                                                  "Residents and partial migrants",
                                                  "All birds")) 

df_all$Coefficient <- factor(df_all$Coefficient, levels = c("Seasonality", "Latitude", "Disturbance")) 



## plot forest plot ##
pd <- position_dodge(0.75)
main_full <- ggplot(data=df_all, aes(x=Coefficient, y=Estimate))+ 
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(ymin=Estimate-CI_interval, ymax=Estimate+CI_interval, 
                     colour = Coefficient, alpha = Sample),lwd = 3, position = pd, show.legend=F) + 
  geom_point(aes(group = Sample, color = Coefficient), fill = "White",  size = 5,  shape = 21,  
             position = pd, key_glyph = draw_key_pointrange, show.legend=F) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  scale_color_manual(values = c("red", "deepskyblue", "gold"), 
                     breaks = c("Disturbance",
                                "Seasonality", "Latitude"),
                     guide = guide_legend(reverse = TRUE)) +
  ylab("Coefficient estimate") +
  xlab("") +
  ylim(-0.3, 0.2) +
  ggtitle("c") +
  geom_vline(xintercept=3.75, lty = 2, size = 1, col = "black", alpha = 0.6, show.legend=F)+
  coord_flip(clip = "off") +
  theme(axis.title = element_text(size = 19),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 17, hjust = 0.4),
        axis.title.x = element_blank(),
        legend.position = "none",# c(0.77,0.75),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold", vjust = -6, hjust = 0.05),
        plot.margin = margin(0,0,0,0, unit = "pt"),
        legend.key.size = unit(1, "cm")) +
  annotate("text",fontface = "bold", x=0.9, y=0.155, label= paste0("AIC: ", round(AIC_all_full, 0)), size = 5) +
  annotate("text", fontface = "bold",x=0.6, y=0.155, label= paste0("R²:  ", round(R2_all_full, 3)), size = 5)



#############################################
# Repeat above but model without Seasonality
#############################################
s <- subset[[1]]
rand_data <- as.data.frame(s)

comp_data <- comparative.data( consensus_tree, rand_data, "Jetz_name")


mod_i_DL <- pgls(nHWI_z ~  Disturbance_mean + Abs_latitude_mean,
                 lambda = "ML",
              data = comp_data)

summary(mod_i_DL)
summary_all_DL <- summary(mod_i_DL)
R2_all_LD <- summary_all_DL$r.squared
AIC_all_LD <- AIC(mod_i_DL)

##############################################
s <- subset[[2]]

rand_data_2 <- as.data.frame(s)
comp_data_2 <- comparative.data( consensus_tree, rand_data_2, "Jetz_name")

mod_no_mig_DL <- pgls(nHWI_z ~ Disturbance_mean + Abs_latitude_mean,
                      lambda = "ML",
                   data = comp_data_2)

summary_no_mig_DL <- summary(mod_no_mig_DL)
R2_no_mig_DL <- summary_no_mig_DL$r.squared
AIC_no_mig_DL <- AIC(mod_no_mig_DL)

###############################################
s <- subset[[3]]

rand_data_3 <- as.data.frame(s %>% distinct(Jetz_name, .keep_all = TRUE))
comp_data_3 <- comparative.data( consensus_tree, rand_data_3, "Jetz_name")

mod_resident_DL <- pgls(nHWI_z ~  Disturbance_mean + Abs_latitude_mean ,
                        lambda = "ML",
                     data = comp_data_3)

summary_resident_DL <- summary(mod_resident_DL)
R2_resident_DL <- summary_resident_DL$r.squared
AIC_resident_DL <- AIC(mod_resident_DL)


##############################################
mod_all_df <- as.data.frame(summary_all_DL$coefficients)
mod_all_df <- mod_all_df[-1,]

colnames(mod_all_df)[2] <- "se"
mod_all_df$CI_interval <- as.numeric(mod_all_df$se) * 1.96
mod_all_df$Coefficient <- c("Disturbance", "Latitude")
mod_all_df$Sample <- "All birds"


mod_no_mig_df <- as.data.frame(summary_no_mig_DL$coefficients)
mod_no_mig_df <- mod_no_mig_df[-1,]

colnames(mod_no_mig_df)[2] <- "se"
mod_no_mig_df$CI_interval <- as.numeric(mod_no_mig_df$se) * 1.96
mod_no_mig_df$Coefficient <- c("Disturbance", "Latitude")
mod_no_mig_df$Sample <- "Residents and partial migrants"

mod_resident_df <- as.data.frame(summary_resident_DL$coefficients)
mod_resident_df <- mod_resident_df[-1,]

colnames(mod_resident_df)[2] <- "se"
mod_resident_df$CI_interval <- as.numeric(mod_resident_df$se) * 1.96
mod_resident_df$Coefficient <- c("Disturbance", "Latitude")
mod_resident_df$Sample <- "Residents only"

df_all <- rbind(mod_all_df, mod_no_mig_df, mod_resident_df)



df_all$Sample <- factor(df_all$Sample, levels = c("Residents only",
                                                  "Residents and partial migrants",
                                                  "All birds")) 

df_all$Coefficient <- factor(df_all$Coefficient, levels = c("Seasonality", "Latitude", "Disturbance")) 



## ggplot
pd <- position_dodge(0.75)
main_dist_lat <- ggplot(data=df_all, aes(x=Coefficient, y=Estimate))+ 
  theme_classic() + #clean graph
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(ymin=Estimate-CI_interval, ymax=Estimate+CI_interval, 
                     colour = Coefficient, alpha = Sample),lwd = 3, position = pd, show.legend=F) + 
  geom_point(aes(group = Sample, color = Coefficient), fill = "White",  size = 5,  shape = 21,  
             position = pd, key_glyph = draw_key_pointrange, show.legend=F) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  scale_color_manual(values = c("red", "gold", "deepskyblue"), 
                     breaks = c("Disturbance", "Latitude",
                                "Seasonality"),
                     guide = guide_legend(reverse = TRUE)) +
  ylab("") +
  xlab("") +
  ylim(-0.3, 0.2) +
  ggtitle("b") +
  geom_vline(xintercept=2.75, lty = 2, size = 1, col = "black", alpha = 0.6) +
  coord_flip(clip = "off") +
  theme(axis.title = element_text(size = 22),
        #axis.text.x = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 17, hjust = 0.4),
        #axis.title.x = element_text(vjust = -1),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        legend.position = "none", #c(0.77,0.75),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"),
        plot.title = element_text(size = 20, face = "bold", vjust = -6, hjust = 0.05),
        plot.margin = margin(0,0,0,0, unit = "pt")) +
  annotate("text", fontface = "bold", x=0.9, y=0.155, label= paste0("AIC: ", round(AIC_all_LD, 0)), size = 5) +
  annotate("text", fontface = "bold", x=0.6, y=0.155, label= paste0("R²:  ", format(round(R2_all_LD, 3), nsmall = 3)), size = 5)



###########################################################################
# Repeat for final time but do for Disturbance only
########################################################################
## run a foreach loop which runs for each subset
s <- subset[[1]]
rand_data <- as.data.frame(s)
comp_data <- comparative.data( consensus_tree, rand_data, "Jetz_name")


mod_i_D <- pgls(nHWI_z ~Disturbance_mean,
                lambda = "ML",
              data = comp_data)

summary(mod_i_D)
summary_all_D <- summary(mod_i_D)
R2_all_dist <- summary_all_D$r.squared
AIC_all_dist <- AIC(mod_i_D)

##############################################
s <- subset[[2]]
rand_data_2 <- as.data.frame(s)
comp_data_2 <- comparative.data( consensus_tree, rand_data_2, "Jetz_name")

mod_no_mig_D <- pgls(nHWI_z ~ Disturbance_mean,
                     lambda = "ML",
                   data = comp_data_2)

summary_no_mig_D <- summary(mod_no_mig_D)
R2_no_mig_dist <- summary_no_mig_D$r.squared
AIC_no_mig_D <- AIC(mod_no_mig_D)


################################################
s <- subset[[3]]
rand_data_3 <- as.data.frame(s)
comp_data_3 <- comparative.data( consensus_tree, rand_data_3, "Jetz_name")

mod_resident_D <- pgls(nHWI_z ~ Disturbance_mean ,
                       lambda = "ML",
                     data = comp_data_3)

summary_resident_D <- summary(mod_resident_D)
R2_resident_dist <- summary_resident_D$r.squared
AIC_resident_D <- AIC(mod_resident_D)

##############################################
mod_all_df <- as.data.frame(summary_all_D$coefficients)
mod_all_df <- mod_all_df[-1,]

colnames(mod_all_df)[2] <- "se"
mod_all_df$CI_interval <- as.numeric(mod_all_df$se) * 1.96
mod_all_df$Coefficient <- c("Disturbance")
mod_all_df$Sample <- "All birds"


mod_no_mig_df <- as.data.frame(summary_no_mig_D$coefficients)
mod_no_mig_df <- mod_no_mig_df[-1,]

colnames(mod_no_mig_df)[2] <- "se"
mod_no_mig_df$CI_interval <- as.numeric(mod_no_mig_df$se) * 1.96
mod_no_mig_df$Coefficient <- c("Disturbance")
mod_no_mig_df$Sample <- "Residents and partial migrants"

mod_resident_df <- as.data.frame(summary_resident_D$coefficients)
mod_resident_df <- mod_resident_df[-1,]

colnames(mod_resident_df)[2] <- "se"
mod_resident_df$CI_interval <- as.numeric(mod_resident_df$se) * 1.96
mod_resident_df$Coefficient <- c("Disturbance")
mod_resident_df$Sample <- "Residents only"

df_all <- rbind(mod_all_df, mod_no_mig_df, mod_resident_df)



df_all$Sample <- factor(df_all$Sample, levels = c("Residents only",
                                                  "Residents and partial migrants",
                                                  "All birds")) 

df_all$Coefficient <- factor(df_all$Coefficient, levels = c("Seasonality", "Latitude", "Disturbance")) 



############ ggplot #################
pd <- position_dodge(0.75)
main_dist <- ggplot(data=df_all, aes(x=Coefficient, y=Estimate))+ 
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(ymin=Estimate-CI_interval, ymax=Estimate+CI_interval, 
                     colour = Coefficient, alpha = Sample),lwd = 3, position = pd, show.legend=F) + 
  geom_point(aes(group = Sample, color = Coefficient), fill = "White",  size = 5,  shape = 21,  
             position = pd, key_glyph = draw_key_pointrange, show.legend=F) +
  scale_color_manual(values = c("red", "gold", "deepskyblue"), 
                     breaks = c("Disturbance", "Latitude",
                                "Seasonality"),
                     guide = guide_legend(reverse = TRUE)) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  ylab("") +
  xlab("") +
  ylim(-0.3, 0.2) +
  ggtitle("a") +
  coord_flip(clip = "off") +
  theme(axis.title = element_text(size = 22),
        #axis.text.x = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 17, hjust = 0.4),
        #axis.title.x = element_text(vjust = -1),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        legend.position = "none", #c(0.77,0.5),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_rect(color = NA, fill = NA),
        plot.title = element_text(size = 20, face = "bold", vjust = -5, hjust = 0.05),
        plot.margin = margin(0,0,0,0, unit = "pt")) +
  annotate("text", fontface = "bold", x=0.9, y=0.155, label= paste0("AIC: ", round(AIC_all_dist, 0)), size = 5) +
  annotate("text", fontface = "bold", x=0.6, y=0.155, label= ifelse(R2_all_dist < 0.001,
                                                   paste0("R² < 0.001"),
                                                   paste0("R²:  ", round(R2_all_dist, 3))), size = 5)



## add an arrow specifying direction of change to dispersal
arrow <- ggplot(data=df_all, aes(x=Coefficient, y=Estimate))+ 
  theme_void() + #clean graph
  geom_linerange(aes(ymin=Estimate-CI_interval, ymax=Estimate+CI_interval),lwd = 1.5, color = "WHITE") + 
  ylim(-0.3, 0.15) +
  coord_flip(clip = "off") +
  geom_segment(
    x = 1.04, y = -0.25,
    xend = 1.04, yend = 0.16,
    lineend = "round", 
    linejoin = "round",
    size = 3.5, 
    arrow = arrow(length = unit(0.4, "inches")),
    colour = "grey80" ) +
  annotate(geom="text", x=1.625, y=-0.04, label="Dispersal limitation (nHWI)",
           color="black", size = 7) +
  theme(axis.title = element_blank(),
        axis.text = element_blank())



################################################################################
# Add hierarchical partitioning
################################################################################

### Check all possible regressions fits to check correct order
a <- all.regs(rand_data$nHWI_z, rand_data %>% 
                dplyr::select(Disturbance_mean, Abs_latitude_mean), print.vars = TRUE)

## Run the different possible PGLS combinations
mod_dist_only <- pgls(nHWI_z ~ Disturbance_mean,
                      lambda = "ML",
                      data = comp_data)

mod_lat_only <- pgls(nHWI_z ~ Abs_latitude_mean,
                     lambda = "ML",
                     data = comp_data)

mod_seas_only <- pgls(nHWI_z ~ Seasonality_mean,
                      lambda = "ML",
                      bounds=list(lambda=c(0.01,0.99)),
                      data = comp_data)

mod_seas_dist <- pgls(nHWI_z ~ Disturbance_mean + Seasonality_mean,
                      lambda = "ML",
                      data = comp_data)

mod_seas_lat <- pgls(nHWI_z ~  Seasonality_mean + Abs_latitude_mean,
                     lambda = "ML",
                     data = comp_data)

mod_lat_dist <- pgls(nHWI_z ~ Disturbance_mean + Abs_latitude_mean,
                     lambda = "ML",
                     data = comp_data)

mod_seas_dist_lat <- pgls(nHWI_z ~  Disturbance_mean + Seasonality_mean + Abs_latitude_mean,
                          lambda = "ML",
                          data = comp_data)


## Combine r.squareds as vectors put in same order as all regs vector
gfcs <- c(0, summary(mod_dist_only)$r.squared, summary(mod_seas_only)$r.squared, summary(mod_lat_only)$r.squared,
          summary(mod_seas_dist)$r.squared, summary(mod_lat_dist)$r.squared, summary(mod_seas_lat)$r.squared,
          summary(mod_seas_dist_lat)$r.squared)

## Pull appropriate data
vals_full_all <- as.data.frame(partition(gfcs, pcan = 3, var.names = c("Disturbance_mean", "Seasonality_mean", "Abs_latitude_mean"))$I.perc)
vals_full_all$Sample <- "Full sample"
vals_full_all$Model <- "Full model"
vals_full_all$variable <- row.names(vals_full_all)

## repeat with approprtiate models 
gfcs_2 <- c(0, summary(mod_dist_only)$r.squared, summary(mod_lat_only)$r.squared,
            summary(mod_lat_dist)$r.squared)

vals_full_DL <- as.data.frame(partition(gfcs_2, pcan = 2, var.names = c("Disturbance_mean",  "Abs_latitude_mean"))$I.perc)
vals_full_DL$Sample <- "Full sample"
vals_full_DL$Model <- "Disturbance & Latitude"
vals_full_DL$variable <- row.names(vals_full_DL)

hp_full <- rbind(vals_full_all, vals_full_DL)
hp_full[6,] <- c(100, "Full sample", "Disturbance only", "Disturbance_mean")



############################################################################################
### REPEAT FOR THE DIFFERENT SUBSETS
################################################################################
mod_dist_only_2 <- pgls(nHWI_z ~ Disturbance_mean,
                        lambda = "ML",
                        data = comp_data_2)

mod_lat_only_2 <- pgls(nHWI_z ~ Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_2)

mod_seas_only_2 <- pgls(nHWI_z ~ Seasonality_mean,
                        lambda = "ML",
                        data = comp_data_2)

mod_seas_dist_2 <- pgls(nHWI_z ~ Disturbance_mean + Seasonality_mean,
                        lambda = "ML",
                        data = comp_data_2)

mod_seas_lat_2 <- pgls(nHWI_z ~  Seasonality_mean + Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_2)

mod_lat_dist_2 <- pgls(nHWI_z ~ Disturbance_mean + Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_2)

mod_seas_dist_lat_2 <- pgls(nHWI_z ~  Disturbance_mean + Seasonality_mean + Abs_latitude_mean,
                            lambda = "ML",
                            data = comp_data_2)

gfcs_nm <- c(0, summary(mod_dist_only_2)$r.squared, summary(mod_seas_only_2)$r.squared, summary(mod_lat_only_2)$r.squared,
             summary(mod_seas_dist_2)$r.squared, summary(mod_lat_dist_2)$r.squared, summary(mod_seas_lat_2)$r.squared,
             summary(mod_seas_dist_lat_2)$r.squared)

vals_nomig_all <- as.data.frame(partition(gfcs_nm, pcan = 3, var.names = c("Disturbance_mean", "Seasonality_mean", "Abs_latitude_mean"))$I.perc)
vals_nomig_all$Sample <- "No full migrants"
vals_nomig_all$Model <- "Full model"
vals_nomig_all$variable <- row.names(vals_nomig_all)

gfcs_nm_2 <- c(0, summary(mod_dist_only_2)$r.squared, summary(mod_lat_only_2)$r.squared,
               summary(mod_lat_dist_2)$r.squared)

vals_nomig_DL <- as.data.frame(partition(gfcs_nm_2, pcan = 2, var.names = c("Disturbance_mean",  "Abs_latitude_mean"))$I.perc)
vals_nomig_DL$Sample <- "No full migrants"
vals_nomig_DL$Model <- "Disturbance & Latitude"
vals_nomig_DL$variable <- row.names(vals_nomig_DL)

hp_nomig <- rbind(vals_nomig_all, vals_nomig_DL)
hp_nomig[6,] <- c(100, "No full migrants", "Disturbance only", "Disturbance_mean")



############################################################################################
# REPEAT FOR THE DIFFERENT SUBSETS
################################################################################

mod_dist_only_3 <- pgls(nHWI_z ~ Disturbance_mean,
                        lambda = "ML",
                        data = comp_data_3)

mod_lat_only_3 <- pgls(nHWI_z ~ Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_3)

mod_seas_only_3 <- pgls(nHWI_z ~ Seasonality_mean,
                        lambda = "ML",
                        data = comp_data_3)

mod_seas_dist_3 <- pgls(nHWI_z ~ Disturbance_mean + Seasonality_mean,
                        lambda = "ML",
                        data = comp_data_3)

mod_seas_lat_3 <- pgls(nHWI_z ~  Seasonality_mean + Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_3)

mod_lat_dist_3 <- pgls(nHWI_z ~ Disturbance_mean + Abs_latitude_mean,
                       lambda = "ML",
                       data = comp_data_3)

mod_seas_dist_lat_3 <- pgls(nHWI_z ~  Disturbance_mean + Seasonality_mean + Abs_latitude_mean,
                            lambda = "ML",
                            data = comp_data_3)

gfcs_rs <- c(0, summary(mod_dist_only_3)$r.squared, summary(mod_seas_only_3)$r.squared, summary(mod_lat_only_3)$r.squared,
             summary(mod_seas_dist_3)$r.squared, summary(mod_lat_dist_3)$r.squared, summary(mod_seas_lat_3)$r.squared,
             summary(mod_seas_dist_lat_3)$r.squared)

vals_resident_all <- as.data.frame(partition(gfcs_rs, pcan = 3, var.names = c("Disturbance_mean", "Seasonality_mean", "Abs_latitude_mean"))$I.perc)
vals_resident_all$Sample <- "Residents only"
vals_resident_all$Model <- "Full model"
vals_resident_all$variable <- row.names(vals_resident_all)


gfcs_rs_3 <- c(0, summary(mod_dist_only_3)$r.squared, summary(mod_lat_only_3)$r.squared,
               summary(mod_lat_dist_3)$r.squared)

vals_resident_DL <- as.data.frame(partition(gfcs_rs_3, pcan = 2, var.names = c("Disturbance_mean",  "Abs_latitude_mean"))$I.perc)
vals_resident_DL$Sample <- "Residents only"
vals_resident_DL$Model <- "Disturbance & Latitude"
vals_resident_DL$variable <- row.names(vals_resident_DL)

hp_resident <- rbind(vals_resident_all, vals_resident_DL)
hp_resident[6,] <- c(100, "Residents only", "Disturbance only", "Disturbance_mean")

##############################
### Combine datasets and plot
#############################

hp <- rbind.fill(hp_full, hp_nomig, hp_resident)
colnames(hp)[1] <- "Percent_independent_effect"
hp$Percent_independent_effect <- as.numeric(hp$Percent_independent_effect)


new_data <- data.frame(Percent_independent_effect = c(0,0,0),
                       Sample = rep("Full sample", 3), 
                       Model = c("Disturbance only", "Disturbance only", "Disturbance & Latitude"),
                       variable = c("Abs_latitude_mean", "Seasonality_mean", "Seasonality_mean"))

hp2 <- rbind(hp, new_data)

hp2$Model <- factor(hp2$Model, levels = c("Disturbance only", 
                                          "Disturbance & Latitude",
                                          "Full model"))

hp2$variable[which(hp2$variable == "Abs_latitude_mean")] <- "Latitude"
hp2$variable[which(hp2$variable == "Disturbance_mean")] <- "Disturbance"
hp2$variable[which(hp2$variable == "Seasonality_mean")] <- "Seasonality"

hp2$variable <- factor(hp2$variable, levels = c("Disturbance",
                                                "Latitude",
                                                "Seasonality"))

hp2$Sample <- factor(hp2$Sample, levels = c("Full sample", 
                                            "No full migrants", 
                                            "Residents only"))


hp_plot_main <- ggplot(data = hp2 %>% dplyr::filter(
  Model == "Full model"), aes(x = variable, y = Percent_independent_effect,
                              fill = variable, alpha = Sample)) +
  ggtitle("f") +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge()) +
  scale_x_discrete(labels = c("Disturbance", "Latitude", "Seasonality")) +
  
  scale_fill_manual(values = c("red", "gold", "deepskyblue"), 
                    breaks = c("Disturbance", "Latitude",
                               "Seasonality"),
                    guide = guide_legend(reverse = TRUE)) +
  
  scale_alpha_discrete(range = c(1, 0.2)) +
  
  ylim(0,100) +
  ylab("") +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", vjust = -1.5, hjust = 0.1),
        axis.text.x = element_text(size = 17, angle = 60, hjust = 1.1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size  = 17),
        axis.title.x = element_blank(),
        plot.margin = margin(0.0,0,0.0,0, "npc"))


hp_plot_DL<- ggplot(data = hp2 %>% dplyr::filter(
  Model == "Disturbance & Latitude"), aes(x = variable, y = Percent_independent_effect, fill = variable, alpha = Sample)) +
  ggtitle("e") +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge()) +
  
  scale_fill_manual(values = c("red", "gold", "deepskyblue"), 
                     breaks = c("Disturbance", "Latitude",
                                "Seasonality"),
                     guide = guide_legend(reverse = TRUE)) +
  
  scale_alpha_discrete(range = c(1, 0.2)) +
  
  ylim(0,100) +
  ylab("Independent effect (%)") +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", vjust = -1.5, hjust = 0.1),
        #axis.title.y = element_text(size  = 16, vjust = 3),
        axis.title = element_blank(),
        plot.margin = margin(0.0,0,0.0,0, "npc"))

hp_plot_dist <- ggplot(data = hp2 %>% dplyr::filter(
  Model == "Disturbance only"), aes(x = variable, y = Percent_independent_effect, fill = variable, alpha = Sample)) +
  ggtitle("d") +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge()) +
  scale_fill_manual(values = c("red", "gold", "deepskyblue"), 
                    breaks = c("Disturbance", "Latitude",
                               "Seasonality"),
                    guide = guide_legend(reverse = TRUE)) +
  
  scale_alpha_discrete(range = c(1, 0.2)) +
  
  ylim(0,100) +
  ylab("") +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", vjust = -1.5, hjust = 0.1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size  = 17),
        axis.title.x = element_blank(),
        plot.margin = margin(0.0,0,0.0,0, "npc"))

################ create pngs and tiffs files using patchwork library ###########
### Figures are adapted outside of R - Figure legends ###
png("nHWI drivers (local).png",
    width = 900, height = 650)

((arrow / main_dist / main_dist_lat / main_full) + plot_layout(heights = c(0.5, 1,2,3)) | plot_spacer() |
    (plot_spacer() /  hp_plot_dist / hp_plot_DL/ hp_plot_main) + plot_layout(heights = c(0.25,1,1,1))) + plot_layout(widths = c(4,0.5, 1))



grid::grid.draw(grid::textGrob(hp_plot_DL$labels$y, gp = gpar(col = "black", fontsize = 20), y = 0.55, x = 0.8, rot = 90))
grid::grid.draw(grid::textGrob(main_full$labels$y, gp = gpar(col = "black", fontsize = 20), y = 0.095, x = 0.45))

dev.off()



tiff("ED - 4 nHWI drivers (local).tif",
    width = 900, height = 650)

((arrow / main_dist / main_dist_lat / main_full) + plot_layout(heights = c(0.5, 1,2,3)) | plot_spacer() |
    (plot_spacer() /  hp_plot_dist / hp_plot_DL/ hp_plot_main) + plot_layout(heights = c(0.25,1,1,1))) + plot_layout(widths = c(4,0.5, 1))



grid::grid.draw(grid::textGrob(hp_plot_DL$labels$y, gp = gpar(col = "black", fontsize = 20), y = 0.55, x = 0.8, rot = 90))
grid::grid.draw(grid::textGrob(main_full$labels$y, gp = gpar(col = "black", fontsize = 20), y = 0.095, x = 0.45))

dev.off()



library(ggExtra)
library(gridExtra)
library(grid)
library(ggplot2)
library(phytools)
library(geiger)
library(tidytree)
library(phytools)
library(phylobase)
library(ggtree)
library(ape)
library(ggnewscale)
library(boot)
library(dplyr)
library(cowplot)
library(broom)
library(readxl)
library(caper)

#phylo_data <- read.csv("Analysis/HPC folder/phylo_data2009_medium.csv")
#Read in the data
site_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "Hand.Wing.Index", "nHWI", "Mass (g)", "Forest_dependency",
                "Disturbance_range_total", "Disturbance_range_breeding",
                "Centroid_latitude_range_total (Behrmann's projection)",
                "Centroid_latitude_range_breeding (Behrmann's projection)",
                "Seasonality_range_total", "Seasonality_range_breeding")

site_sp_list <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)

##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## edit species name 
data$Jetz_name <- gsub(" ", "_", data$Species_name)

## create response variables
## restricted
data$dep_sense <- 0
data$dep_sense[data$Forest_dependency %in% c("High") & data$Classification == "Forest Core"] <- 1

## expanded
data$dep_sense_medium <- 0
data$dep_sense_medium[data$Forest_dependency %in% c("High", "Medium") & data$Classification == "Forest Core"] <- 1

## transform data to take -logs
data$nHWI_transformed <- -log(-data$nHWI)

## add suffix
data$Disturbance <- paste(data$Disturbance, "disturbance")

## summarise the tings
data1 <- data %>%
  group_by(Dataset_ID) %>%
  summarise(mean_HWI = mean(nHWI),
            mean_HWI_logged = mean(nHWI_transformed),
            abs_lat = abs(Latitude),
            Disturbed_factor = Disturbance[1],
            n = n(),
            cores = length(which(Classification == "Forest Core")),
            sensitive = length(which(dep_sense_medium == 1))) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## create numerical disturbance score
data1$Disturbed <- ifelse(data1$Disturbed_factor == "High disturbance", 1, 0)

## calculate proportions of sensitivity
data1$prop_core <- data1$cores/data1$n
data1$prop_sensitive <- data1$sensitive/data1$n


### model ###
form = as.formula("prop_sensitive ~ mean_HWI")
mod <- glm(prop_sensitive ~ mean_HWI_logged, data=data1, family = "quasibinomial")

##preditc
pred = NULL
n_pred = 100
newdata = data.frame(mean_HWI_logged = seq(min(data1$mean_HWI_logged), max(data1$mean_HWI_logged), length=n_pred))
pred_info = predict(mod, newdata=newdata, allow.new.levels=TRUE, se.fit=TRUE)

## convert output
newdata$y = 0
X = model.matrix(y ~ mean_HWI_logged, data=newdata) # model matrix (no random effects)
se_logit = sqrt(diag(X %*% vcov(mod) %*% t(X)))
pred_info[[2]] = se_logit

## df
pred = data.frame(newdata,
                  EI.sensitivity = inv.logit(pred_info[[1]]),
                  lower = inv.logit(pred_info[[1]] - qnorm(0.975)*pred_info[[2]]),
                  upper = inv.logit(pred_info[[1]] + qnorm(0.975)*pred_info[[2]]),
                  p.value = tidy(mod)$p.value[tidy(mod)$term == "mean_HWI_logged"]
)         
out = list(mod_summary=data.frame(n=nrow(data1), tidy(mod, effect="fixed"), AIC=AIC(mod)),
           pred=pred)

pred = out$pred

### plot
pmain <- ggplot() +
  geom_ribbon(data=pred, aes(x=mean_HWI_logged, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE,alpha=0.2) +
  geom_point(data = data1,  aes(y = prop_sensitive, x = mean_HWI_logged,  color = Disturbed_factor), size = 7) +
  scale_color_manual(values = c("red", "deepskyblue")) +
  #ggtitle("a") +
  geom_line(data=pred, aes(x=mean_HWI_logged, y=EI.sensitivity), color="purple4", size=1.5) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 10) ) ) +
  ylab("Fragmentation sensitivity") +
  xlab(paste0("Mean assemblage dispersal limitation (nHWI)")) +
  #ylim(c(-0.01, 0.4)) +
  #xlim(c(0.29,0.37)) +
  theme_classic() +
  #scale_y_discrete(position = "right") +
  theme(legend.position=  c(0.27, 0.875),
        #legend.direction='vertical',
        #legend.key.size = unit(0.1, "cm"),
        legend.title = element_blank(),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.text=element_text(size= 30, colour = "black",),
        axis.title.x=element_text(size= 35, vjust = 0.5),
        axis.title.y=element_text(size= 35, vjust = 2),
        axis.line = element_line(size = 3),
        plot.margin = margin(10,10,10,10, unit = "pt"),
        plot.title = element_blank(),#(size = 60, face = "bold"),
        legend.text = element_text(size = 40))

## add summary stats
pmain_annotated <-  pmain + 
  annotate("text", x=-2.8865, y=0.035, label= paste0("Coefficient: ", format(round(as.numeric(mod$coefficients[2]), 3),  nsmall = 3)), size = 12) +
  annotate("text", x=-2.87, y=-0.006, label= paste0("p-value: ", round(summary(mod)$coefficients[2,4], 3)), size = 12)

## add marginal boxplots
ybox <- axis_canvas(pmain_annotated, axis = "y") + 
  geom_boxplot(data = data1, size = 1.5, width = 2, aes(y = prop_sensitive, x = factor(Disturbed_factor), color = factor(Disturbed_factor))) + 
  scale_x_discrete() +
  scale_color_manual(values = c("red", "deepskyblue"))
xbox <- axis_canvas(pmain_annotated, axis = "x") + 
  geom_boxplot(data = data1, size = 1.5, width = 2, aes(y = mean_HWI, x = factor(Disturbed_factor), color = factor(Disturbed_factor))) +
  scale_x_discrete() +
  scale_color_manual(values = c("red", "deepskyblue")) +
  coord_flip()

p1 <- insert_xaxis_grob(pmain_annotated, xbox, grid::unit(1, "in"), position = "top")
p2 <- insert_yaxis_grob(p1, ybox, grid::unit(1, "in"), position = "right") 

# save image
jpeg("NEE main files/Figures/Extended Data/ED_figure_4.jpg",
    width = 1073, height = 894)

plot(p2)

dev.off()

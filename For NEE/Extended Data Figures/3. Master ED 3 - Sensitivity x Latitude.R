library(tidyverse)
library(ggExtra)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggplotify)
library(cowplot)
library(phytools)
library(geiger)
library(tidytree)
library(phytools)
library(tidyverse)
library(ape)
library(geiger)
library(phytools)
library(ggtree)
library(ggnewscale)
library(ggplot2)
library(gridExtra)
library(phylobase)
library(boot)
library(broom)
library(ggpubr)

## read in data
site_data <- read_excel("Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Seasonality")

species_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Order", "nHWI", "Hand.Wing.Index", "Mass (g)", "Forest_dependency")

site_sp_list <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 4) %>%
  filter(Species_name %in% species_data$Species_name)

## combine sheets
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## add suffix
data$Disturbance <- paste(data$Disturbance, "disturbance")
data$abs_lat <- abs(data$Latitude)

## create a sensitiviry score
data$Sensitive <- 0
data$Sensitive[data$Classification == "Forest Core"] <- 1

## create response variables
## restricted
data$dep_sense <- 0
data$dep_sense[data$Forest_dependency %in% c("High") & data$Sensitive == 1] <- 1

## Expanded
data$dep_sense_medium <- 0
data$dep_sense_medium[data$Forest_dependency %in% c("High", "Medium") & data$Sensitive == 1] <- 1

## Take the negative logs
## Same as log of 1/HWI
data$HWI <- log(1/(data$Hand.Wing.Index))

## Summarise the tings
data1 <- data %>%
  group_by(Dataset_ID) %>%
  summarise(mean_HWI = mean(HWI),
            abs_lat = abs_lat,
            n = n(),
            Disturbance = Disturbance, 
            cores = length(which(Sensitive == 1)),
            sensitive = length(which(dep_sense_medium == 1))) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## Calculate sensitivity proporitions
data1$prop_core <- data1$cores/data1$n
data1$prop_sensitive <- data1$sensitive/data1$n

## rename a column
data1$Disturbed_factor <- paste(data1$Disturbance)


### run model
form = as.formula("prop_sensitive ~ abs_lat")
mod <- glm(prop_sensitive ~ abs_lat, data=data1, family = "quasibinomial")

## predit model
pred = NULL
n_pred = 100
newdata = data.frame(abs_lat = seq(min(data1$abs_lat), max(data1$abs_lat), length=n_pred))
pred_info = predict(mod, newdata=newdata, allow.new.levels=TRUE, se.fit=TRUE)

## convert back to real data
newdata$y = 0
X = model.matrix(y ~ abs_lat, data=newdata) # model matrix (no random effects)
se_logit = sqrt(diag(X %*% vcov(mod) %*% t(X)))
pred_info[[2]] = se_logit

## create df
pred = data.frame(newdata,
                  EI.sensitivity = inv.logit(pred_info[[1]]),
                  lower = inv.logit(pred_info[[1]] - qnorm(0.975)*pred_info[[2]]),
                  upper = inv.logit(pred_info[[1]] + qnorm(0.975)*pred_info[[2]]),
                  p.value = tidy(mod)$p.value[tidy(mod)$term == "abs_lat"]
)         

## plot
out = list(mod_summary=data.frame(n=nrow(data1), tidy(mod, effect="fixed"), AIC=AIC(mod)),
           pred=pred)

pred = out$pred


pmain_expanded = ggplot() +
  geom_ribbon(data=pred, aes(x=abs_lat, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE,alpha=0.2) +
  geom_line(data=pred, aes(x=abs_lat, y=EI.sensitivity), color="purple4", size=1.5) +
  ylab("") +
  xlab(paste0("Absolute latitude")) +
  #ylim(c(0, 0.90)) +
  #xlim(c(-35,-15)) +
  ggtitle("b") +
  theme_classic() +
  #scale_y_discrete(position = "right") +
  theme(legend.position= "bottom",
        #legend.direction='vertical',
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank(),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.text=element_text(size= 40, colour = "black",),
        axis.title.x=element_text(size= 40, vjust = 0.5),
        axis.title.y=element_text(size= 40, vjust = 2),
        axis.line = element_line(size = 3),
        plot.margin = margin(10,10,10,10, unit = "pt"),
        plot.title = element_text(size = 50, face = "bold", vjust = -7, hjust = 0.1),
        legend.text = element_text(size = 40)) + guides(shape = guide_legend(override.aes = list(size = 10)))


pmain_expanded <-  pmain_expanded + geom_point(data = data1,  aes(y = prop_sensitive, x = abs_lat,  color = Disturbed_factor), size = 7) +
  scale_color_manual(values = c("red", "deepskyblue")) +
  annotate("text", x=50, y=0.68, label= paste0("Coefficient: ", round(as.numeric(mod$coefficients[2]), 3)), size = 15) +
  annotate("text", x=52.2, y=0.63, label= paste0("P-value: ", round(summary(mod)$coefficients[2,4], 3)), size = 15)


######################################################
## Repeat above but look for a different level of sensitivity
## Summarise our tings
data2 <- data %>%
  group_by(Dataset_ID) %>%
  summarise(mean_HWI = mean(HWI),
            abs_lat = abs_lat,
            n = n(),
            Disturbance = Disturbance, 
            cores = length(which(Sensitive == 1)),
            sensitive = length(which(dep_sense == 1))) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## clculate proportions
data2$prop_core <- data2$cores/data2$n
data2$prop_sensitive <- data2$sensitive/data2$n

## Disturbance new column name
data2$Disturbed_factor <- paste(data1$Disturbance)

## model ##
form = as.formula("prop_sensitive ~ abs_lat")
mod <- glm(prop_sensitive ~ abs_lat, data=data2, family = "quasibinomial")

## predict ##
pred = NULL
n_pred = 100
newdata = data.frame(abs_lat = seq(min(data2$abs_lat), max(data2$abs_lat), length=n_pred))
pred_info = predict(mod, newdata=newdata, allow.new.levels=TRUE, se.fit=TRUE)

## convert
newdata$y = 0
X = model.matrix(y ~ abs_lat, data=newdata) # model matrix (no random effects)
se_logit = sqrt(diag(X %*% vcov(mod) %*% t(X)))
pred_info[[2]] = se_logit

## df
pred = data.frame(newdata,
                  EI.sensitivity = inv.logit(pred_info[[1]]),
                  lower = inv.logit(pred_info[[1]] - qnorm(0.975)*pred_info[[2]]),
                  upper = inv.logit(pred_info[[1]] + qnorm(0.975)*pred_info[[2]]),
                  p.value = tidy(mod)$p.value[tidy(mod)$term == "abs_lat"]
)         

##pull output
out = list(mod_summary=data.frame(n=nrow(data2), tidy(mod, effect="fixed"), AIC=AIC(mod)),
           pred=pred)
pred = out$pred


### plot
pmain_restricted = ggplot() +
  geom_ribbon(data=pred, aes(x=abs_lat, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE,alpha=0.2) +
  geom_line(data=pred, aes(x=abs_lat, y=EI.sensitivity), color="purple4", size=1.5) +
  ylab("Fragmentation sensitivity") +
  xlab(paste0("Absolute latitude")) +
  #ylim(c(0, 0.90)) +
  #xlim(c(-35,-15)) +
  ggtitle("a") +
  theme_classic() +
  #scale_y_discrete(position = "right") +
  theme(legend.position= "bottom",
        #legend.direction='vertical',
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank(),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.text=element_text(size= 40, colour = "black",),
        axis.title.x=element_text(size= 40, vjust = 0.5),
        axis.title.y=element_text(size= 40, vjust = 2),
        axis.line = element_line(size = 3),
        plot.margin = margin(10,10,10,10, unit = "pt"),
        plot.title = element_text(size = 50, face = "bold", vjust = -7, hjust = 0.1),
        legend.text = element_text(size = 40)) 

## add points, color scheme & summary statistics
pmain_restricted <-  pmain_restricted + 
  geom_point(data = data2,  aes(y = prop_sensitive, x = abs_lat,  color = Disturbed_factor), size = 7) +
  scale_color_manual(values = c("red", "deepskyblue")) +
  annotate("text", x=50, y=0.35, label= paste0("Coefficient: ", round(as.numeric(mod$coefficients[2]), 3)), size = 15) +
  annotate("text", x=52.2, y=0.32, label= paste0("P-value: ", round(summary(mod)$coefficients[2,4], 3)), size = 15) + 
  guides(colour = guide_legend(override.aes = list(size = 15)))


## save plot
jpeg("NEE main files/Extended data figure 3.jpg", width = 2546, height = 1175)
ggarrange(pmain_restricted, pmain_expanded, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()

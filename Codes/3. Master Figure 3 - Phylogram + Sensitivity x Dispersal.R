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


#Read in the data
site_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Anthro_dist", "Natural_dist",
                "Seasonality")

species_data <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "Hand.Wing.Index", "nHWI", "Mass (g)", "Forest_dependency",
                "Disturbance_range_total", "Disturbance_range_breeding",
                "Centroid_latitude_range_total (Behrmann's projection)",
                "Centroid_latitude_range_breeding (Behrmann's projection)",
                "Seasonality_range_total", "Seasonality_range_breeding")

site_sp_list <- read_excel("Data/Supplementatry Dataset 1.xlsx", sheet = 4)

##Join them up
data <- left_join(site_sp_list, site_data)
data <- left_join(data, species_data)

## Edit species name
data$Jetz_name <- gsub(" ", "_", data$Species_name)

## Create my response variables
## Restricted
data$dep_sense <- 0
data$dep_sense[data$Forest_dependency %in% c("High") & data$Classification == "Forest Core"] <- 1

## Expanded
data$dep_sense_medium <- 0
data$dep_sense_medium[data$Forest_dependency %in% c("High", "Medium") & data$Classification == "Forest Core"] <- 1


#######################################
## Take the negative logs of HWI
data$nHWI_transformed <- -log(-data$nHWI)

## add "Disturbance" suffix
data$Disturbance <- paste(data$Disturbance, "disturbance")

## Get means per Landscape
data1 <- data %>%
  group_by(Dataset_ID) %>%
  summarise(mean_HWI = mean(nHWI),
            mean_HWI_logged = mean(nHWI_transformed),
            abs_lat = abs(Latitude),
            Disturbed_factor = Disturbance[1],
            n = n(),
            cores = length(which(Classification == "Forest Core")),
            sensitive = length(which(dep_sense == 1))) %>%
  distinct(Dataset_ID, .keep_all = TRUE) %>%
  ungroup()

## Create a numerical disturbance value
data1$Disturbed <- ifelse(data1$Disturbed_factor == "High disturbance", 1, 0)

## Calculate proportions of Core and Sensitive species
data1$prop_core <- data1$cores/data1$n
data1$prop_sensitive <- data1$sensitive/data1$n



#########
## Predict the model for plotting #
mod <- glm(prop_sensitive ~ mean_HWI_logged, data=data1, family = "quasibinomial")
pred = NULL
n_pred = 100
newdata = data.frame(mean_HWI_logged = seq(min(data1$mean_HWI_logged), max(data1$mean_HWI_logged), length=n_pred))
pred_info = predict(mod, newdata=newdata, allow.new.levels=TRUE, se.fit=TRUE)

newdata$y = 0
X = model.matrix(y ~ mean_HWI_logged, data=newdata) # model matrix (no random effects)
se_logit = sqrt(diag(X %*% vcov(mod) %*% t(X)))
pred_info[[2]] = se_logit

pred = data.frame(newdata,
                  EI.sensitivity = inv.logit(pred_info[[1]]),
                  lower = inv.logit(pred_info[[1]] - qnorm(0.975)*pred_info[[2]]),
                  upper = inv.logit(pred_info[[1]] + qnorm(0.975)*pred_info[[2]]),
                  p.value = tidy(mod)$p.value[tidy(mod)$term == "mean_HWI_logged"]
)         
out = list(mod_summary=data.frame(n=nrow(data1), tidy(mod, effect="fixed"), AIC=AIC(mod)),
           pred=pred)


pred = out$pred

## Plot model
pmain <- ggplot() +
  geom_ribbon(data=pred, aes(x=mean_HWI_logged, y=NULL, ymin=lower, ymax=upper),  fill= "purple4", show.legend=FALSE,alpha=0.2) +
  geom_point(data = data1,  aes(y = prop_sensitive, x = mean_HWI_logged,  color = Disturbed_factor), size = 7) +
  scale_color_manual(values = c("red", "deepskyblue")) +
  ggtitle("b") +
  geom_line(data=pred, aes(x=mean_HWI_logged, y=EI.sensitivity), color="purple4", size=1.5) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 10) ) ) +
  ylab("Fragmentation sensitivity") +
  xlab(paste0("Mean assemblage dispersal limitation (nHWI)")) +
  ylim(c(-0.01, 0.4)) +
  theme_classic() +
  theme(legend.position=  c(0.27, 0.875),
        legend.title = element_blank(),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length.x = unit(rel(0.5), "cm"),
        axis.ticks.length.y = unit(rel(0.5), "cm"),
        axis.text=element_text(size= 40, colour = "black",),
        axis.title.x=element_text(size= 50, vjust = 0.5),
        axis.title.y=element_text(size= 50, vjust = 2),
        axis.line = element_line(size = 3),
        plot.margin = margin(10,10,10,10, unit = "pt"),
        plot.title = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 50))

## Add the statistics
pmain_annotated <-  pmain + 
  annotate("text", x=-2.8865, y=0.014, label= paste0("Coefficient: ", round(as.numeric(mod$coefficients[2]), 3)), size = 12) +
  annotate("text", x=-2.87, y=-0.006, label= paste0("p-value: ", round(summary(mod)$coefficients[2,4], 3)), size = 12)

## Add marginal boxplots
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

plot(p2)

################################################################################
## Create new data file edit Jetz_name and create Restricted response variable
phylo_data <- data
phylo_data$Hand.Wing.Index <- -phylo_data$nHWI

phylo_data$Jetz_name <- chartr(" ", "_", phylo_data$Jetz_name)
phylo_data$dep_sense <- 0
phylo_data$dep_sense[phylo_data$Forest_dependency %in% c("High", "Medium") & phylo_data$Classification == "Forest Core"] <- 1

## Read consensus tree
tree <- read.tree("Trees/consensus_tree_2021_current.nex")

## select appropriate columns
data_df <- phylo_data %>%
  distinct(Jetz_name, .keep_all = TRUE) %>%
  dplyr::select(Jetz_name, Hand.Wing.Index, Classification, Family, Order, dep_sense)

## Create a genus column
for (i in 1:nrow(data_df)) {
  data_df$Genus[i] <- strsplit(data_df$Jetz_name, "_")[[i]][1]
}

## Summarise Genus and Family sensitivity scores and HWI
data_df2 <- data_df %>%
  group_by(Genus) %>%
  summarise(Jetz_name = Jetz_name,
            mean_HWI = mean(Hand.Wing.Index),
            prop_sensitive = mean(dep_sense),
            Family = trimws(Family),
            dep_sense = dep_sense,
            order = Order) %>%
  ungroup() %>%
  group_by(Family) %>%
  summarise(Jetz_name = Jetz_name,
            mean_HWI = mean_HWI,
            prop_sensitive_genus = prop_sensitive,
            prop_sensitive_family = mean(dep_sense),
            dep_sense = dep_sense,
            order = order,
            Genus = Genus) %>%
  ungroup()

## Order by family
data_df2 <- data_df2[with(data_df2, order(Family)),]

## Drop to Genus Tips and reorder dataframe in order of consensus tree
data_df2 <- data_df2 %>%
  distinct(Genus, .keep_all = TRUE) %>%
  dplyr::select(Jetz_name, mean_HWI, prop_sensitive_genus, prop_sensitive_family, order, Family)

Tips <- data_df2$Jetz_name
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, Tips))
data_df2<- data_df2[match(tree$tip.label, data_df2$Jetz_name),]



#--------------------------------------------------------------------------
#Create map of the HWI across a phylogeny
#--------------------------------------------------------------------------
## Convert tree format
phylo_4tree <- as(tree, "phylo4")

## add node numbers to the df for each tip
data_df2$node <- getNode(phylo_4tree, data_df2$Jetz_name)

## set row names
data_df <- as.data.frame(data_df2)
row.names(data_df) <- data_df$Jetz_name
data_df <- data_df[,-1]

## take HWI data and give species rownames
HWI_data <- as.data.frame(data_df[,1])
rownames(HWI_data) <- row.names(data_df)
colnames(HWI_data) <- "HWI" 

## Do same with orders
orders <- as.data.frame(data_df[,4])
row.names(orders) <- row.names(data_df)
colnames(orders) <- "Order"

## Do same with sensitivity score
Sensitive_data <- as.data.frame(data_df[,3])
rownames(Sensitive_data) <- row.names(data_df)
colnames(Sensitive_data) <- "Sensitive" 

## Log the HWI value for visualization
HWI_logs <- as.matrix(log(1/(data_df[,1])))

## Put HWI value next to the node id
nodes_HWI <- data.frame(node = nodeid(tree, row.names(data_df)), HWI_logs = HWI_logs)

## Fit HWI values to the tree
fit2 <- phytools::fastAnc(tree, HWI_logs, vars=TRUE, CI=TRUE)

## create new df
df <- data.frame(node = names(fit2$ace), HWI_logs = fit2$ace)

## bind nodes HWI to the df
d.1 <- rbind(nodes_HWI, df)

## change classes
d.1$node <- as.numeric(d.1$node)
d.1$HWI_logs <- as.numeric(d.1$HWI_logs)

## Fully join the tree
tree.2 <- dplyr::full_join(tree, d.1)

## Run tree with color scheme
tree_plot_colored <- ggtree(tree.2, size = 1, layout = "circular", branch.length = "none", aes(color = HWI_logs), ladderize = FALSE, continuous = FALSE)  +
  ggplot2::scale_color_gradientn(name = paste0("Dispersal limitation (nHWI)"), colors = c("dodgerblue4", "dodgerblue3", "dodgerblue2","deepskyblue", "tomato", "indianred1", "red1", "red2", "red3", "red4", "firebrick4")) 

## Create constants for plotting
k=0.01
os = 6

## Plot with tree tip images
############################################################################################################################
## Tip images will not be supplied with online code 
## Permission was granted from Birds of the World  https://birdsoftheworld.org/ - Cornell University
############################################################################################################################
tree_plot_2 <- tree_plot_colored #+
  # geom_tiplab2(aes(subset=node == jpeg_matrix$clade_node[1], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[1],".png")),  inherit.aes = FALSE, geom="image", offset=1.2*os, align=2, size= 11*k) +
  # geom_tiplab2(aes(subset=node == jpeg_matrix$clade_node[2], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[2],".png")),  inherit.aes = FALSE, geom="image", offset=1.2*os, align=2, size= 10*k) +
  # geom_tiplab2(aes(subset=node == jpeg_matrix$clade_node[3], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[3],".png")),  inherit.aes = FALSE, geom="image", offset=1.2*os, align=2, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[4], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[4],".png")), inherit.aes = FALSE, geom="image", offset=1.65*os, align=7, size= 9*k, hjust = -0.5) +
  # #geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[5], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[5],".png")), inherit.aes = FALSE, geom="image", offset= 1.6*os, align=2, size= 11*k, hjust = -1) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[6], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[6],".png")), inherit.aes = FALSE, geom="image", offset= 1.6*os, align=2, size= 9*k, hjust = -1.2) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[7], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[7],".png")), inherit.aes = FALSE, geom="image", offset= 1.5*os,align=7, size= 10*k, hjust = 0.9) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[8], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[8],".png")), inherit.aes = FALSE, geom="image", offset= 1.5*os, align=2, size= 10*k, hjust = 0.9) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[9], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[9],".png")), inherit.aes = FALSE, geom="image", offset=1.65*os, align=2, size= 12*k, hjust = 0.7) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[10], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[10],".png")), inherit.aes = FALSE, geom="image", offset=1.3*os, align=2, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[11], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[11],".png")), inherit.aes = FALSE, geom="image", offset=1.2*os, align=7, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[12], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[12],".png")), inherit.aes = FALSE, geom="image", offset=1.5*os, align=7, size= 6*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[13], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[13],".png")), inherit.aes = FALSE, geom="image", offset=1.5*os, align=7, size= 6*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[14], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[14],".png")), inherit.aes = FALSE, geom="image", offset=1.2*os, align=2, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[15], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[15],".png")), inherit.aes = FALSE, geom="image", offset=1.3*os, align=2, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[16], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[16],".png")), inherit.aes = FALSE, geom="image", offset=1.5*os, align=2, size= 10*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[17], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[17],".png")), inherit.aes = FALSE, geom="image", offset=1.9*os, align=2, size= 13*k, hjust = -0.8) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[18], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[18],".png")), inherit.aes = FALSE, geom="image", offset=2.2*os, align=7, size= 13*k, hjust = -1.1) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[19], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[19],".png")), inherit.aes = FALSE, geom="image", offset=1.75*os, align=7, size= 13*k, hjust = -0.8) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[20], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[20],".png")), inherit.aes = FALSE, geom="image", offset=1.5*os, align=7, size= 8*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[21], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[21],".png")), inherit.aes = FALSE, geom="image", offset=1.2*os, align=7, size= 8*k) +
  # geom_tiplab(aes(subset=node == jpeg_matrix$clade_node[22], image=paste0("Figures/BOW images/Phylo/",jpeg_matrix$photo_clade[22],".png")), inherit.aes = FALSE, geom="image", offset=1.2*os, align=7, size= 15*k)


## Plot heatmap
gheat_Sensitivity <- gheatmap(p = tree_plot_2, data=Sensitive_data, offset = -2, width=0.1, colnames = FALSE, color = NA) +
  new_scale(aes(color = Sensitive_data)) +
  scale_fill_gradientn(name = paste("Fragmentation sensitivity"), colors = c("grey40", "grey50","grey60", "grey70", "grey80", "yellow2","yellow1","yellow1", "yellow", "yellow"), breaks = c(0.0, 1.0)) +
  ggtitle("a") +
  theme(legend.position = c(0.5, 0.5),
        plot.title = element_text(size = 60, face = "bold"),
        legend.direction='horizontal',
        legend.box = "veritcal",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.key.height  = unit(15, "cm"),
        legend.key.width  = unit(2.5, "cm"),
        plot.margin = margin(10, 100, 80, 50, unit = "pt"),
        legend.spacing.y = unit(0.00, "cm")) +
  guides(fill = guide_colourbar(title.position="top", barheight =  1, title.hjust = 0.5, title.vjust = 1.5,  label.vjust = -0.4, order = 1),
         color = guide_colourbar(title.position="top",  barheight = 1, title.hjust = 0.5, title.vjust = 1.5, label.vjust = -0.4, order = 2))

plot(gheat_Sensitivity)


## SAVE

png("Figures/Figure 3.png",
    width = 2943/1.25, height = 1408/1.25)
grid.arrange(
  grobs = list(gheat_Sensitivity, p2),
  widths = c(2, 2),
  layout_matrix = rbind(c(1, 2))
)
dev.off()


tiff("Figures/Figure_3.tif",
     width = 2943/1.25, height = 1408/1.25)
grid.arrange(
  grobs = list(gheat_Sensitivity, p2),
  widths = c(2, 2),
  layout_matrix = rbind(c(1, 2))
)
dev.off()


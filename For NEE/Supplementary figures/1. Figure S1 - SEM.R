library(brms)
library(dplyr)
library(ape)
library(phylopath)
library(readxl)

## Read in and merge Supplementary dataset 1
site_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Disturbance", "Seasonality", "Latitude", "Longitude", "Seasonality")

species_data <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Family", "Order", "Migration", "nHWI", "Mass (g)", "Forest_dependency")

site_sp_list <- read_excel("NEE main files/Second Review/Data/Supplementatry Dataset 1_oct20.xlsx", sheet = 4)

phylo_data <- left_join(site_sp_list, site_data)
phylo_data <- left_join(phylo_data, species_data)

## Edit name so it links later down
phylo_data$Jetz_name <- chartr(" ", "_",  phylo_data$Species_name)

## Create my Restricted sensitivity score
phylo_data$dep_sense <- 0
phylo_data$dep_sense[phylo_data$Classification == "Forest Core" & phylo_data$Forest_dependency == "High"] <- 1

## Create my Expanded sensitivity score
phylo_data$dep_sense_medium <- 0
phylo_data$dep_sense_medium[phylo_data$Classification == "Forest Core" & phylo_data$Forest_dependency %in%c("Medium", "High")] <- 1

## Create a numerical Disturbance score
phylo_data$Disturbance_number <- ifelse(phylo_data$Disturbance == "High", 1, 0)

## Group by species name and get average values for Disturbance, Seasonality, Latitude, Mass and Dispersal
data_species <- phylo_data %>% dplyr::group_by(Jetz_name) %>%
  dplyr::summarise(Jetz_name = Jetz_name[1],
    prop_sensitive = mean(dep_sense),
            Disturbance = mean(Disturbance_number),
            Seasonality = mean(Seasonality),
            centroid_local = abs(mean(Latitude)),
            nHWI = mean(nHWI))

data_species <- as.data.frame(data_species)

## Create a tips vector to check against tree
Tips <- data_species$Jetz_name

## Scale all variables by sd
data_species$nHWI_z <- scale(data_species$nHWI)
data_species$Seasonality_z <- scale(data_species$Seasonality)
data_species$Disturbance_z <- scale(data_species$Disturbance)
data_species$centroid_local_z <- scale(data_species$centroid_local)

## read in phylo consensus tree
my_tree <- read.tree("Analysis/trees/consensus_tree_2021_current.nex")

## define path routes
m <- define_model_set(
  
  direct = c(prop_sensitive ∼ Disturbance_z + nHWI_z + Seasonality_z),

  indirect = c(prop_sensitive ∼ nHWI_z, 
               nHWI_z ~ Disturbance_z  + Seasonality_z),

  combo = c(nHWI_z ~ Disturbance_z +  Seasonality_z, 
             prop_sensitive ∼ Disturbance_z + nHWI_z + Seasonality_z)
  
  )

## make appropriate row names
row.names(data_species) <- data_species$Jetz_name

## run phylopath analysis
p <- phylo_path(m, data_species, my_tree)
s <- summary(p)
b <- best(p)


############### Plot results ###################
## pull necessary columns 
s_plot <- s %>% dplyr::select(model, w)

## relevel
s_plot$model <- factor(s_plot$model, levels = c("combo", "indirect", "direct"))

## plot model weightings
ggplot(s_plot, aes(x = model, y = w)) +
  geom_bar(stat = "identity", aes(fill = w), width = 0.8) +
  ylab("Model weight") + 
  scale_x_discrete(labels = c("Combination", "Indirect", "Direct")) +
  scale_fill_gradientn(name = "Model weight", colors = c("#FFF9E3", "gold")) + 
  theme_classic() +
  guides(fill = guide_colorbar(title.position = "right",
                               label.position = "left"))  +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, hjust = 1.1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.position = "none", # "right",
        legend.key.height = unit(0.148, "npc"),
        legend.key.width = unit(0.025, "npc"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, angle = 270, hjust = 0.5))

##### Model set plot is edited outside of R prior to publishing
plot_model_set(m)

##### Most parsimonious plot is edited outside of R prior to publishing
plot(b)



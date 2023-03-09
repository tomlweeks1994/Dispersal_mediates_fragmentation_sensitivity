library(brms)
library(dplyr)
library(ape)
library(geiger)
library(phytools)
library(readxl)

setwd("../data")

#setwd("C:/Users/tlw119/Dropbox/My drive/PhDizzle/MyBioFrag/Fragmentation_Project/Analysis")

#Read in the data
#phylo_data_og <- read.csv("Analysis/data/analysis_dataset_final_current_1110.csv") %>%
#  dplyr::select(PID_all, Jetz_name, Category)

site_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 2) %>%
  dplyr::select("Dataset_ID", "Latitude", "Disturbance", "Seasonality",
                "Interim_tree_cover_change", "Mean forest change", "Pixel_change_15")

species_data <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 3) %>%
  dplyr::select("Species_name", "Order", "nHWI", "Mass (g)", "Forest_dependency")

site_sp_list <- read_excel("NEE main files/First Review/Supplementatry Dataset 1_oct20.xlsx", sheet = 4) %>%
  filter(Species_name %in% species_data$Species_name)

phylo_data <- left_join(site_sp_list, site_data)
phylo_data <- left_join(phylo_data, species_data)

largest_change_dataset_IDs <- site_data$Dataset_ID[which(site_data$Interim_tree_cover_change > median(site_data$Interim_tree_cover_change))]
smallest_change_dataset_IDs <- site_data$Dataset_ID[-which(site_data$Dataset_ID %in% largest_change_dataset_IDs)]

phylo_data <- phylo_data %>% filter(Dataset_ID %in% largest_change_dataset_IDs)

#########################################
##### Standardized and center variates #########
## Log transform where appropriate (Mass and Seasonality)
## Scale predictor variables
## Scale by 2sd as comparing against binary variable ##
## Relevel the disturbance factor
#########################################
phylo_data$nHWI <- phylo_data$nHWI - mean(phylo_data$nHWI)
phylo_data$nHWI_z <- phylo_data$nHWI /(2*sd(phylo_data$nHWI))

phylo_data$Mass_log <- log(phylo_data$`Mass (g)`)
phylo_data$Mass_log <- phylo_data$Mass_log - mean(phylo_data$Mass_log)
phylo_data$Mass_log_z <- phylo_data$Mass_log / (2*sd(phylo_data$Mass_log))

phylo_data$abs_lat <- abs(phylo_data$Latitude)
phylo_data$abs_lat <- phylo_data$abs_lat - mean(phylo_data$abs_lat)
phylo_data$abs_lat_z <- phylo_data$abs_lat / (2*sd(phylo_data$abs_lat))

phylo_data$Seasonality_log <- log(phylo_data$Seasonality)
phylo_data$Seasonality_log <- phylo_data$Seasonality_log - mean(phylo_data$Seasonality_log)
phylo_data$Seasonality_z <- phylo_data$Seasonality_log / (2*sd(phylo_data$Seasonality_log))

phylo_data$Disturbance <- factor(phylo_data$Disturbance, levels = c("Low",
                                                                    "High"))

############################################
##### Create response variables ###########
############################################
## Restricted analysis=
phylo_data$dep_sense <- 0
phylo_data$dep_sense[phylo_data$Forest_dependency == "High" & phylo_data$Classification == "Forest Core"] <- 1

## Expanded analysis ##
phylo_data$dep_sense_medium <- 0
phylo_data$dep_sense_medium[phylo_data$Forest_dependency %in% c("High", "Medium") & phylo_data$Classification == "Forest Core"] <- 1

##Create a vector of expected tree tips
phylo_data$Jetz_name <-  gsub(" ", "_", phylo_data$Species_name)
Tips <- phylo_data$Jetz_name

#Read in all available Hackett back trees
bird_tree1 <- read.tree("Useful info/treesets/BirdzillaHackett1.tre")
bird_tree2 <- read.tree("Useful info/treesets/BirdzillaHackett2.tre")
bird_tree3 <- read.tree("Useful info/treesets/BirdzillaHackett3.tre")
bird_tree4 <- read.tree("Useful info/treesets/BirdzillaHackett4.tre")
bird_tree5 <- read.tree("Useful info/treesets/BirdzillaHackett5.tre")
bird_tree6 <- read.tree("Useful info/treesets/BirdzillaHackett6.tre")
bird_tree7 <- read.tree("Useful info/treesets/BirdzillaHackett7.tre")
bird_tree8 <- read.tree("Useful info/treesets/BirdzillaHackett8.tre")
bird_tree9 <- read.tree("Useful info/treesets/BirdzillaHackett9.tre")
bird_tree10 <- read.tree("Useful info/treesets/BirdzillaHackett10.tre")

#Create a list of all these treesets
alltree_list <- list(bird_tree1, bird_tree2, bird_tree3, bird_tree4, bird_tree5, bird_tree6, bird_tree7, bird_tree8, bird_tree9, bird_tree10)

#########################################CREATE THE DATA FOR FULL MODEL#############################################
#Create empty lists for brms to run through
my_phylo_list <- list()

#Populate these lists with 50 random trees and 50 of the phylo data
#set.seed(0)
for (i in 1:50) {
  print(i) #print out
  
  alltree_index <- round(runif(1, 1, length(alltree_list))) #random index
  
  treeset_i <- alltree_list[[alltree_index]] #random treeset
  
  tree_index <- round(runif(1, 1, length(treeset_i))) #random index
  
  my_tree_i <- treeset_i[[tree_index]] #random tree from random treeset
  
  my_tree_i <- ape::drop.tip(my_tree_i, setdiff(my_tree_i$tip.label, Tips)) #droptips in the tree so it is the same as my data
  
  A_i <-ape::vcv.phylo(my_tree_i) #convert to vcv matrix (brms format for phylogeny)
  
  my_phylo_list[[i]] <- list(A = A_i) #populate
}


################################################### MODELLLING ##########################################################
## These models will use the script to create 200 chains for each model, using 50 different trees and 50 identical datasets 
## in the lists created above lists this therefore requires 200 cores to run in parallel as brm. If less cores are 
## available, split the phylo_lists and data_lists into n = cores/(4*50). Run these smaller models sequentially and 
## combine using combine_models.
################################################## ALL BIRDS ############################################################

priors<- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(0, 10)", coef = "nHWI_z", class = "b"),
  set_prior("normal(0, 10)", coef = "DisturbanceHigh", class = "b"),
  set_prior("normal(0, 10)", coef = "Seasonality_z", class = "b"),
  set_prior("normal(0, 10)", coef = "abs_lat_z", class = "b")
)


priors_BM <- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(0, 10)", coef = "Mass_log_z", class = "b"),
  set_prior("normal(0, 10)", coef = "nHWI_z:Mass_log_z", class = "b"),
  set_prior("normal(0, 10)", coef = "nHWI_z", class = "b"),
  set_prior("normal(0, 10)", coef = "DisturbanceHigh", class = "b"),
  set_prior("normal(0, 10)", coef = "Seasonality_z", class = "b"),
  set_prior("normal(0, 10)", coef = "abs_lat_z", class = "b")
)



array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))

#################### BODY MASS ##############################################
full_model_medium <- brm(
  dep_sense ~ nHWI_z * Mass_log_z + Disturbance + abs_lat_z  + Seasonality_z + (1 | Dataset_ID) + (1 | Species_name) + (1 | gr(Jetz_name, cov = A)),
  data = phylo_data,
  family = bernoulli(),
  data2 = my_phylo_list[[array_number]],
  warmup = 2000,
  iter = 10000,
  backend = "cmdstan",
  chains = 4,
  threads = threading(8),
  cores = 32,
  prior = priors_BM,
  control=list(adapt_delta=0.99, max_treedepth = 12),
  file = paste0("../output/main_models_largest_change/full_model_high", array_number)
)

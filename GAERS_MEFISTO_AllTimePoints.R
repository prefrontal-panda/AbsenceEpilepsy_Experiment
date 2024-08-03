# This file is to run MEFISTO on the proteomics and metabolomics datasets across all three time points.

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#----Preparing Files----
# Loading
# proteomics
prot_16 <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Aim 1 - GAERS/Prot_SCx_16wk_noOut_timpt.csv", 
                    stringsAsFactors = F, sep = ",")
colnames(prot_16)[1] <- "sample"
prot_16 <- prot_16 %>%
  add_column(group = c(rep("GAERS",10), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(16, 21), .after = "group")

prot_7 <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Aim 1 - GAERS/Prot_SCx_7wk_noOut_timpt.csv", 
                   stringsAsFactors = F, sep = ",")
colnames(prot_7)[1] <- "sample"
prot_7 <- prot_7 %>%
  add_column(group = c(rep("GAERS",9), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(7, 20), .after = "group")

prot_3 <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Aim 1 - GAERS/Prot_SCx_3wk_noOut_timpt.csv", 
                   stringsAsFactors = F, sep = ",")
colnames(prot_3)[1] <- "sample"
prot_3 <- prot_3 %>%
  add_column(group = c(rep("GAERS",11), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(3, 22), .after = "group")

# metabolomics
met_16 <- read.csv("Aim 1 - GAERS/Met_Tha_16wk_noOut_timpt.csv",
                   stringsAsFactors = F, sep = ",")
colnames(met_16)[1] <- "sample"
met_16 <- met_16 %>%
  add_column(group = c(rep("GAERS",10), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(16, 21), .after = "group")

met_7 <- read.csv("Aim 1 - GAERS/Met_Tha_7wk_noOut_timpt.csv",
                  stringsAsFactors = F, sep = ",")
colnames(met_7)[1] <- "sample"
met_7 <- met_7 %>%
  add_column(group = c(rep("GAERS",9), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(7,20), .after = "group")

met_3 <- read.csv("Aim 1 - GAERS/Met_Tha_3wk_noOut_timpt.csv",
                  stringsAsFactors = F, sep = ",")
colnames(met_3)[1] <- "sample"
met_3 <- met_3 %>%
  add_column(group = c(rep("GAERS",10), rep("NEC",10)), .after = "sample") %>%
  add_column(week = rep(3,20), .after = "group")

# Removing missing samples from the individual datasets
# Row binding all the time points together
prot_all <- rbind(prot_3, prot_7, prot_16)
met_all <- rbind(met_3, met_7, met_16)
# Finding intersection
prot_common <- prot_all[prot_all$sample %in% met_all$sample,]
met_common <- met_all[met_all$sample %in% prot_all$sample,]

# Unit variance scaling the data (make samples rows and features columns)
prot_scale <- 
  data.frame(scale(prot_common[,-c(1:3)], center = T, scale = T)) %>%
  add_column(., prot_common[,1:3], .before = 1) %>%
  arrange(week)

met_scale <- 
  data.frame(scale(met_common[,-c(1:3)], center = T, scale = T)) %>%
  add_column(., met_common[,1:3], .before = 1) %>%
  arrange(week) 

# Creating a long data frame
# Proteomics
prot_long <- prot_scale %>% # Converting to a long dataframe
  pivot_longer(!c("sample", "group", "week"), # all columns except strain
               names_to = "feature", # proteins to go into column; 'features'
               values_to =  "value") %>% # values to go into coumn; 'value'
  arrange(feature) %>% # arranging by proteins (ascending order) 
  mutate(view = "Proteomic", .after = feature) # add 'view' column (or in this case, omic)
# Metabolomics
met_long <-met_scale %>%
  pivot_longer(!c("sample", "group", "week"),
               names_to = "feature",
               values_to = "value") %>%
  arrange(feature) %>%
  mutate(view = "Metabolomic", .after = feature)
# Concatenating both data frames
join <- rbind(prot_long, met_long)

# Saving for future use
write.csv(join, ".csv", row.names = F)

#----Preparing and Training the MOFA model----
# Loading files
df <- read.csv("MEFISTO_Input_New_Scx_20240424.csv",  stringsAsFactors = F, sep = ",")
# Creating MOFA object
# Ignoring grouping factor as we don't want to run multi-group inferences
df_nogrp <- df_filt[,-2]
# Making a new MEFISTO object
MEFobj <- create_mofa(df_nogrp)
MEFobj <- set_covariates(MEFobj,
                         covariates = "week")
plot_data_overview(MEFobj, # Viewing overview of data object
                   show_covariate = T,
                   show_dimensions = T)

# Training object
# Data options
data_opts <- get_default_data_options(MEFobj)
# If you want to change something:
data_opts$scale_views <- F
head(data_opts)

# Model options
model_opts <- get_default_model_options(MEFobj)
# Changing options
model_opts$num_factors <- 14
model_opts$spikeslab_weights <- TRUE # introduces sparsity
#model_opts$ard_weights <- TRUE
head(model_opts)

# Training options
train_opts <- get_default_training_options(MEFobj)
# Modifying training options
train_opts$maxiter <- 1500
train_opts$convergence_mode <- "medium"

# MEFISTO options
mefis_opts <- get_default_mefisto_options(MEFobj)
head(mefis_opts)
# Modifying
mefis_opts$scale_cov <- F
mefis_opts$warping <- F
head(mefis_opts)

# Now we build the final MOFA object
model_train <- prepare_mofa(object = MEFobj,
                            data_options = data_opts,
                            model_options = model_opts,
                            training_options = train_opts,
                            mefisto_options = mefis_opts)
# Training the object. We may need to tweak the above options according to the final training result
# Ideally we will select the model with the highest Evidence Lower Bound (ELBO)
outfile = file.path(getwd(),"MEFISTO_Input_New_Scx_TransFilt_20240425.hdf5")
MEFAobj.trained <- run_mofa(model_train, outfile, save_data = T, use_basilisk = T)

#----Interrogating the data----
# Loading the model
filepath <- file.path(getwd(), "MEFISTO_Input_New_Tha_TransFilt_20240425.hdf5")
trained_model <- load_model(filepath)

# Adding metadata
# Reading in file
samp_metadata <- read.csv("Aim 1 - GAERS/Sample metadata.csv", stringsAsFactors = F, sep = ",")
# Changing the s to t for thalamic data
samp_metadata$sample <- gsub("s$", "t", samp_metadata$sample)
# Excluding samples that are not in the MEFISTO model
metadata <- samp_metadata[samp_metadata$sample %in% samples_metadata(trained_model)$sample,]
# Rearranging rows according to trained_model$sample column
metadata <- metadata[order(match(metadata$sample, samples_metadata(trained_model)$sample)),]
# Rename column names
colnames(metadata) <- c("sample", "strain", "strainCode", "week", "3 week", "7 week", 
                        "16 week", "Time in Open Field Centre", "% Sucrose Preference")

# Attaching to model
samples_metadata(trained_model) <- metadata
# Checking
samples_metadata(trained_model)

# Plotting data overview
plot_data_overview(trained_model,
                   show_covariate = T,
                   show_dimensions = T)
# Checking correlation between factors
plot_factor_cor(trained_model,
                col=colorRampPalette(c("blue", "white", "red"))(200))

# Factor Exploration
# Now we can quantify the amount of variance explained by individual data modules
# We can look at either the total variance per view or per group
head(trained_model@cache$variance_explained$r2_total[[1]])
head(trained_model@cache$variance_explained$r2_per_factor[[1]])
# Plotting the amount of variance explained
# Total variance
plot_variance_explained(trained_model, x="group", y="factor", plot_total = T)[[2]]
# Per factor
plot_variance_explained(trained_model, x="view", y="factor")

# Visualising sample separation and assocciation with metadata
# Visualising factors vs covariates
plot_factors_vs_cov(trained_model,
                    factors = c(1:3),
                    covariates = "week",
                    color_by = "strain",
                    legend = T)
# Association with metadata
# Correlation plot
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode", "3 week", "7 week",
                                                 "16 week", "Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "r",
                                  col=colorRampPalette(c("blue", "white", "red"))(200),
                                  mar = c(0,0,1,0),
                                  win.asp = 0.5)
# Heatmap
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode", "3 week", "7 week",
                                                 "16 week", "Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "log_pval",
                                  alpha = 0.05,
                                  win.asp = 0.5)

# Visualisation of feature weights
# We can change the name of the features to commonly used gene name and KEGG IDs
feat_name <- features_names(trained_model) # Getting the name of features for the model

# Proteomics
p_lab <- read.csv("Prot_ED1_Genes_ProteinID_ForMEFISTO.csv", stringsAsFactors = F)
p_trained <- feat_name[["Proteomic"]] # Getting proteomics data from the trained_model
p_labMatch <- p_lab[p_lab$Protein.IDs %in% p_trained,] # Subset the main label dataframe 
p_labMatch <- p_labMatch[order(match(p_labMatch$Protein.IDs, p_trained)),] # Ordering against trained_model

# Metabolomics
m_lab <- read.csv("Met_SCx_KEGG_CompoundName.csv", stringsAsFactors = F, header = 1)
m_trained <- feat_name[["Metabolomic"]]
m_labMatch <- m_lab[m_lab$KEGG.ID %in% m_trained,]
m_labMatch <- m_labMatch[order(match(m_labMatch$KEGG.ID, m_trained)),] # Checking what is missing
m_labMatch$Compound <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", m_labMatch$Compound) # Changing some names
m_labMatch$Compound <- gsub("_",",", m_labMatch$Compound) 

# Making a list
feat_list <- list(Metabolomic = m_labMatch$Compound,
                  Proteomic = p_labMatch$Gene.Name,
                  Transcriptomic = feat_name$Transcriptomic)
# Changing feature name
features_names(trained_model) <- feat_list

# Highlighting features found in common between single-omics analysis
# Proteomics
p <- plot_top_weights(trained_model,
                      view = "Proteomic",
                      factors = 3,
                      nfeatures = 20,
                      scale = T)
# Read in the table with the common findings
#p_common <- read.csv("Prot_All_ForMEFISTO.csv", stringsAsFactors = F, row.names = 1)
p_common <- read.csv("Prot_SCx_CommonSigDE_ForMEFISTOPlot.csv", stringsAsFactors = F, row.names = 1)
p_match <- rev(unlist(p[["data"]][["feature_id"]]) %in% row.names(p_common)) # Save a logical vector of matches. Need to reverse it as MOFA starts from the bottom
p[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(p_match == T, "darkorchid2", "black") #deeppink1 # Changing colour of the y axis ticks
plot(p)

# Metabolomics
m <- plot_top_weights(trained_model,
                      view = "Metabolomic",
                      factors = 3,
                      nfeatures = 20,
                      scale = T)
m_common <- read.csv("Met_SCx_CommonSigDE_ForMEFISTOPlot.csv", stringsAsFactors = F, row.names = 1)
#read.csv("Met_All_ForMEFISTO.csv", stringsAsFactors = F, row.names = 1)
row.names(m_common) <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", row.names(m_common)) # Changing the names
row.names(m_common) <- gsub("_",",", row.names(m_common))  
m_match <- rev(unlist(m[["data"]][["feature_id"]]) %in% row.names(m_common))
m[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(m_match == T, "darkorchid2", "black")
plot(m)

# Saving graphs
save(t_scx, t_scx_commonDE, p_scx, p_scx_commonDE, p_all_scx_commonDE, m_scx, m_scx_commonDE, m_all_scx_commonDE,
     file = "GAERS_Scx_MedConverge_WeightPlots.RData")

# Plotting all graphs together
load("MEFISTO w transcripts/GAERS_Scx_MedConverge_WeightPlots.RData")
pdf("GAERS_MEFSingleComm_Plots.pdf", width = 10, height = 8)
ggarrange(MEF_Scx_Prot, MEF_Scx_met, MEF_Tha_Prot, MEF_Tha_Met,
          labels = c("A","B","C","D"),
          ncol = 2, nrow = 2, widths = c(1,1.25))
dev.off()

# Extracting the weights for the relevant factor 
weights_f2 <- get_weights(trained_model,  
                          views = "all", factors = 2, 
                          abs = F, scale = T,
                          as.data.frame = T) 
head(weights_f2, n=3L)
write.csv(weights_f2, "GAERS_Scx_MEFISTO_weights_MedConverge_20240425.csv", row.names = F)
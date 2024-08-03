# This script is for running MOFA for the GAERS 16 week data

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(MOFA2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(GGally)

## Reading in and preparing the data
# Transcriptomics
trans_16 <- read.csv("GAERS Transcriptomics/GAERS_Tha_Transcriptomics_MOFA.csv", stringsAsFactors = F, sep = ",")
colnames(trans_16)[1] <- "sample"
trans_16  <- trans_16 %>%
  add_column(group = c(rep("GAERS", 10), rep("NEC",12)), .after = "sample") %>%
  add_column(week = rep(16,22), .after = "group")
# proteomics
prot_16 <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Aim 1 - GAERS/Prot_Tha_16wk_noOut_timpt.csv", 
                    stringsAsFactors = F, sep = ",")
colnames(prot_16)[1] <- "sample"
prot_16 <- prot_16 %>%
  add_column(group = c(rep("GAERS",10), rep("NEC",12)), .after = "sample") %>%
  add_column(week = rep(16, 22), .after = "group")
# metabolomics
met_16 <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/r/Aim 1 - GAERS/Met_Tha_16wk_noOut_timpt.csv",
                   stringsAsFactors = F, sep = ",")
colnames(met_16)[1] <- "sample"
met_16 <- met_16 %>%
  add_column(group = c(rep("GAERS",10), rep("NEC",11)), .after = "sample") %>%
  add_column(week = rep(16, 21), .after = "group")
# Finding intersection
prot_common <- prot_16[prot_16$sample %in% met_16$sample,]
met_common <- met_16[met_16$sample %in% prot_16$sample,]
trans_common <- trans_16[-15,] # Removing NEC 10 from the Scx common samples, NEC 2 from Tha common samples
# Unit variance scaling the data (make samples rows and features columns). Unit variance scaling uses individual columns
prot_scale <- 
  data.frame(scale(prot_common[,-c(1:3)], center = T, scale = T)) %>%
  add_column(., prot_common[,1:3], .before = 1) %>%
  arrange(week) 
met_scale <- 
  data.frame(scale(met_common[,-c(1:3)], center = T, scale = T)) %>%
  add_column(., met_common[,1:3], .before = 1) %>%
  arrange(week)
trans_scale <- 
  data.frame(scale(trans_common[,-c(1:3)], center = T, scale = T)) %>%
  add_column(., trans_common[,1:3], .before = 1) %>%
  arrange(week)

## Creating a long data frame
# Proteomics
prot_long <- prot_scale %>% # Converting to long dataframe
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
# Transcriptomics
trans_long <- trans_scale %>%
  pivot_longer(!c("sample", "group", "week"),
               names_to = "feature",
               values_to = "value") %>%
  arrange(feature) %>%
  mutate(view = "Transcriptomic", .after = feature)
# Concatenating both data frames
join <- rbind(trans_long, prot_long, met_long)
# Saving for future use
write.csv(join, "GAERS_16wk_Tha_MOFA_20240502.csv", row.names = F)

## Creating and training the MOFA object
# Reading in data
df <- read.csv("GAERS_16wk_Tha_MOFA_20240502.csv",  stringsAsFactors = F, sep = ",")
# Creating MOFA object
# Removing some of the transcriptomic features
ordered <- read.csv("MEFISTO w transcripts/GAERS_Tha_CPMord.csv", stringsAsFactors = F, sep = ",", row.names = 1)
ordered <- as.data.frame(t(ordered))
#Removing 
remove_ids <- ordered %>% slice_tail(n = 4994)
df_filt <- df[!df$feature %in% rownames(remove_ids),]
# Ignore grouping the find factors that separates groups
df_nogrp <- df_filt[,-2]
# Making a new MEFISTO object
MOFobj <- create_mofa(df_nogrp)
plot_data_overview(MOFobj)
                   
## Training
# Data options
data_opts <- get_default_data_options(MOFobj)
data_opts$scale_views <- F
head(data_opts)
# Model options
model_opts <- get_default_model_options(MOFobj)
# Changing options
model_opts$num_factors <- 5
model_opts$spikeslab_weights <- TRUE # introduces sparsity
head(model_opts)
# Training options
train_opts <- get_default_training_options(MOFobj)

# Build the final MOFA object
model_train <- prepare_mofa(object = MOFobj,
                            data_options = data_opts,
                            model_options = model_opts,
                            training_options = train_opts)
## Saving
outfile = file.path(getwd(),"GAERS_MOFA_Input_Tha_TransFilt_20240502.hdf5")
MEFAobj.trained <- run_mofa(model_train, outfile, save_data = T, use_basilisk = T)

## Reading in the trained model
filepath <- file.path(getwd(), "GAERS 16wk MOFA/GAERS_MOFA_Input_Tha_TransFilt_20240502.hdf5")
# print(filepath)
trained_model <- load_model(filepath)

## Adding metadata to model
samp_metadata <- read.csv("Aim 1 - GAERS/Sample metadata.csv", stringsAsFactors = F, sep = ",")
samp_metadata$sample <- gsub("s$", "t", samp_metadata$sample) # Changing the s to t for thalamic data
metadata <- samp_metadata[samp_metadata$sample %in% samples_metadata(trained_model)$sample,] # Excluding samples that are not in the MEFISTO model
metadata <- metadata[order(match(metadata$sample, samples_metadata(trained_model)$sample)),] # Rearranging rows according to trained_model$sample column
metadata <- metadata[,-c(4:7)] # Removing columns that are not needed
colnames(metadata) <- c("sample", "strain", "strainCode", "Time in Open Field Centre", "% Sucrose Preference") # Rename column names
samples_metadata(trained_model) <- metadata # Attaching to model
samples_metadata(trained_model) # Checking

## Plotting data overview
plot_data_overview(trained_model)
## Checking correlation between factors
plot_factor_cor(trained_model,
                col=colorRampPalette(c("blue", "white", "red"))(200))

## Plotting the amount of variance explained
# Total variance
plot_variance_explained(trained_model, x="group", y="factor", plot_total = T)[[2]]
head(trained_model@cache$variance_explained$r2_total[[1]])
# Per factor
plot_variance_explained(trained_model, x="view", y="factor")
head(trained_model@cache$variance_explained$r2_per_factor[[1]])

## Plotting factor values
# Visualisation of single factors
plot_factor(trained_model,
            factors = c(1:2),
            color_by = "strain",
            legend = T)

## Association analysis
# Correlation plot
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode", "Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "r",
                                  col=colorRampPalette(c("blue", "white", "red"))(200),
                                  mar = c(0,0,1,0),
                                  win.asp = 0.5)
# Heatmap
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode","Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "log_pval",
                                  alpha = 0.05,
                                  win.asp = 0.5)

## Plotting Feature weights
# Changing the labels of the molecules
feat_name <- features_names(trained_model) # Getting the name of features for the model

# Proteomics
p_lab <- read.csv("Prot_ED1_Genes_ProteinID_ForMEFISTO.csv", stringsAsFactors = F)
p_trained <- feat_name[["Proteomic"]] # Getting proteomics data from the trained_model
p_labMatch <- p_lab[p_lab$Protein.IDs %in% p_trained,] # Subset the main label dataframe 
p_labMatch <- p_labMatch[order(match(p_labMatch$Protein.IDs, p_trained)),] # Ordering against trained_model
# Metabolomics
m_lab <- read.csv("Met_Tha_KEGG_CompoundName.csv", stringsAsFactors = F, header = 1)
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

# Colouring in the features that were found to be common in the single-omics
# Proteomics
p <- plot_top_weights(trained_model,
                      view = "Proteomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# Highlight proteins of interest
p_common <- c("Nit2", "Ak1", "Tmlhe")
p_match <- rev(unlist(p[["data"]][["feature_id"]]) %in% p_common)
p[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(p_match == T, "red", "black")
plot(p)

# Metabolomics
m <- plot_top_weights(trained_model,
                      view = "Metabolomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# Highlighting
m_common <- c("2-Keto-glutaramic acid", "Isopentenyl pyrophosphate", "Cystathionine")
m_match <- rev(unlist(m[["data"]][["feature_id"]]) %in% m_common)
m[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(m_match == T, "red", "black")
plot(m)

# Transcriptomics
t <- plot_top_weights(trained_model,
                      view = "Transcriptomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# Highlighting NIT2
t_common <- c("Nit2")
t_match <- rev(unlist(t[["data"]][["feature_id"]]) %in% t_common)
t[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(t_match == T, "red", "black")
plot(t)

# Getting weights
weights <- get_weights(trained_model,  
                       views = "all", factors = 1, 
                       abs = F, scale = T,
                       as.data.frame = T) 
write.csv(weights, "GAERS_Scx_MOFA_weights_20240429.csv", row.names = F)


# Saving Plots
pdf("GAERS_16wk_Tha_MOFA.pdf")
plot_data_overview(trained_model)
plot_factor_cor(trained_model,
                col=colorRampPalette(c("blue", "white", "red"))(200))
plot_variance_explained(trained_model, x="view", y="factor")
plot_factor(trained_model,
            factors = c(1:2),
            color_by = "strain",
            legend = T)
plot_factor(trained_model,
            factors = 1,
            color_by = "strain",
            legend = T)
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode", "Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "r",
                                  col=colorRampPalette(c("blue", "white", "red"))(200),
                                  mar = c(0,0,1,0),
                                  win.asp = 0.5)
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("strainCode","Time in Open Field Centre", 
                                                 "% Sucrose Preference"),
                                  plot = "log_pval",
                                  alpha = 0.05,
                                  win.asp = 0.5)
plot_top_weights(trained_model,
                 view = "Proteomic",
                 factors = 1,
                 nfeatures = 20,
                 scale = T)
plot_top_weights(trained_model,
                 view = "Metabolomic",
                 factors = 1,
                 nfeatures = 20,
                 scale = T)
plot_top_weights(trained_model,
                 view = "Transcriptomic",
                 factors = 1,
                 nfeatures = 20,
                 scale = T)
dev.off()

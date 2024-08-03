# This file is for the pre-processing of the proteomics and metabolomics data obtained from the Monash Proteomics and Metabolomics Platform (MPMP)

# Setting up Working Directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

#----- Proteomics Processing ----
# As the data has come from the MPMP group rather than direclty from the mass-spectrometer, there is no need to filer for high confidence proteins and exclude contaminants.
# Import Libraries 
library(DEP)
library(tidyverse)
# Re-labelling file 
# Loading in dataframe for analysis
prot <- read.csv("P22_0548_Exp1_ED1_unimputed_tmt_data.csv", stringsAsFactors = F, sep = ",")
# sample labels
label <- read.csv("Proteomics_SampleLabels.csv", stringsAsFactors = F, header = F)
# Changing column names for the dataset
colnames(prot) <- label[,1]
# Removing all proteins with missing values
prot_in <- na.omit(prot)

# DEP requires data to be in a SummarisedExperiment format 
# Making dataframe. Sample information should have columns; 'label', 'condition' and 'replicate'
# The label should correspond to the column names of the original dataframe.
exp <- data.frame(label = colnames(prot_in)[3:156],
                  batch = label[3:156,3],
                  condition = c(rep("G_ref",10), 
                                rep("NEC_16_scx",12), rep("GAERS_16_scx",12),
                                rep("NEC_3_scx",12), rep("GAERS_3_scx",12),
                                rep("NEC_7_scx",12), rep("GAERS_7_scx",12),
                                rep("NEC_16_tha",12), rep("GAERS_16_tha",12),
                                rep("NEC_3_tha",12), rep("GAERS_3_tha",12),
                                rep("NEC_7_tha",12), rep("GAERS_7_tha",12)),
                  replicate = c(1:10, rep(1:12, times = 12)))
# Making unique names for proteins with no Gene names - use the protein ID
data_unique <- make_unique(prot_in, "Gene Name", "Protein IDs", delim = ";")
# Getting the column number of the samples
samp_cols <- c(grep("G_ref", colnames(data_unique)), 
               grep("scx", colnames(data_unique)),
               grep("tha", colnames(data_unique)))
# Making the final SummarisedExperiment dataframe
data_se <- make_se(data_unique, samp_cols, exp)
# Missing value check
plot_frequency(data_se) # plotting barchart of protein identification overlap between samples
plot_numbers(data_se) # plot number of proteins per samples. Compare this with the one from MPMF
# Internal Reference Standard (IRS) Normalisation 
  # Adapted from: https://pwilmart.github.io/IRS_normalization/understanding_IRS.html and Plubell et al. Molecular & Cellular Proteomics, 2017
dt_long <- get_df_long(data_se)
dt_long <- dt_long[,-c(3:4, 7:9)]
# Changing the long data frame to wide data frame and pasting the values from the first two columns together
dt_wide <- dt_long %>% 
  pivot_wider(names_from = name, values_from = intensity) %>%
  unite("sample", label:batch, sep = "_", remove = T) %>%
  column_to_rownames(., var = "sample")
# Transposing. Make proteins in rows and samples in columns
dt <- as.data.frame(t(dt_wide))
# Performing exponential function on the values as they were originally log2 transformed by DEP and will be log2 transformed again
dt <- data.frame(lapply(dt, function(x) 2^x))
# Separation out reference channels
refs <- dt[,1:10]
# Finding geometric means
refs$geomean <- apply(refs, 1, function(x) exp(mean(log(x))))
# Scaling factors
refs$GRef1 <- refs$geomean / refs$G_ref_1_A
refs$GRef2 <- refs$geomean / refs$G_ref_2_B
refs$GRef3 <- refs$geomean / refs$G_ref_3_C
refs$GRef4 <- refs$geomean / refs$G_ref_4_D
refs$GRef5 <- refs$geomean / refs$G_ref_5_E
refs$GRef6 <- refs$geomean / refs$G_ref_6_F
refs$GRef7 <- refs$geomean / refs$G_ref_7_G
refs$GRef8 <- refs$geomean / refs$G_ref_8_H
refs$GRef9 <- refs$geomean / refs$G_ref_9_I
refs$GRef10 <- refs$geomean / refs$G_ref_10_J
# Applying scaling factors
prot_irs <- select(dt, ends_with("_A")) * refs$GRef1
prot_irs <- cbind(prot_irs, select(dt, ends_with("_B")) * refs$GRef2)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_C")) * refs$GRef3)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_D")) * refs$GRef4)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_E")) * refs$GRef5)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_F")) * refs$GRef6)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_G")) * refs$GRef7)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_H")) * refs$GRef8)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_I")) * refs$GRef9)
prot_irs <- cbind(prot_irs, select(dt, ends_with("_J")) * refs$GRef10)
# Putting the columns back into order according to the data frame from DEP
match_idx <- match(colnames(dt), colnames(prot_irs))
irs_reord <- prot_irs[ , match_idx]
# Putting back into DEP required object
# We need to add back the name and/or ID columns as DEP requires this
cols <- get_df_wide(data_se) # Getting the columns required
cols <- cols[,156:157] 
irs_reord <- add_column(irs_reord, cols, .before = "G_ref_1_A") # Adding column
# Cleaning
colnames(irs_reord) <- colnames(prot)
# Making unique and putting it back into DEP format
dt_unq <- make_unique(irs_reord, "Gene Name", "Protein IDs", delim = ";") 
data_irs <- make_se(dt_unq, samp_cols, exp)
#Plotting boxplot
plot_normalization(data_irs)
plot_pca(data_irs, x = 1, y = 2, indicate =  "condition")
# Variance stablising normalisation (VSN)
  # This is useful for the purposes of differential expression.
data_norm <- normalize_vsn(data_irs)
plot_normalization(data_norm) # Visualising samples before and after normalisation
plot_pca(data_norm, x = 1, y = 2, indicate =  "condition")

# Convert back to a dataframe for future analysis.
data_prepro <- get_df_wide(data_norm)
# Save a .RData object for each of future use
save(data_se, data_irs, data_norm, data_prepro, file = "Prot_SCx_MPMF_ED1.RData")
# Loading .RData file
load("Prot_Scx_noNA_ED1_2023_05_18.RData")

# Saving by time points
# Renaming the column names
colnames(data_prepro)[12:155]  <- colnames(prot)[13:156]
# Splitting by time points
prot_3wk_scx <- data_prepro[,c(157,36:59)]
prot_7wk_scx <- data_prepro[,c(157,60:83)] 
prot_16wk_scx <- data_prepro[,c(157,12:35)]

prot_3wk_tha <- data_prepro[,c(157,108:131)]
prot_7wk_tha <- data_prepro[,c(157,132:155)]
prot_16wk_tha <- data_prepro[,c(157,84:107)]
# Cleaning up the dataframes for final saving
data_final <- prot_16wk_tha %>% 
  setNames(gsub("tha_", "", names(.))) %>%
  setNames(gsub("_16 weeks", "_16t", names(.))) %>% # Including timepoint information for future analysis
  setNames(gsub(".", " ", names(.), fixed = T)) %>%
  column_to_rownames(., var = "Protein IDs") %>%
  select(order(colnames(.)))
# Transpose the data frame as we need the samples in rows and features in columns for future analysis
data_final <- as.data.frame(t(data_final))
# saving final pre-processed file
write.csv(data_final, "Prot_tha_16wk_ED1_noNA_timpt.csv", row.names = T)

#----- Metabolomics Processing ----
# For this, we will run our own pre-processing.
  # The qmtools package will be used. Briefly, the steps are:
  # (i) Features with > 50% missing values within a group are removed.
  # (ii) Imputation using quantile regression imputation of left-censored data (QRILC) 
  # (iii) Normalisation using Probababilistic Quotient Normalization (PQN) 
  # (iv) Pareto scaling and log10 transformation
# Samples are required to be in columns and features in rows.
# KEGG ID matching is performed and columns with no matches are removed.
m <- read.csv("//ad.monash.edu/home/User028/dcho0009/Desktop/Tha_Dupes_NAPresent.csv",
              stringsAsFactors = F)
# Removing columns
m1 <- m %>% select(-contains(c("X", "NA")))
# Saving
write.csv(m1, "Met_tha_noNA.csv", row.names = F)

#Now we run the pre-processing
# Loading libraries
library(qmtools)
library(tidyverse)
library(ggplot2)
library(naniar)
library(DEP)
# Reading files
met <- read.csv("Met_tha_noNA_noDupes.csv", stringsAsFactors = F)
# Extracting the required matrices
# Data matrix. should include both sample names as rows and variable names as columns
peak <- met %>%
  column_to_rownames(., var = "sample") %>%
  select(-label) 
# Converting NA# to NA so it can be read as missing properly. We can also coerce NAs so any non-numerics become NA
peak=data.frame(lapply(peak,as.numeric),check.names = FALSE)
rownames(peak) <- met[,1]
# transposing the file so that the samples are in columns
peak_t <- as.data.frame(t(peak))
# log10 transformation as qmtools is unable to do it
# log10 transformation
pk_lg <- log10(peak_t)
# Intensity 
pk_int <- as.matrix(pk_lg)
# Now we make the SummarisedExperiment class object for qmtools
# Feature metadata
features <- data.frame(Metabolites = colnames(peak))
# Sample metadata
sm <- data.frame(label = rownames(peak),
                 strain = c(rep("GAERS", 34), rep("NEC",36), rep(" ", 11)),
                 timept = c(rep("16 weeks", 12), rep("3 weeks",10), rep("7 weeks", 12),
                            rep("16 weeks", 12), rep("3 weeks",12), rep("7 weeks", 12),
                            rep(" ", 11)),
                 id = c(rep("sample", 70), rep("QC",11)))
# Converting the columns to factors            
sm[sapply(sm, is.character)] <- lapply(sm[sapply(sm, is.character)], 
                                       as.factor)
str(sm)
# Making SummarisedExperiment class
data <- SummarizedExperiment(assays = pk_int,
                             rowData = features,
                             colData = sm)
assayNames(data) <- "raw"
data
# Visualising density distribution of missing values to decide on an imputation method
  # For these plots, change the x-axis label for the DEP density plot.
  # Then copy paste the plot code into thet console for saving.
  # Otherwise, the DEP plot will be printed twice and the x-axis label will be wrong.
pdf("Met_SCx_Missing Distribution.pdf")
# Missing value plot
visdat::vis_miss(peak_t) + 
  theme(axis.text.x = element_text(angle = 90)) # add ggplot2 aesthetics to the plot 
# Density plot for missing values
# placing it into a variable to view the plot
p <- plot_detect(data) 
# changing the x-axis label, Check the GRID.text number!
p[["grobs"]][[1]]$grobs[[12]]$children$GRID.text.95$label <- expression(log[10] ~ "Intensity") 
# original: expression(log[2] ~ "Intensity")
# do the same for the CumSum plot
p[["grobs"]][[2]]$grobs[[12]]$children$GRID.text.153$label <- expression(log[10] ~ "Intensity")
# original: expression(log[2] ~ "Intensity")
# Plotting the modified plot
plot(p)
dev.off()
# Perform imputation and normalisation
# Visualising the missing values in the dataset
# Sample group information that we want to plot
g <- factor(colData(data)$strain)
# Filtering
# We can check the dimensions of the dataset after filtering
dim(removeFeatures(data, i = "raw",
                   group = colData(data)$strain,
                   cut = 0.50))
# returns 578 76 which is the original dimension of the dataset for scx
# for tha: 576 81
# Imputation
se <- imputeIntensity(data, i = "raw",
                      method = "QRILC",
                      name = "QLC")
# PQN
se <- normalizeIntensity(se, i = "QLC",
                         method = "pqn", ref_samples = c(71:81),
                         name = "QLC_PQN")
# We can also visualise the results before and after normalisation
pdf("Met_SCx_DataDistribution.pdf")
# before normalization
plotBox(se, i = "QLC", group = g, log2 = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "Group", labels = c("QC", "GAERS", "NEC")) +
  ggtitle("Before Normalisation")
# after normalization
plotBox(se, i = "QLC_PQN", group = g, log2 = F)  +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "Group", labels = c("QC", "GAERS", "NEC")) +
  ggtitle("After Normalisation")
dev.off()
# Pareto scaling
se <- normalizeIntensity(se, i = "QLC_PQN",
                         method = "feature.scale",
                         type = "pareto",
                         name = "QLC_PQN_PS")
# Saving as an RData file
saveRDS(se, "Met_Scx_QLC.rds")
se <- readRDS("Met_Scx_QLC.rds")
# Exporting file
dt_fin <- as.data.frame(assay(se, "QLC_PQN_PS"))
# Changing sample labels
samp <- read.csv("scx_labels.csv")
colnames(dt_fin) <- samp[,2]
# Splitting by time point
met_16s <- as.data.frame(dt_fin[,c(1:12,35:46)])
met_7s <- as.data.frame(dt_fin[,c(24:34,58:69)])
met_3s <- as.data.frame(dt_fin[,c(13:23,47:57)])
# Thalamus
met_16t <- as.data.frame(dt_fin[,c(1:12,35:46)])
met_7t <- as.data.frame(dt_fin[,c(23:34,59:70)])
met_3t <- as.data.frame(dt_fin[,c(13:22,47:58)])
# Making into final dataset
data_final <- met_3t %>%
  setNames(gsub("_tha", "", names(.))) %>%
  setNames(gsub("3 weeks", "3t", names(.))) %>% # Including timepoint information for future analysis
  select(order(colnames(.)))
# Transposing
data_final <- as.data.frame(t(data_final))
# Saving
write.csv(data_final, "Met_Tha_3wk_timpt.csv")
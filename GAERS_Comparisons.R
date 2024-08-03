# This is to perform comparison of the single-omics results across all time points and groups.

# Setting wrorking directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(ComplexHeatmap)
library(ggbeeswarm)
library(viridis)

# Loading all files
# Proteomics
# Scx
SO_p16 <- read.csv("Prot_Scx_16wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_p7 <- read.csv("Prot_Scx_7wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_p3 <- read.csv("Prot_Scx_3wks_DEAll_FC1.5.csv", stringsAsFactors = F)
# Filter out non-significant proteins
SO_p16 <- SO_p16 %>% filter(SO_p16$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_p7 <- SO_p7 %>% filter(SO_p7$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_p3 <- SO_p3 %>% filter(SO_p3$Sig.DE!="Not Sig") %>% rename(feature = X)
# Thalamus
SO_p16 <- read.csv("Prot_Tha_16wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_p7 <- read.csv("Prot_Tha_7wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_p3 <- read.csv("Prot_Tha_3wks_DEAll_FC1.5.csv", stringsAsFactors = F)
# Filter out non-significant proteins
SO_p16 <- SO_p16 %>% filter(SO_p16$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_p7 <- SO_p7 %>% filter(SO_p7$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_p3 <- SO_p3 %>% filter(SO_p3$Sig.DE!="Not Sig") %>% rename(feature = X)

# Metabolomics
# Scx
SO_m16 <- read.csv("Met_SCx_16wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_m7 <- read.csv("Met_SCx_7wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_m3 <- read.csv("Met_SCx_3wks_DEAll_FC1.5.csv", stringsAsFactors = F)
# Filter out non-significant proteins
SO_m16 <- SO_m16 %>% filter(SO_m16$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_m7 <- SO_m7 %>% filter(SO_m7$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_m3 <- SO_m3 %>% filter(SO_m3$Sig.DE!="Not Sig") %>% rename(feature = X)
# Thalamus
SO_m16 <- read.csv("Met_Tha_16wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_m7 <- read.csv("Met_Tha_7wks_DEAll_FC1.5.csv", stringsAsFactors = F)
SO_m3 <- read.csv("Met_Tha_3wks_DEAll_FC1.5.csv", stringsAsFactors = F)
# Filter out non-significant proteins
SO_m16 <- SO_m16 %>% filter(SO_m16$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_m7 <- SO_m7 %>% filter(SO_m7$Sig.DE!="Not Sig") %>% rename(feature = X)
SO_m3 <- SO_m3 %>% filter(SO_m3$Sig.DE!="Not Sig") %>% rename(feature = X)

#----Finding Common Features----
# Proteomics
# Making a list of the dataframes
list_p <- list(SO_p3, SO_p7, SO_p16)
names(list_p) <- c("3 week", "7 week", "16 week")
# Finding common proteins
common_p <- list_p %>% reduce(inner_join, by = join_by(feature))
names(common_p) <- gsub("\\.x", "_3wk", names(common_p))
names(common_p) <- gsub("\\.y", "_7wk", names(common_p))
write.csv(common_p, "Comparisons/Prot_Tha_CommonSigDE.csv", row.names = F)
# Metabolomics
# List
list_m <- list(SO_m3, SO_m7, SO_m16)
names(list_m) <- c("3 week", "7 week", "16 week")
# Commonalities
common_m <- list_m %>% reduce(inner_join, by = join_by(feature))
names(common_m) <- gsub("\\.x", "_3wk", names(common_m))
names(common_m) <- gsub("\\.y", "_7wk", names(common_m))
write.csv(common_m, "Comparisons/Met_Scx_CommonSigDE.csv", row.names = F)

#----Barplot----
# Loading files
prot_bar <- read.csv("Comparisons/Prot_Scx_CommonSigDE.csv", header = 1)
feat <- read.csv("Prot_ED1_Genes_ProteinID.csv")
# Matching the Accession ID with the protein name
match <- prot_bar %>% 
  left_join(feat, by = c("feature" = "Protein.IDs")) %>%
  mutate(feature = coalesce(Gene.Name,feature), .keep = "none")
prot_bar[,1] <- match
# Making a grouped bar chart
p_3wk <- prot_bar %>%
  select(c("feature", "logFC_3wk")) %>%
  mutate(`Time Point (week)` = "3", .before = 2) %>%
  rename(`Log fold-change` = `logFC_3wk`)
p_7wk <- prot_bar %>%
  select(c("feature", "logFC_7wk")) %>%
  mutate(`Time Point (week)` = "7", .before = 2) %>%
  rename(`Log fold-change` = `logFC_7wk`)
p_16wk <- prot_bar %>%
  select(c("feature", "logFC")) %>%
  mutate(`Time Point (week)` = "16", .before = 2) %>%
  rename(`Log fold-change` = `logFC`)
p_bar_a <- rbind(p_3wk, p_7wk, p_16wk) %>%
  mutate(ID = c(rep(1,18), rep(2,18), rep(3,18)))
# Making bar plot
ggplot(p_bar_a, 
       aes(fill = `Time Point (week)`, group = -ID,
           y = `Log fold-change`, x = reorder(feature, `Log fold-change`))) +
  ylab(expression("Log fold-change")) + 
  xlab(expression("Significantly Differentially Expressed Proteins Across Time Points (SCx)")) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  scale_fill_viridis(discrete = T, option = "viridis", breaks = c("3","7","16")) +
  geom_bar(position = "dodge", stat = "identity") + # position used to make bars next to each other
  coord_flip()

# Metabolomics
met_bar <- read.csv("Met_Scx_CommonSigDE.csv", header = 1)
kegg <- read.csv("Met_Scx_KEGG_CompoundName.csv")
kegg$Compound <- gsub("_",",", kegg$Compound)
kegg$Compound <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", kegg$Compound)
# Matching the Accession ID with the protein name
matchM <- met_bar %>% 
  left_join(kegg, by = c("feature" = "KEGG.ID")) %>%
  mutate(feature = coalesce(Compound,feature), .keep = "none")
met_bar[,1] <- matchM 
# Making a grouped bar chart
m_3wk <- met_bar %>%
  select(c("feature", "logFC_3wk")) %>%
  mutate(`Time Point (week)` = "3", .before = 2) %>%
  rename(`Log fold-change` = `logFC_3wk`)
m_7wk <- met_bar %>%
  select(c("feature", "logFC_7wk")) %>%
  mutate(`Time Point (week)` = "7", .before = 2) %>%
  rename(`Log fold-change` = `logFC_7wk`)
m_16wk <- met_bar %>%
  select(c("feature", "logFC")) %>%
  mutate(`Time Point (week)` = "16", .before = 2) %>%
  rename(`Log fold-change` = `logFC`)
m_bar_a <- rbind(m_3wk, m_7wk, m_16wk) %>%
  mutate(ID = c(rep(1,10), rep(2,10), rep(3,10))) 
# Making bar plot
ggplot(m_bar_a, 
       aes(fill = `Time Point (week)`, group = -ID,
           y = `Log fold-change`, x = reorder(feature, `Log fold-change`))) +
  ylab(expression("Log fold-change")) + 
  xlab(expression("Significantly Differentially Expressed Metabolites Across Time Points (Scx)")) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  scale_fill_viridis(discrete = T, option = "viridis", breaks = c("3","7","16")) +
  geom_bar(position = "dodge", stat = "identity") + # position used to make bars next to each other
  coord_flip()

#----UpSet Plot----
# Proteomics
# First, we need to combine the data we want into a single matrix
# Making a binary matrix from the a list of interest
intList <- list(Scx_3wk = SO_p3$feature,
                Scx_7wk = SO_p7$feature,
                Scx_16wk = SO_p16$feature,
                Tha_3wk = SO_pt3$feature,
                Tha_7wk = SO_pt7$feature,
                Tha_16wk = SO_pt16$feature)
# Making the combination matrix
plot_mat <- make_comb_mat(intList, mode = "distinct") # Same as: make_comb_mat(list_to_matrix(intList))
plot_mat # Looking at the details
# Plot
pdf("Prot_All_CommonDE.pdf", width = 9, height = 5)
UpSet(plot_mat,
      set_order = c("Scx_3wk", "Tha_3wk", "Scx_7wk", "Tha_7wk", "Scx_16wk", "Tha_16wk"),
      comb_order = order(comb_degree(plot_mat), decreasing = T),
      top_annotation = upset_top_annotation(plot_mat, add_numbers = T),
      right_annotation = upset_right_annotation(plot_mat, add_numbers = T),
      column_title = "Common Proteins (Brain Region and Timepoint)")
dev.off()
# Compiling the commonalities
allList <- list(SO_p3, SO_p7, SO_p16, SO_pt3, SO_pt7, SO_pt16)
comm <- allList %>% 
  reduce(inner_join, by = join_by(feature)) %>%
  column_to_rownames("feature")
# Renaming columns
names(comm) <- gsub("\\.y\\.y\\.y", "_tha_16wk", names(comm))
names(comm) <- gsub("\\.x\\.x\\.x", "_tha_7wk", names(comm))
names(comm) <- gsub("\\.y\\.y", "_tha_3wk", names(comm))
names(comm) <- gsub("\\.x\\.x", "_scx_16wk", names(comm))
names(comm) <- gsub("\\.y", "_scx_7wk", names(comm))
names(comm) <- gsub("\\.x", "_scx_3wk", names(comm))
# Getting columns we want out
comm_all <- comm %>%
  select(grep("logFC", names(comm)), grep("AveExpr", names(comm)),
         grep("Sig.DE", names(comm)))
# Saving
write.csv(comm_all, "Prot_All_CommonDE.csv", row.names = T)

# Metabolomics
# Making the combination matrix
intList <- list(Scx_3wk = SO_m3$feature,
                Scx_7wk = SO_m7$feature,
                Scx_16wk = SO_m16$feature,
                Tha_3wk = SO_mt3$feature,
                Tha_7wk = SO_mt7$feature,
                Tha_16wk = SO_mt16$feature)
plot_mat <- make_comb_mat(intList, mode = "distinct") 
plot_mat # Looking at the details
# Plot
pdf("Met_All_CommonDE.pdf", width = 9, height = 5)
UpSet(plot_mat,
      set_order = c("Scx_3wk", "Tha_3wk", "Scx_7wk", "Tha_7wk", "Scx_16wk", "Tha_16wk"),
      comb_order = order(comb_degree(plot_mat), decreasing = T),
      top_annotation = upset_top_annotation(plot_mat, add_numbers = T),
      right_annotation = upset_right_annotation(plot_mat, add_numbers = T),
      column_title = "Common Metabolites (Brain Region and Timepoint)")
dev.off()
# Compiling the commonalities
allList <- list(SO_m3, SO_m7, SO_m16, SO_mt3, SO_mt7, SO_mt16)
comm <- allList %>% 
  reduce(inner_join, by = join_by(feature)) %>%
  column_to_rownames("feature")
# Renaming columns
names(comm) <- gsub("\\.y\\.y\\.y", "_tha_16wk", names(comm))
names(comm) <- gsub("\\.x\\.x\\.x", "_tha_7wk", names(comm))
names(comm) <- gsub("\\.y\\.y", "_tha_3wk", names(comm))
names(comm) <- gsub("\\.x\\.x", "_scx_16wk", names(comm))
names(comm) <- gsub("\\.y", "_scx_7wk", names(comm))
names(comm) <- gsub("\\.x", "_scx_3wk", names(comm))
# Getting columns we want out
comm_all <- comm %>%
  select(grep("logFC", names(comm)), grep("AveExpr", names(comm)),
         grep("Sig.DE", names(comm)))
# Saving
write.csv(comm_all, "Comparisons/Met_All_CommonDE.csv", row.names = T)
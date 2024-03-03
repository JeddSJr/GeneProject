---
title: "Updated_LUSC_CPTAC_GSEA"
author: "GSEA"
date: "2024-03-2"
output: html_document
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)


```
# Install BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of CRAN packages
cran_packages <- c("tidyverse", "ggplot2", "magrittr", "dplyr", "tidyr", "tibble", 
                   "stringr", "ggpubr", "gridExtra", "reshape2", "RColorBrewer", 
                   "circlize", "ggtext", "magick", "grid", "multipanelfigure", 
                   "corrplot", "msigdbr", "pathview", "enrichplot", "GOSemSim", 
                   "ggdendro", "ggridges", "cowplot", "ggplot2", "readr")

# Install CRAN packages that are not already installed
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# List of Bioconductor packages
bioc_packages <- c("ComplexHeatmap", "EnrichmentBrowser", "clusterProfiler", 
                   "GSEAmining", "AnnotationDbi", "BiocStyle", "MutationalPatterns", 
                   "maftools", "BSgenome.Hsapiens.UCSC.hg19", "DelayedArray", "NMF")

# Install Bioconductor packages that are not already installed
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}



```{r}
library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ggtext)
library(magick)
library(grid)
library(multipanelfigure)
library(corrplot)
library(EnrichmentBrowser)
library(msigdbr)
library(clusterProfiler)
library(pathview)
library(GSEAmining)
library(enrichplot)
library(RColorBrewer)
library(GOSemSim)
library(ggdendro)
library(ggridges)
library(cowplot)
library(AnnotationDbi)
library(ggplot2)
library(readr)
library(BiocStyle)
library(MutationalPatterns)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DelayedArray)
library(NMF)
library(dplyr)

```


```{r}
# Set working directory.
# Adjust as needed
setwd("/Users/awsalmirahmad/Downloads/TMM30092/lusc_cptac_2021")
file_loc <- "/Users/awsalmirahmad/Downloads/TMM30092/lusc_cptac_2021"



# Read in the csv file
df <- read.csv('Final_LUSC_CPTAC_data_mrna_seq_tpm.csv', stringsAsFactors=FALSE)

# Create new column X by pasting Hugo_Symbol and Entrez_Gene_Id together
df$X <- paste(df$Hugo_Symbol, "(", df$Entrez_Gene_Id, ")", sep = "")

# Make X the first column
df <- df[, c(ncol(df), 1:(ncol(df)-1))]

df <- df[,-c(2,3)]

# Write the updated dataframe back to the csv file
write.csv(df, 'CPTAC_TransformedX_data_mrna_seq_tpm123.csv', row.names = FALSE)



# Read the CSV file
#data <- read.csv("TransformedX_data_mrna_seq_tpm.csv")

# Save the column names before transposing
col_names <- colnames(df)

# Transpose the dataframe
transposed_data <- as.data.frame(t(df))

# Replace the row names (which are now "V1", "V2", etc.) with the original column names
rownames(transposed_data) <- col_names

dim(transposed_data)

# Set the column names of the transposed data to the first row
colnames(transposed_data) <- transposed_data[1, ]

# Then remove the first row
transposed_data12 <- transposed_data[-1, ]

class(transposed_data12)
dim(transposed_data12)
```


```{r}

write.csv(transposed_data12, "CPTAC_transposed123_data_tpm.csv")
matrix_original77 <- read.csv("CPTAC_transposed123_data_tpm.csv")
dim(matrix_original77)


```
```{r}
SBS_table <- read.csv("SBS_relative_cptac_contribution.csv")

SBS_table <- SBS_table %>% filter(Sample_ID %in% matrix_original77$X)

rownames(SBS_table) <- SBS_table$Sample_ID
```
```{r}
##################
### check correlation between the SBS values
res <- cor(SBS_table[,2:7], method="pearson")

# same but saving as multi panel figure
my_cor_plot_m <- capture_base_plot(corrplot(	res, 
                                             type = "upper", 
                                             order = "hclust", 
                                             tl.col = "black", 
                                             tl.srt = 45,
                                             tl.cex = 0.5))

montage <- multi_panel_figure(
  width = 8,
  height = 8,
  rows = 1,
  columns = 1,
  row_spacing = 0.0,
  column_spacing = 0.0,
  unit = "in",
  figure_name = "Correlation",
  panel_label_type = "none")

montage %<>% 
  fill_panel(my_cor_plot_m, row = 1, column = 1)

# Try adding a title
montage <- annotate_figure(montage, top = text_grob(paste0("Correlation between Cancer samples\nNormalized exposure values for SBS signatures"), color = "black", face = "bold", size = 16))

# save as pdf
ggsave(paste("Correlation_plot.Normalized_exposures_to_SBS_signatures", ".02.pdf", sep=""), plot=montage, width=9, height=9, unit="in")

```
```{r}

#######################
### Find genes most correlated with SBS1 expo values

# Cleaning up the expression matrix
# Reformat the expression data matrix
matrix_temp1 <- matrix_original77
# set row names as the cell line
rownames(matrix_temp1) <- matrix_temp1[,1]
# Remove column with cell line name
matrix_temp1$X <- NULL

library(stringr)

# Splitting by "." and taking all but the last element as the gene symbol
gene_symbol <- sapply(strsplit(colnames(matrix_temp1), "\\."), function(x) {paste(head(x, -1), collapse = ".")})

# Splitting by "." and taking the last element as the Entrez ID
entrez_id <- sapply(strsplit(colnames(matrix_temp1), "\\."), function(x) {tail(x, n=1)})

gene_info <- data.frame(ID=colnames(matrix_temp1), Symbol=gene_symbol, Entrez=entrez_id)
rownames(gene_info) <- gene_info$ID


# Retain only cancer samples (rows) present in both studies 
allGenes_allSamples <- filter(matrix_temp1, rownames(matrix_temp1) %in% SBS_table$Sample_ID)

# Ensure allGenes_allSamples is a data frame or matrix
# If it's not, convert it. This example assumes it's a matrix:
# allGenes_allSamples <- data.frame(allGenes_allSamples)

# Check if allGenes_allSamples has row names
if(is.null(rownames(allGenes_allSamples))) {
  # If not, you need to assign row names. 
  # This step depends on your specific data structure.
  # For example, if you have a vector of names:
  # rownames(allGenes_allSamples) <- vector_of_names
}

# Now, ensure that SBS_table$Sample_ID contains values that match the row names of allGenes_allSamples
# You can check for mismatches like this:
mismatches <- !SBS_table$Sample_ID %in% rownames(allGenes_allSamples)
if(any(mismatches)) {
  cat("The following Sample_IDs do not match row names in allGenes_allSamples:\n")
  print(SBS_table$Sample_ID[mismatches])
  # You may need to correct these mismatches before proceeding
}

# Assuming now that allGenes_allSamples has proper row names and they match with SBS_table$Sample_ID
allGenes_allSamples <- allGenes_allSamples[SBS_table$Sample_ID, ]


# Make in the same order as the SBS_table
#allGenes_allSamples <- allGenes_allSamples[SBS_table$Sample_ID, ]





```
```{r}
# Calculate correlation for all genes with SBS1, , "SBS4"
cor_genes_to_sbs1 <- as.data.frame(cor(allGenes_allSamples, SBS_table[, c("SBS1")], method="pearson"))
colnames(cor_genes_to_sbs1) <- c("Corr_to_SBS1")


# sort in descending order
cor_genes_to_sbs1 <- arrange(cor_genes_to_sbs1, desc(Corr_to_SBS1)) 

write.csv(cor_genes_to_sbs1, "Main_CPTAC_Cor_genes_to_SBS1.csv")
```

```{r}
# make plot and save
histo_cor <- capture_base_plot(hist(	cor_genes_to_sbs1$Corr_to_SBS1, 
                                     breaks=100, 
                                     main="Distribution of correlation coefficients\nbetween gene expression and SBS1 values,\nin 484 cancer samples for 20531 genes",
                                     xlab="Pearson correlation coefficient",
                                     ylab="Number of genes"))

montage_histogram <- multi_panel_figure(
  width = 8,
  height = 9,
  rows = 1,
  columns = 1,
  row_spacing = 0.2,
  column_spacing = 0.2,
  unit = "in",
  figure_name = "SBS1 correlation with gene expression in cancer samples",
  panel_label_type = "none")

montage_histogram %<>% 
  fill_panel(histo_cor, row = 1, column = 1)



# save as pdf
ggsave(paste("Histogram_SBS1_correlation_distribution.", "pdf", sep=""), plot=montage_histogram, width=9, height=10, unit="in")



```
```{r}
# Load necessary libraries
library(dplyr)      # For data manipulation
library(msigdbr)    # For MSigDB gene sets

# Assuming cor_genes_to_sbs1 is your dataset with correlations
geneList_unsorted <- cor_genes_to_sbs1$Corr_to_SBS1

# Assuming you have a reason to match genes to entrez_id, but the code was commented out
# If you need to extract entrez_id and it's located as described, here's a more concise way:
# Just ensure your 'cor_genes_to_sbs1' data frame has proper rownames or adjust accordingly
names(geneList_unsorted) <- sapply(strsplit(rownames(cor_genes_to_sbs1), "\\."), function(x) x[2])

# Sort geneList_unsorted in decreasing order
my_geneList <- sort(geneList_unsorted, decreasing = TRUE)

# Write the sorted gene list to a CSV file
write.csv(my_geneList, "my_CPTAC_geneList_SBS1.csv")

# Plot a histogram of the gene-level statistics
hist(my_geneList, main="Histogram of SBS1 Gene-Level Statistics", xlab="Statistic Value")

# Detect duplicated genes based on their values
# Note: This might not work as expected because `duplicated()` checks for duplicated indices, not values.
# If you're trying to find duplicated gene names (assuming names are in `names(my_geneList)`), you should use:
duplicated_genes <- names(my_geneList)[duplicated(names(my_geneList))]
head(duplicated_genes)

# Prepare MSigDB gene sets for Homo sapiens
all_gene_sets <- msigdbr(species = "Homo sapiens")

# Optional: Subset to specific categories if desired
# my_sets <- all_gene_sets %>% filter(gs_cat %in% c("H", "C2", "C3", "C5", "C7"))

# Keep all sets for this example
my_sets <- all_gene_sets

# Flatten to dataframe and ensure distinct gene-set pairs
my_sets_t2g <- my_sets %>% distinct(gs_name, entrez_gene) %>% as.data.frame()

# Determine the number of unique gene sets
num_sets <- length(unique(my_sets_t2g$gs_name))
cat("Number of unique gene sets:", num_sets, "\n")

# OPTIONAL: Filter the gene sets to remove any gene not part of your dataset
# This step assumes your geneList has entrez IDs as names
# sets_noMissing <- my_sets_t2g %>% filter(entrez_gene %in% names(my_geneList))

# Identify and retain gene sets with 10 to 1000 genes present in our data (adjust according to your needs)
# This part was commented out in your code, but here's how you could do it:
# filtered_gsNames <- sets_noMissing %>%
#   group_by(gs_name) %>%
#   summarize(genes_in_set = n()) %>%
#   filter(genes_in_set >= 10, genes_in_set <= 1000) %>%
#   select(gs_name)

# Proceed with sets_to_Use_t2g as per your preference (filtered or unfiltered)
sets_to_Use_t2g <- my_sets_t2g

# Stats on how many gene sets will be tested
gene_set_stats <- sets_to_Use_t2g %>%
  group_by(gs_name) %>%
  summarize(genes_in_set = n()) %>%
  filter(genes_in_set >= 10, genes_in_set <= 1000)

```



```{r}
# Checking for duplicates within each gene set
duplicated_in_sets <- sets_to_Use_t2g %>%
  group_by(gs_name) %>%
  filter(duplicated(entrez_gene) | duplicated(entrez_gene, fromLast = TRUE))

# Print a message depending on whether duplicates are present or not
if (nrow(duplicated_in_sets) > 0) {
  print("Duplicates found within gene sets.")
  # Optionally, print the duplicated entries
  print(duplicated_in_sets)
} else {
  print("No duplicates found within gene sets.")
}


```






```{r}

library(clusterProfiler)

# run GSEA
set.seed(4441) # You can choose any number for the seed and every time I want to reproduce the same results, set the seed to the same value before running the GSEA() function.
Results_gsea <- GSEA(	my_geneList, 
                      minGSSize = 10,
                      maxGSSize = 1000,
                      TERM2GENE=sets_to_Use_t2g, 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH",
                      by="fgsea",
                      eps=0)


```

```{r}


# Write Results_gsea to a CSV file

#write.csv(Results_gsea, file = "Results_2_gseaSBS1.csv", row.names = FALSE)
write.csv(Results_gsea, file = "Log_Results_CPTAC_gseaSBS1.csv", row.names = FALSE)

#Results_gsea <- read.csv("Results_CPTAC_gseaSBS1.csv")


# Make a summary of GSEA results
library(dplyr)  # Ensure dplyr is loaded for the pipe (%>%) and other functions
library(stringr)  # Ensure stringr is loaded for str_split()
library(readr)

my_gsea_summary <- head(Results_gsea, n=nrow(Results_gsea)) %>%
  rowwise() %>%
  mutate(core_size=length(str_split(core_enrichment, "\\/")[[1]])) %>%
  dplyr::select(-c(Description, leading_edge, core_enrichment)) %>%
  as.data.frame()
# print results summary to file
gsea_summ_filename <- paste("GSEA_MSigDB_summary.", "Genes_correlated_to_SBS1", ".01.txt", sep="")
write_tsv(my_gsea_summary, gsea_summ_filename, quote="none")

# Examine summary
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install clusterProfiler
BiocManager::install("clusterProfiler")

library(clusterProfiler)


my_gsea_summary %>% filter(NES < -3)
gseaplot(Results_gsea, geneSetID = "HALLMARK_GLYCOLYSIS", title = "HALLMARK_GLYCOLYSIS")

```
```{r}
#Exploring a Gene of Interest potential Biological function


library(clusterProfiler)

gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

library(clusterProfiler)
goplot(Results_gsea, geneSetID = "HALLMARK_GLYCOLYSIS", title = "HALLMARK_GLYCOLYSIS")



```
```{r}



nes_values <- Results_gsea$NES

# Calculate frequencies
num_directly_correlated <- sum(nes_values >= 2)
num_inversely_correlated <- sum(nes_values <= -2)
num_in_between <- sum(nes_values > -2 & nes_values < 2)

# Create a histogram of NES values, without coloring bars
hist_info <- hist(nes_values, plot=FALSE)

# Determine the colors for each bar
colors <- ifelse(hist_info$mids >= 2, "red", ifelse(hist_info$mids <= -2, "blue", "gray"))

# Plot the histogram with colors
plot(hist_info, col=colors, main="Histogram of SBS1 correlated CPTAC LSCC gene sets NES Values", xlab="NES", ylab="Frequency of gene sets", border="black")

# Calculate the mean and standard deviation of the NES values
mean_nes <- mean(nes_values)
std_dev_nes <- sd(nes_values)

# Calculate how many standard deviations away from the mean the threshold values are
directly_correlated_z_score <- (2 - mean_nes) / std_dev_nes
inversely_correlated_z_score <- (-2 - mean_nes) / std_dev_nes

# Add legend with calculated frequencies and Z-scores
legend(x = "topright", y = max(hist_info$counts) - 1, legend=c(
  paste("Directly correlated gene sets (NES >= 2): ", num_directly_correlated,
        " (Z-score ~", round(directly_correlated_z_score, 2), ")"),
  paste("Inversely correlated gene sets (NES <= -2): ", num_inversely_correlated,
        " (Z-score ~", round(inversely_correlated_z_score, 2), ")"),
  paste("No correlation gene sets (NES in between -2 and 2): ", num_in_between)),
  col=c("red", "blue", "gray"), pch=15, cex=0.6, box.col="white")

```

```{r}
# Read the CSV file into a dataframe
data45 <- read.csv("Log_Results_CPTAC_gseaSBS1.csv")

# Sort the data by the NES column in descending order
sorted_data <- data45[order(-data45$NES),]

# Optionally, print the sorted data to console
print(sorted_data)

# Write the sorted data to a new CSV file
write.csv(sorted_data, "Log_Sorted_CPTAC_gseaSBS1.csv", row.names = FALSE)


```








```{r}
# Data for heatmap for core enrichment genes, expression are already log2(TPM+1)
# need to transpose sinc this has genes as columns
my_exp <- as.data.frame(t(allGenes_allSamples))
normDataForPlot_medCenter <- my_exp - apply(my_exp, 1, median)
normDataForPlot_scaled <- t(apply(my_exp, 1, scale))
colnames(normDataForPlot_scaled) <- colnames(my_exp)
plotted_samples <- SBS_table[colnames(my_exp),]

library(RColorBrewer)
set.seed(4359)
long_color_vector <- c(brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 10, name = "RdYlBu"), brewer.pal(n = 10, name = "PuOr"))




legend_annotation <- 	list(	title = "Median-centered log2TPM", 
                            at = c(-4, -2, 0, 2, 4), 
                            labels = c(-4, -2, 0, 2, 4),
                            direction = "horizontal")















```

```{r} 
#source of ERRORS
# make report for sets of interest
soi <- c(
  
  
  #Upregulated
  
  
  "GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_DN",
  "chr19p13",
  "GSE22886_NAIVE_TCELL_VS_DC_UP",
  "GINESTIER_BREAST_CANCER_20Q13_AMPLIFICATION_DN",
  "GSE6269_HEALTHY_VS_FLU_INF_PBMC_UP",
  "chr17q11",
  "chr5q35",
  "TRAVAGLINI_LUNG_CD4_NAIVE_T_CELL",
  "JISON_SICKLE_CELL_DISEASE_DN",
  "NIKOLSKY_BREAST_CANCER_17Q11_Q21_AMPLICON",
  
  
  #Downregulated
  
  "GSE22886_NAIVE_BCELL_VS_NEUTROPHIL_DN",
  "TRAVAGLINI_LUNG_EREG_DENDRITIC_CELL",
  "WINTER_HYPOXIA_UP",
  "HOWARD_PBMC_INACT_MONOV_INFLUENZA_A_INDONESIA_05_2005_H5N1_AGE_19_39YO_AS03_ADJUVANT_VS_BUFFER_1DY_UP",
  "RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS",
  "HALLMARK_MTORC1_SIGNALING",
  "GSE6269_FLU_VS_STAPH_AUREUS_INF_PBMC_DN",
  "TRAVAGLINI_LUNG_MACROPHAGE_CELL",
  "TRAVAGLINI_LUNG_PROLIFERATING_MACROPHAGE_CELL",
  "TRAVAGLINI_LUNG_TREM2_DENDRITIC_CELL"
  
  
)

library(stringr)
library(clusterProfiler)
library(grid)
library (tibble)


list_running_plots <- list()
list_stat_boxes <- list()
list_heatmap_plots <- list()
list_set_description <- list()
gseaResToUse <- Results_gsea %>% filter(ID %in% soi)
contrast_name <- "SBS1_corr"

for (i in 1:length(soi)){
  my_soi <- gseaResToUse$ID[i]
  
  # running ES plot
  list_running_plots[[i]] <- gseaplot(gseaResToUse, geneSetID = i, title = paste(gseaResToUse$ID[i], contrast_name, sep=" --- "))
  
  # stat box
  my_stats=paste(	paste0( "Set ID: ", gseaResToUse$ID[i]),
                  paste0(	"Set size: ", gseaResToUse$setSize[i]),
                  paste0(	"Enrichment score: ", signif(gseaResToUse$enrichmentScore[i], digits=3)),   
                  paste0(	"Normalized enrichment score: ", signif(gseaResToUse$NES[i], digits=3)), 	 
                  paste0(	"p value: ", signif(gseaResToUse$pvalue[i], digits=3)),	
                  paste0(	"p value adjusted: ", signif(gseaResToUse$p.adjust[i], digits=3)),
                  paste0(	"q value: ", signif(gseaResToUse$qvalue[i], digits=3)),
                  paste0(	"Core enrichment genes: ", length(str_split(gseaResToUse$core_enrichment[i], "\\/")[[1]])),
                  sep="\n")
  
  #if(gseaResToUse$NES[i] >0){
  
  #if(!is.na(gseaResToUse$NES[i]) && gseaResToUse$NES[i] > 0){
  
  if(gseaResToUse$NES[i] >0){
    text_color <- "red"
  } else{
    text_color <- "blue"				
  }
  list_stat_boxes[[i]] <- textGrob(my_stats, just="left", gp=gpar(col=text_color, fontsize=14), x = unit(0, "npc"), y = unit(0.5, "npc"))
  
  # description
  set_description <- strwrap(as.character((all_gene_sets %>% filter(gs_name %in% gseaResToUse$ID[i]) %>%  dplyr::select(gs_description))[1,]), width = 100, simplify = FALSE)
  set_description_formatted <- sapply(set_description, paste, collapse = "\n")
  list_set_description[[i]] <- textGrob(set_description_formatted, just="left", gp=gpar(col="black", fontsize=12), x = unit(0, "npc"), y = unit(0.5, "npc"))
  
  # heatmap
  # get ENTREZID
  core_ID <- str_split(gseaResToUse$core_enrichment[i], "\\/")[[1]]
  # Get rownames for the genes
  core_rowname <- gene_info %>% filter(Entrez %in% core_ID) %>%  dplyr::select(ID)
  # Get expression data
  data_to_plot <- normDataForPlot_medCenter[core_rowname$ID, ]
  
  
  # Determine row labels (symbols)
  row_labels_to_label <- gene_info[core_rowname$ID, "Symbol"]
  labels_DF <- data.frame(Gene=row_labels_to_label, Term=rep(gseaResToUse$ID[i], times=length(row_labels_to_label)))
  write_tsv(labels_DF,
            paste("GSEA_MSigDB_summary.", contrast_name, ".core_genes.my_soi", ".txt", sep=""),
            quote="none")
  
  ha_row <- rowAnnotation(	Genes = anno_text(	row_labels_to_label, 
                                              just = "center", 
                                              location = unit(0.5, "npc"),
                                              gp = gpar(fontsize = 3)), 
                           annotation_name_rot = 0)
  
  hm_temp <-	Heatmap(data_to_plot, # col=color_scale, 
                     cluster_columns = TRUE, 
                     cluster_rows = TRUE,	# column_order = c(4:18, 1:3),
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     column_names_gp = gpar(fontsize = 8), 
                     width = unit(30, "cm"), 
                     height = unit(15, "cm"), 
                     heatmap_legend_param = legend_annotation, 
                     #top_annotation= column_ha, 
                     column_title = paste("Core set of genes", gseaResToUse$ID[i],  sep="\n"),
                     clustering_distance_rows="euclidean", 
                     clustering_method_rows="complete", 
                     right_annotation = ha_row,
                     use_raster = TRUE,
                     raster_device = "png",
                     raster_quality = 10,
                     raster_resize_mat = FALSE,
                     raster_by_magick = requireNamespace("magick", quietly = TRUE))
  
  list_heatmap_plots[[i]] <- draw(hm_temp, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  
  # Assuming list_running_plots[[i]] is your list of ggplot objects
  plots_list <- list_running_plots[[i]]
  
  # Combine all ggplot objects in the list into a single grob
  combined_plot <- plot_grid(plotlist = plots_list)
  
  # make plot and save
  montage_gsea <- multi_panel_figure(
    width = 18,
    height = 18,
    rows = 7,
    columns = 6,
    row_spacing = 0.2,
    column_spacing = 0.2,
    unit = "in",
    figure_name = "GSEA results summary",
    panel_label_type = "upper-alpha")
  
  montage_gsea %<>% 
    fill_panel(list_stat_boxes[[i]], row = 1, column = 1:3, label="A") %<>%  
    fill_panel(list_set_description[[i]], row = 2, column = 1:3, label="B") %<>% 
    fill_panel(combined_plot, row = 1:2, column = 4:6, label="C") %<>% 
    fill_panel(list_heatmap_plots[[i]], row = 3:7, column = 1:6, label="D")
  
  # Try adding a title
  header1 <- "GSEA results with complete MSigDB"
  header2 <- "Gene sorted by Pearson correlation between SBS1 exposure values and expression in 484 cancer samples"
  montage_gsea <- annotate_figure(montage_gsea, top = text_grob(paste(header1, header2, sep="\n"), color = "black", face = "bold", size = 16))
  

  
  # save as pdf
  ggsave(paste("GSEA_MSigDB_summary.", contrast_name, ".Running_ES+heatmap.", my_soi, ".pdf", sep=""), plot=montage_gsea, width=19, height=19, unit="in")
  
}					
```


```{r}
# Load the necessary library
library(dplyr)
```
```{r}
# Read the CSV file
data10 <- read.csv("Log_Sorted_CPTAC_gseaSBS1.csv")

# Filter the rows where NES values are >= 2 or <= -2
up_filtered_data <- subset(data10, NES >= 2)

# Optionally, write the filtered data to a new CSV file
write.csv(up_filtered_data, "Upregulated_Filtered_Results_gseaSBS1.csv", row.names = FALSE)
```
```{r}
# Read the CSV file
data10 <- read.csv("Log_Sorted_CPTAC_gseaSBS1.csv")

# Filter the rows where NES values are >= 2 or <= -2
dn_filtered_data <- subset(data10, NES <= -2)

# Optionally, write the filtered data to a new CSV file
write.csv(dn_filtered_data, "Downregulated_Filtered_Results_gseaSBS1.csv", row.names = FALSE)
```

#Compare Cluster
#data10 <- read.csv("Log_Sorted_CPTAC_gseaSBS1.csv")
#up_filtered_data <- subset(data10, NES >= 2)
#dn_filtered_data <- subset(data10, NES <= -2)
#formula_res <- compareCluster(up_filtered_data+dn_filtered_data, data=data10, fun="enrichKEGG")
#head(formula_res)


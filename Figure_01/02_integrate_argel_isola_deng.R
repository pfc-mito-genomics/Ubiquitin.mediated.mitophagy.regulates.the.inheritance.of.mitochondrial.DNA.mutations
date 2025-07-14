# Author: Cameron Ryall
# Project: MF_selection
# Description: This script merges and integrates Seurat objects from Isola et al. (2024), Argel et al., and Deng datasets. 
#              It includes normalization, PCA, Harmony integration, clustering, and plotting of marker gene expression (e.g., Zp3).
#              Outputs include intermediate and final Seurat objects and a series of UMAP/feature plots for downstream interpretation.
# Runtime: ~2 hours
library(MouseGastrulationData)
library(Seurat)
library(harmony)
library(glmGamPoi)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
library(plotly)

PROJECT = "MF_selection"
HOME <- "/suffolk/WorkGenomics/cdr42/MF_selection"
#################################################################
# Loading in Argel multiome from package: MouseGastrulationData #
#################################################################
argel_multiome <- RAMultiomeData(type = "rna", samples = NULL)

# Exclude CRISPR genotype (this dataset has WT and CRISPR knockout, only want WT)
meta_cleaned_argel_multiome <- argel_multiome[,colData(argel_multiome)$genotype == "WT"]

# Collect metadata from SCE
argel_meta_data <- as.data.frame(colData(meta_cleaned_argel_multiome))

# Convert to Seurat object and add in metadata
argel_seu <- CreateSeuratObject(counts(meta_cleaned_argel_multiome))
argel_seu <- AddMetaData(argel_seu, argel_meta_data)

# Define assay as counts
counts <- SeuratObject::GetAssayData(object = argel_seu, slot = "counts")
argel_seu <- SetAssayData(object = argel_seu, slot = "data", new.data = counts)

# Save to RDS
saveRDS(argel_seu, "argel/seurats/argel_wt_seurat.RDS")

# Add cell cycle scores to Argel dataset
argel_seu <- CellCycleScoring(argel_seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
Idents(argel_seu) <- "orig.ident"
gc()

######################################################################
# Load Isola et al 2024 Seurat object from Granulosa cell subset RDS #
######################################################################
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE232309&format=file&file=GSE232309%5FGranulosa%5FCells%2Erds%2Egz", "/suffolk/WorkGenomics/cdr42/MF_selection/isola/merged_seurats/GSE232309_Granulosa_Cells.rds")
R.utils::gunzip(destfile, remove = FALSE)  # remove=FALSE keeps original gz file

isola_seu=readRDS("/suffolk/WorkGenomics/cdr42/MF_selection/isola/merged_seurats/GSE232309_Granulosa_Cells.rds")

########################################################################
# Load Deng et al 2014, processed in 01_deng_prepare_expression_matrix #
########################################################################
deng_seurat <- readRDS("/suffolk/WorkGenomics/cdr42/MF_selection/deng/deng_seurat.RDS")


# Merge datasets into a single Seurat object
data.filt <- merge(x=isola_seu, y= c(argel_seu, deng_seurat), 
                   add.cell.ids = c("isola", "argel", "deng"), project="ovary",
                   merge.data = TRUE)


# Normalize, identify variable features, and scale data
data.filt <- data.filt %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 200000) %>% 
  ScaleData()

# Run PCA on merged dataset
data.filt <- RunPCA(data.filt, assay = "RNA", npcs = 50, reduction.name = "merged_pca")

# Extract sample ID from cell names and store in metadata
cell_names <- colnames(data.filt)
sample_meta <- substr(cell_names, 1, 5)
data.filt <- AddMetaData(data.filt, sample_meta, col.name = "sample")

# Combine potential celltype label sources into a unified label
data.filt@meta.data$celltype_combined <- coalesce(
  data.filt@meta.data$cluster.names, 
  data.filt@meta.data$celltype,
  data.filt@meta.data$cell_type)

# Save the merged, normalized (but not integrated) Seurat object
saveRDS(data.filt, "/suffolk/WorkGenomics/cdr42/MF_selection/deng/merged_seurats/merged.RDS")
gc()

# --------------------- Integration with Harmony --------------------- #

# Run Harmony integration based on sample identity
harmonized_seurat <- RunHarmony(data.filt, 
                                group.by.vars = "sample", 
                                reduction = "merged_pca", assay.use = "RNA", reduction.save = "harmony")


# Dimensionality reduction and clustering on integrated data
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:40, reduction.name = "harmonised.umap")
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 1, cluster.name = "harmonised.clusters")

# Save Harmony-integrated Seurat object
saveRDS(harmonized_seurat, "/suffolk/WorkGenomics/cdr42/MF_selection/deng/merged_seurats/harmonised.RDS")

# --------------------- Expression and QC Plots --------------------- #

# Set working directory for figure output
setwd(HOME)

# Plot Zp3 expression in Isola sample only
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/isola_zp3_expression.pdf", width = 6, height = 8)
FeaturePlot(isola_seu, features = "Zp3")
dev.off()

# Cleanup
rm(isola_seu, argel_seu)
gc()

# Load previously saved merged (but not integrated) object
merged_seu <- readRDS("/suffolk/WorkGenomics/cdr42/MF_selection/deng/merged_seurats/merged.RDS")
merged_seu <- RunUMAP(merged_seu, reduction = "merged_pca", assay = "RNA", dims = 1:40, reduction.name = "merged.umap")
merged_seu <- FindNeighbors(object = merged_seu, reduction = "merged_pca")
merged_seu <- FindClusters(merged_seu, resolution = 1, cluster.name = "merged.clusters")

# UMAP colored by age
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/merged_UMAP.pdf", width = 6, height = 6)
DimPlot(merged_seu, reduction = "merged_pca", group.by = "age", raster= FALSE)
dev.off()


# Violin plots for Zp gene expression per sample
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/merged_sample_specific_zp3_violin.pdf", width = 6, height = 6)
VlnPlot(
  merged_seu,
  group.by="sample",
  features = c("Zp1", "Zp2", "Zp3"),
  ncol = 3, pt.size = 0.1, raster = FALSE)
dev.off()

# FeaturePlot of Zp3
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/merged_zp3_expression.pdf", width = 6, height = 6)
FeaturePlot(merged_seu, features = "Zp3", raster=FALSE, reduction = "merged.umap")
dev.off()

rm(merged_seu)
gc()

# --------------------- Post-integration plots --------------------- #
# Unify age/stage into a single metadata column
data.filt@meta.data$age_combined <- ifelse(!is.na(data.filt@meta.data$stage),
                                                       data.filt@meta.data$stage,
                                                       data.filt@meta.data$age)

# UMAPs by age, sample, and cell type
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/integrated_age_heatmap.pdf", width = 6, height = 6)
DimPlot(data.filt, reduction = "harmonised.umap", group.by = "age_combined", raster= FALSE)
dev.off()

pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/integrated_sample_heatmap.pdf", width = 6, height = 6)
DimPlot(data.filt, reduction = "harmonised.umap", group.by = "sample", raster= FALSE)
dev.off()

pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/integrated_sample_heatmap.pdf", width = 6, height = 6)
DimPlot(data.filt, reduction = "harmonised.umap", group.by = "celltype_combined", raster= FALSE)
dev.off()

# Zp3 expression after integration
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/integrated_zp3_expression.pdf", width = 6, height = 6)
FeaturePlot(data.filt, reduction = "harmonised.umap", features = "Zp3", raster= FALSE)
dev.off()

# Violin plots of Zp genes by Harmony clusters
pdf("/suffolk/WorkGenomics/cdr42/MF_selection/figures/integration_figures/deng/integrated_cluster-specific_zp1-2-3_violin.pdf", width = 6, height = 6)
VlnPlot(
  data.filt,
  group.by="harmonised.clusters",
  features = c("Zp1", "Zp2", "Zp3"),
  ncol = 3, pt.size = 0.1, raster = FALSE)
dev.off()
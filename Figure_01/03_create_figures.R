# Supplementary Script for Mitochondrial Pathway Analysis in Germline scRNA-seq Data

# Load required libraries
library(SingleCellExperiment)
library(scran)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scater)
library(readr)
library(gplots)
library(viridis)
library(Seurat)
library(FSA) 
library(rstatix)
library(presto)
library(plotly)
library(ComplexHeatmap)
library(biomaRt)
library(ggpubr)
library(circlize)

# Set working directory (change as needed)
setwd("~/Mount/Suffolk/WorkGenomics/cdr42/MF_selection/")

# Load harmonised Seurat object
young_integrated <- readRDS("deng/merged_seurats/harmonised.RDS")

# Combine ages into one metadata column
young_integrated@meta.data$age_combined <- ifelse(!is.na(young_integrated@meta.data$stage),
                                                       young_integrated@meta.data$stage,
                                                       young_integrated@meta.data$age)

# UMAP plot of cell types (interactive)
p <- DimPlot(young_integrated, reduction = "harmonised.umap", group.by = "celltype_combined", raster = FALSE)
ggplotly(p)

# Load curated mitochondrial function gene list
curated_MF_mito_genes_df <- read_csv("MF_input_genes/MF_KEGG_genes.csv")
curated_MF_mito_genes_df <- curated_MF_mito_genes_df[!duplicated(curated_MF_mito_genes_df), ]
curated_MF_mito_genes <- curated_MF_mito_genes_df$GeneName
table(curated_MF_mito_genes_df$KEGG_Pathway_Name)

# Filter to genes present in dataset
filtered_genes <- rownames(young_integrated$RNA)[rownames(young_integrated$RNA) %in% curated_MF_mito_genes_df$GeneName]
curated_MF_mito_genes_df <- curated_MF_mito_genes_df[curated_MF_mito_genes_df$GeneName %in% filtered_genes, ]

# Subset for germline and PGC/oocyte cells, remove zygotes
deng_indices <- which(young_integrated@meta.data$sample == "deng_" & 
                        young_integrated@meta.data$celltype_combined != 'zy1' &
                        young_integrated@meta.data$celltype_combined != 'zy2' &
                        young_integrated@meta.data$celltype_combined != 'zy3' &
                        young_integrated@meta.data$celltype_combined != 'zy4')

set.seed(123)  # Setting a seed for reproducibility

# Create a logical vector for other cell types (PGC, Oocytes, Luteal) and combine
other_indices <- which(young_integrated@meta.data$celltype_combined %in% c("PGC", "Oocytes"))
germline_indices <- c(deng_indices, other_indices)

# Subset the Seurat object with the selected indices
germline_seu <- young_integrated[, germline_indices]


# Remove Seurat cluster 33 oocytes
filtered_integrated <- subset(germline_seu, !(Idents(germline_seu) == 33 & celltype_combined == "Oocytes"))

# Order cell types
cell_order <- c("PGC", "Oocytes","early2cell", "mid2cell", "late2cell", 
                "4cell", "8cell", "16cell", "earlyblast", 
                "midblast", "lateblast")
filtered_integrated@meta.data$celltype_combined <- factor(filtered_integrated@meta.data$celltype_combined, levels = cell_order)

# Prepare scaled expression matrix
mat <- filtered_integrated[filtered_genes, ]@assays$RNA$scale.data %>% as.matrix()

# Annotations
cluster_anno<- filtered_integrated@meta.data$celltype_combined[complete.cases(mat)]
quantile(mat, c(0.1, 0.95), na.rm = TRUE)
col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
default_width <- unit(1, "cm")
pgc_width <- unit(3, "cm")
valid_indices <- complete.cases(mat)
cluster_anno <- filtered_integrated@meta.data$celltype_combined[valid_indices]

# Convert annotations to match the dimensions of `mat`
phase <- filtered_integrated@meta.data$Phase[valid_indices]
sample <- filtered_integrated@meta.data$sample[valid_indices]
age_combined <- filtered_integrated@meta.data$age_combined[valid_indices]
size_factor <- filtered_integrated@meta.data$sizeFactor[valid_indices]
seurat_clusters <- filtered_integrated@meta.data$seurat_clusters[valid_indices]

# Heatmap annotations
heatmap_annotations <- HeatmapAnnotation(
  Phase = phase,
  Sample = sample,
  Age = age_combined,
  Size_Factor = size_factor,
  Clusters = seurat_clusters,
  foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))  # You can keep your original annotation here
)

# KEGG Pathway annotation
pathway_order <- c("Ubiquitin-mediated Proteolysis", "Mitophagy", "Proteasome",
                   "Autophagy")
pathway_factors <- factor(curated_MF_mito_genes_df$KEGG_Pathway_Name[match(rownames(mat), curated_MF_mito_genes_df$GeneName)],
                          levels = pathway_order)
# Ensure the pathway factors are set as a factor with the desired levels
pathway_factors <- factor(pathway_factors, levels = unique(pathway_order))

# Create pathway annotation for the rows
pathway_annotation <- rowAnnotation(
  KEGG_Pathway = pathway_factors    # Set the font size for the text
)


# ComplexHeatmap plot
pdf(file = "figures/single_cell/all_pathways.pdf", width = 15, height = 15)
Heatmap(mat, name = "Expression",  
        column_split = filtered_integrated@meta.data$celltype_combined, 
        cluster_columns = FALSE,
        cluster_column_slices = TRUE, 
        row_split = pathway_factors,
        cluster_rows = TRUE,
        col = circlize::colorRamp2(c(-3, 0, 3), c("#FF00FF", "black", "#FFFF00")),
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        show_column_names = FALSE,
        show_row_names = FALSE,
       left_annotation = pathway_annotation
)
dev.off()

# Seurat DoHeatmap version (optional)
my_palette <- c("white", viridis(n = 255, option = "D"))
pdf(file = "figures/single_cell/all_pathways.pdf", width = 15, height = 15)
DoHeatmap(filtered_integrated, filtered_genes, group.by = "celltype_combined",  assay = "RNA", slot = "scale.data")
dev.off()

# Calculate average expression per cell per pathway
pathway_genes <- curated_MF_mito_genes_df %>%
  filter(GeneName %in% rownames(mat)) %>%
  group_by(KEGG_Pathway_Name) %>%
  summarise(genes = list(GeneName))
pathway_activation <- list()

for(i in seq_along(pathway_genes$KEGG_Pathway_Name)) {
  pathway <- pathway_genes$KEGG_Pathway_Name[i]
  genes_in_pathway <- pathway_genes$genes[[i]]
  sub_mat <- mat[genes_in_pathway, ]
  avg_expression <- colMeans(sub_mat, na.rm = TRUE)
  pathway_activation[[i]] <- data.frame(
  Cell = colnames(mat),
  CellType = cluster_anno,
  Pathway = pathway,
  AvgExpression = avg_expression
  )
}
pathway_activation_df <- do.call(rbind, pathway_activation)
pathway_activation_df <- pathway_activation_df[pathway_activation_df$Pathway %in% pathway_order, ]

# Plot pathway activation per cell type
ggplot(pathway_activation_df, aes(x = CellType, y = AvgExpression, color = Pathway)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.2) + 
  geom_boxplot(aes(fill = Pathway), outlier.shape = NA, alpha = 0.3, position = position_dodge(width = 0.75)) + 
  labs(x = "Cell Type", y = "Average Pathway Activation", title = "Pathway Activation Across Cell Types") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = scales::hue_pal()(length(unique(pathway_activation_df$Pathway)))) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(pathway_activation_df$Pathway))))

# Optional: Custom function to plot individual pathway heatmaps
heatmapper <- function(pathway, font){
  genes <- curated_MF_mito_genes_df$GeneName[curated_MF_mito_genes_df$KEGG_Pathway_Name == pathway]
  genes_in_mat <- genes[genes %in% rownames(mat)]
  
  if (length(genes_in_mat) == 0) {
    message("No genes found in the matrix for pathway: ", pathway)
    return(NULL)
  }
  
  print(genes_in_mat)  
  
  # Generate heatmap only for the genes that are found in mat
  Heatmap(mat[genes_in_mat, ], name = "Expression",  
          column_split = factor(cluster_anno),
          cluster_columns = TRUE,
          show_column_dend = FALSE,
          cluster_column_slices = TRUE,
          column_title_gp = gpar(fontsize = 8),
          column_gap = unit(0.5, "mm"),
          cluster_rows = TRUE,
          show_row_dend = FALSE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 4),
          column_title_rot = 90,
          top_annotation = heatmap_annotations,  # Use the new annotations
          show_column_names = FALSE)
}

pdf(file = "figures/single_cell/autophagy.pdf", width = 15, height = 15)
heatmapper("Autophagy factors", 0.4)
dev.off()
gc()

pdf(file = "figures/single_cell/mitophagy.pdf", width = 15, height = 15)
heatmapper("Mitophagy", 0.4)
dev.off()
gc()

pdf(file = "figures/single_cell/proteasome.pdf", width = 15, height = 15)
heatmapper("Proteasome", 0.4)
dev.off()
gc()

pdf(file = "figures/single_cell/replication.pdf", width = 15, height = 15)
heatmapper("Mitochondrial DNA replication factors", 0.6)
dev.off()
gc()

pdf(file = "figures/single_cell/dynamics.pdf", width = 15, height = 15)
heatmapper("Mitochondrial dynamics", 0.6)
dev.off()
gc()

pdf(file = "figures/single_cell/ubiquitin.pdf", width = 15, height = 15)
heatmapper("Ubiquitin-mediated Proteolysis", 0.4)
dev.off()
gc()

# Visualize Usp30 expression across all cells using a feature plot
FeaturePlot(filtered_integrated, features = "Usp30", reduction = "harmonised.umap", cols = c("lightgrey", "blue"))

# Look at subsets of oocytes, identified on heatmap
oocytes_subset <- subset(filtered_integrated, celltype_combined == "Oocytes")



# Plot USP30 (scaled data)
scaled_expression <- GetAssayData(filtered_integrated, assay = "RNA", slot = "scale.data")
usp30_expression <- scaled_expression["Usp30", ]
filtered_integrated@meta.data$USP30 <- usp30_expression
filtered_integrated@meta.data$CellOrder <- factor(filtered_integrated@meta.data$celltype_combined, levels = cell_order)

# Create dataframe for single-cell barplot
plot_data <- data.frame(
  CellID = rownames(filtered_integrated@meta.data),
  USP30 = filtered_integrated@meta.data$USP30,
  CellOrder = filtered_integrated@meta.data$CellOrder,
  Celltype = filtered_integrated@meta.data$celltype_combined
)

# Sort data by USP30 expression within ordered cell types
plot_data <- plot_data[order(plot_data$USP30, plot_data$CellOrder), ]

# Plot scaled expression
pdf(file = "figures/pseudo/usp30_germline_single_cell.pdf", width = 15, height = 15)
ggplot(plot_data, aes(x = reorder(CellID, USP30), y = USP30, fill = Celltype)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for better readability
  labs(title = "USP30 Expression in Germline Cells (Normalised)",
       x = "Cell ranked by Usp30 expression",
       y = "USP30 Expression") +
  theme_classic() + theme(axis.text.y = element_blank())
dev.off()

# Pseudobulk Boxplot and Statistics (scaled data)

# Prepare data
expression_data <- data.frame(
  CellType = filtered_integrated@meta.data$celltype_combined,
  USP30 = filtered_integrated@meta.data$USP30
)

expression_data <- expression_data %>%
  filter(!is.na(CellType) & !is.na(USP30))
expression_data$CellType <- factor(expression_data$CellType, levels = cell_order)

# Plot pseudobulk (scaled)
pdf(file = "figures/pseudo/usp30_germline_pseudo.pdf", width = 15, height = 15)
ggplot(expression_data, aes(x = CellType, y = USP30)) +
#  geom_violin(fill = "steelblue", position = position_nudge(-0.3), alpha = 0.5) +
  geom_boxplot(fill = "steelblue", width = 0.3, position = position_nudge(-0.3), outliers = FALSE) +
  geom_jitter(, size = 1, width = 0.05, height = 0, alpha=0.4) +
  labs(x = "Cell Type",
       y = "USP30 Expression") +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",  text = element_text(face = "bold"),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.15, colour = "red")  +
  labs(x = "Cell Type", y = "Normalised Usp30 Expression",
       title = "Usp30 expression - early development")

dev.off()

write.csv(dunn_test_result, "figures/pseudo/dunn_test_Usp30_germline.csv")


# Calculate RPKM values


# Get gene lengths using biomaRt
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_lengths_df <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'transcript_length'),
                         mart = mart)
gene_lengths_df$transcript_length_kb <- gene_lengths_df$transcript_length / 1000

# Define function to calculate RPKM
find_rpkm <- function(seurat, gene_lengths_df){
  # Extract counts for each layer
  counts_1 <- GetAssayData(seurat, assay = "RNA", slot = "counts.1")
  counts_2 <- GetAssayData(seurat, assay = "RNA", slot = "counts.2")
  
  # Ensure gene identifiers are consistent across layers
  common_genes <- Reduce(intersect, list(rownames(counts_1), rownames(counts_2)))
  
  # Subset counts matrices to only include common genes
  counts_1 <- counts_1[common_genes, ]
  counts_2 <- counts_2[common_genes, ]
  
  # Combine counts into one matrix
  combined_counts <- cbind(counts_1, counts_2)
  
  # Ensure gene symbols or identifiers match
  gene_symbols <- rownames(combined_counts)
  
  # Filter gene_lengths_df to have only genes present in the counts
  gene_lengths_df <- gene_lengths_df %>%
    filter(external_gene_name %in% gene_symbols)
  
  # Calculate total mapped reads per cell (library size)
  total_mapped_reads <- colSums(combined_counts)
  
  # Convert library size to millions
  total_mapped_reads_millions <- total_mapped_reads / 1e6
  
  # Match gene lengths with counts using gene symbols
  gene_lengths <- gene_lengths_df$transcript_length_kb[match(gene_symbols, gene_lengths_df$external_gene_name)]
  
  # Check for missing gene lengths
  if(any(is.na(gene_lengths))){
    warning("Some gene lengths are NA")
  }
  
  # Divide counts by gene lengths (in kilobases)
  rpk <- sweep(combined_counts, 1, gene_lengths, "/")
  
  # Divide by total mapped reads (in millions) for each sample
  rpkm <- sweep(rpk, 2, total_mapped_reads_millions, "/")
  
  # If you want to create separate assays for each count layer:
  seurat[["RPKM_counts_1"]] <- CreateAssayObject(counts = rpk[, 1:ncol(counts_1)])
  seurat[["RPKM_counts_2"]] <- CreateAssayObject(counts = rpk[, (ncol(counts_1) + 1):(ncol(counts_1) + ncol(counts_2))])
  
  # If you prefer to combine all RPKM values into a single assay:
  seurat[["RPKM_combined"]] <- CreateAssayObject(counts = rpk)
  
  return(seurat)
}

# Apply function
seu <- find_rpkm(filtered_integrated, gene_lengths_df)
print(seu[["RPKM_combined"]])

# Plot USP30 (RPKM)

# Extract RPKM and update metadata
usp30_expression <- FetchData(seu, vars = "Usp30", assay = "RPKM")
seu@meta.data$USP30 <- usp30_expression$Usp30
seu@meta.data$CellOrder <- factor(seu@meta.data$celltype_combined, levels = cell_order)

# Create dataframe for barplot
plot_data <- data.frame(
  CellID = rownames(seu@meta.data),
  USP30 = seu@meta.data$USP30,
  CellOrder = seu@meta.data$CellOrder,
  Celltype = seu@meta.data$celltype_combined
)

# Plot single-cell barplot (RPKM)
plot_data <- plot_data[order(plot_data$USP30, plot_data$CellOrder), ]
pdf(file = "figures/pseudo/usp30_RPKM_germline_single_cell.pdf", width = 15, height = 15)
ggplot(plot_data, aes(x = reorder(CellID, USP30), y = USP30, fill = Celltype)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "USP30 Expression in Germline Cells (RPKM)",
       x = "Cell ranked by Usp30 expression",
       y = "USP30 Expression (RPKM)") +
  theme_classic() + theme(axis.text.y = element_blank())


# Prepare data for pseudobulk plot (RPKM)
expression_data <- data.frame(
  CellType = seu@meta.data$celltype_combined,
  USP30 = seu@meta.data$USP30
)
pseudobulk_data <- expression_data %>%
  group_by(CellType) %>%
  summarise(
    MeanUSP30 = mean(USP30, na.rm = TRUE),
    SDUSP30 = sd(USP30, na.rm = TRUE),  # Standard deviation
    CellCount = n()  # Number of cells
  ) %>%
  mutate(CellType = factor(CellType, levels = cell_order))  # Ensure the correct order

# Plot pseudobulk (RPKM)
pdf(file = "figures/pseudo/usp30_RPKM_germline_pseudo.pdf", width = 15, height = 15)
ggplot(expression_data, aes(x = CellType, y = USP30)) +
  #geom_violin(fill = "steelblue", position = position_nudge(-0.3), alpha = 0.5) +
  geom_boxplot(fill = "steelblue", width = 0.3, position = position_nudge(-0.3), outliers = FALSE) +
  geom_jitter(, size = 1, width = 0.05, height = 0, alpha=0.4) +
  labs(x = "Cell Type",
       y = "USP30 RPKM") +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",  text = element_text(face = "bold"),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.15, colour = "red")  +
  labs(x = "Cell Type", y = "RPKM -  Usp30 Expression",
       title = "Usp30 expression - early development")

dev.off()

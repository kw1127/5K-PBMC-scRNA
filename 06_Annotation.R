# ---- 06. Cell annotation -----------------------------------------------------
# 
# Aim:
#   - Visualise expression of HVGs on the UMAP plot
#   - Identify differentially expressed genes (DEGs) defining cell types
#   - Visualise distribution of gene expression across cell types
#
# ------------------------------------------------------------------------------

# ---- Set optimal resolution ----
Idents(seurat) <- "RNA_snn_res.0.5"

# ---- Identify DEGs ----
positive.markers <- FindAllMarkers(
  seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# ---- Find top DEGs by cluster ----
top6 <- positive.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_max(avg_log2FC, n = 6) %>%   # best marker per cluster
  ungroup()
  
# ---- Visualise top 6 markers ----
DotPlot(seurat, features = unique(top6$gene)) +
  RotatedAxis()

DoHeatmap(seurat, features = top6$gene) + 
  NoLegend()

# ---- Identify top DEGs in each cluster ----
top_per_cluster <- positive.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_max(avg_log2FC, n = 1) %>%   # best marker per cluster
  ungroup() %>%
  pull(gene)

# ---- Visualise top DEGs ----
FeaturePlot(seurat, 
            features = top_per_cluster, 
            order = TRUE,
            min.cutoff = "q10",
            max.cutoff = "q90"
)

# ---- List canonical PBMC markers ----
markers <- list(
  "CD4 T" = c("CD3D", "CD4", "IL7R", "CCR7"),
  "CD8 T" = c("CD8A", "CD8B"),
  "NK" = c("GNLY", "NKG7", "KLRD1"),
  "B" = c("MS4A1", "CD79A", "CD19"),
  "CD14 Monocyte" = c("CD14", "LYZ", "S100A8"),
  "CD16 Monocyte" = c("FCGR3A", "MS4A7"),
  "DC" = c("FCER1A", "CST3"),
  "Platelet" = c("PPBP", "PF4")
)

DotPlot(
  seurat, 
  features = markers) + 
  RotatedAxis()

# ---- T-cell subtyping ----
t_markers <- c(
  # Pan-T
  "CD3D", "CD3E",
  # CD4
  "CD4", "IL7R", "CCR7", "LEF1",
  # CD8
  "CD8A", "CD8B", "CD8B2",
  # Cytotoxic / effector
  "GZMK", "GZMB", "GZMA", "PRF1", "NKG7", "GNLY",
  # Memory/naive
  "SELL", "TCF7", "S100A4",
  # Treg
  "FOXP3", "IL2RA",
  # MAIT
  "KLRB1", "SLC4A10"
)

DotPlot(seurat, features = t_markers,
        idents = c("0", "1", "2", "5", "8")) + RotatedAxis()

# ---- Final cluster annotation ----
cluster.ids <- c(
  "0"  = "CD4 T naive",
  "1"  = "CD4 T memory",
  "2"  = "CD4 T",
  "3"  = "CD14 Monocyte",
  "4"  = "CD14 Monocyte",
  "5"  = "NK",
  "6"  = "B",
  "7"  = "B",
  "8"  = "CD4 T naive",
  "9"  = "CD16 Monocyte",
  "10" = "DC",
  "11" = "Platelet"
)

seurat <- RenameIdents(seurat, cluster.ids)
seurat$cluster_celltype <- Idents(seurat)

# ---- Final UMAP ----
DimPlot(
  seurat,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.5
) + 
  NoLegend()

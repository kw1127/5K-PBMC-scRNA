# ---- 06. Cell annotation ----------------------------------------------------------
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

# ---- Filter and sort DEGs by cluster association ----
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

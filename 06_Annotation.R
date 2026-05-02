# ---- 06. Cell annotation ----------------------------------------------------------
# 
# Aim:
#   - Visualise expression of HVGs on the UMAP plot
#   - Identify differentially expressed genes (DEGs) defining cell types
#   - Visualise distribution of gene expression across cell types
#
# ------------------------------------------------------------------------------

# ---- Plot 6 top HVGs ----
top10features <- FeaturePlot(
  seurat,
  features = top10[1:6],
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
)
top10features

# ---- Find DEGs defining cell types ----
positive.markers <- FindAllMarkers(
  seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# ---- Filter and sort DEGs by cluster association ----
positive.markers <- positive.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1)


# ---- Visualise distribution of some genes ----
VlnPlot(seurat, features = c("IL32", "IL7R", "GLA"))

# ---- 05. Clustering ----------------------------------------------------------
# 
# Aim:
#   - Find k.param neighbours using FindNeighbours()
#   - Compute clusters using FindClusters() using different methods and resolutions
#   - Run UMAP to visualise clusters
#
# ------------------------------------------------------------------------------

# ---- Set seed ----
set.seed(42)

# ---- Find KNN ----
seurat <- FindNeighbors(seurat, dims = 1:optimal_pc, k.param = 20)

# ---- Compute clusters ----
seurat <- FindClusters(seurat, resolution = c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0))

# ---- UMAP ----
seurat <- RunUMAP(seurat, dims = 1:optimal_pc)

# ---- Visualise ----
DimPlot(seurat, reduction = "umap", label = TRUE)


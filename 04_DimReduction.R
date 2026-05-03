# ---- 04. Dimension reduction -------------------------------------------------
# 
# Aim:
#   - Run Principle Component Analysis (PCA) using HVGs
#   - Identify optimal number of PCs using an elbow plot
# 
# ------------------------------------------------------------------------------

# ---- Run PCA ----
seurat <- RunPCA(seurat, npcs = 50)

# ---- Compute suggested PC cutoff ----
pct  <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100 # Calculate % variation explained by each PC
cumu <- cumsum(pct)

# Heuristic 1: cumulative variance > 90% + individual PC < 5%
co1 <- which(cumu > 90 & pct < 5)[1]

# Heuristic 2: where the difference between consecutive PCs drops below 0.1%
co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
            decreasing = TRUE)[1] + 1

optimal_pc <- min(co1, co2)

# ---- Elbow plot ----
ElbowPlot(seurat, ndims = 50, reduction = "pca") +
  geom_vline(xintercept = optimal_pc,
             linetype = "dashed", colour = "blue") +
  annotate("text",
           x = optimal_pc + 1, y = max(seurat[["pca"]]@stdev),
           label = paste0("Suggested: PC", optimal_pc),
           hjust = 0, colour = "blue", size = 4) +
  labs(title = "PCA Elbow Plot",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ---- Visualise genes defining PC1 and PC2 ----
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

# ---- Visualise the PCA output on a 2D scatter plot ----
DimPlot(seurat, reduction = "pca") + NoLegend()
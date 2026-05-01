# ---- 03. Normalisation -------------------------------------------------
# 
# Aim:
#   - Log-transform raw counts using NormalizeData (scale.factor = 10,000)
#   - Identify highly variable genes (HVGs) using vst (nfeatures = 2000)
#   - Scale the log-transformed data using ScaleData
# 
# ------------------------------------------------------------------------------

# ---- Normalise ----
seurat <- NormalizeData(seurat)

# ---- Find HVGs ----
seurat <- FindVariableFeatures(
  seurat,
  selection.method = "vst",
  nfeatures = 2000
)

# ---- Identify top 10 HVGs ----
top10 <- head(VariableFeatures(seurat), 10)

# ---- Plot HVGs unlabeled ----
hvg <- VariableFeaturePlot(seurat)

# ---- Plot HVGs with the top 10 label ----
hvg_labelled <- LabelPoints(
  plot = hvg,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

# ---- Plot side-by-side ----
hvg / hvg_labelled

# ---- Scale data ----
seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

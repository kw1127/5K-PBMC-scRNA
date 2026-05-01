# ---- 02. Quality Control (QC) -------------------------------------------------
# 
# Aim:
#   - Compute QC metrics (percent.mt, percent.ribo, percent.hb)
#   - Visualise the distribution in the data
#   - Study the relationship between QC metrics
#   - Subset the data for downstream analysis.
# 
# ------------------------------------------------------------------------------

# ---- QC ----
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")

head(seurat@meta.data)

# ---- Visualise distribution ----
VlnPlot(seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 4, 
        pt.size = 0.1)

qc1 <- seurat@meta.data %>%
  ggplot(aes(nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "skyblue") +
  geom_vline(xintercept = c(200, 5000), linetype = "dashed", color = "red") +
  scale_x_log10() +
  theme_minimal()

qc2 <- seurat@meta.data %>%
  mutate(ident = Idents(seurat)) %>%
  ggplot(aes(ident, percent.mt)) +
  geom_violin() +
  geom_hline(yintercept = c(5, 15), linetype = "dashed", colour = "red") +
  theme_minimal()

qc3 <- seurat@meta.data %>%
  mutate(ident = Idents(seurat)) %>%
  ggplot(aes(ident, percent.ribo)) +
  geom_violin() +
  geom_hline(yintercept = c(5, 15), linetype = "dashed", colour = "red") +
  theme_minimal()

qc1 | qc2 | qc3

# ---- Visualise the relationship ----
p1 <- FeatureScatter(seurat, feature1 = "nCount_RNA",   feature2 = "percent.mt")
p2 <- FeatureScatter(seurat, feature1 = "nCount_RNA",   feature2 = "nFeature_RNA")
p3 <- FeatureScatter(seurat, feature1 = "percent.mt",   feature2 = "percent.ribo")
p4 <- FeatureScatter(seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt")

(p1 | p2) / (p3 | p4)

# ---- Filter cells based off QC visualisations ----
seurat <- subset(
  seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15
)

# ----01. Loading the data and creating the Seurat object ---------------------------------------------
# 
# Aim:
#   - Load packages for scRNA-seq analysis
#   - Reads the 10X filtered matrix for the 5k healthy-donor PBMC dataset.
#   - Build a Seurat object and tag each cell will QC metrics.
#   - Save Seurat object as an RDS file for future work.
#   
#
# Data: 10x Genomics public dataset | single healthy donor
# =============================================================================

# ---- Load packages ----
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# ---- Read the counts matrix ----
counts <- Read10X(data.dir = "Counts_matrix_PBMC")

seurat <- CreateSeuratObject(
  counts = counts,
  project = "5k_pbmc",
  min.cells = 3,
  min.features = 200
)

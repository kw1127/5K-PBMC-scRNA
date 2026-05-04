# ---- 07. Footprint analysis --------------------------------------------------
# 
# Aim:
#   - Use the PROGENy network to analyse pathway activity
#   - Use the CollecTRI network to analyse TF activity
# ------------------------------------------------------------------------------

# ---- Libraries ----
library(decoupleR)
library(pheatmap)

# ---- Extract PROGENy network ----
progeny <- get_progeny(organism = "human", top = 500)

# ---- Matrix of normalised counts ----
matrix <- as.matrix(seurat[["RNA"]]$data)

# ---- Univariate linear model ----
progeny_ulm <- run_ulm(
  mat = matrix,
  network = progeny,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 5
)

# ---- Stash as new assay ----
seurat[["pathwaysulm"]] <- progeny_ulm %>%
  tidyr::pivot_wider(id_cols = "source",
                     names_from = "condition",
                     values_from = "score") %>%
  tibble::column_to_rownames(var = "source") %>%
  Seurat::CreateAssayObject()

DefaultAssay(seurat) <- "pathwaysulm"

FeaturePlot(
  seurat, 
  features = c("TNFa", "NFkB", "JAK-STAT", "MAPK")
)

# ---- Further exploration ----


# ---- 07. Differential expression (DE) ----------------------------------------
# 
# Aim:
#   - Find DE features between cell types using FindMarkers() 
# ------------------------------------------------------------------------------

# ---- Libraries ----
library(msigdbr)
library(clusterProfiler)

# ---- Find DE features between NK and CD4+ T cells ----
de.features <- FindMarkers(
  seurat,
  ident.1 = "NK",
  ident.2 = "CD4 T"
)

# ---- Find positive DE features between NK cells and all other cell types ----
de.pos.features <- FindMarkers(
  seurat,
  ident.1 = "NK",
  ident.2 = NULL,
  only.pos = TRUE
)

# ---- Extract hallmark gene sets ----
hallmark <- msigdbr(db_species = "HS", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# ---- Over representation analysis (ORA) on NK vs CD4 T DEGs ----
nk_genes <- rownames(de.features)[de.features$avg_log2FC > 1 &
                                  de.features$p_val_adj < 0.05]

ora_nk <- enricher(nk_genes, TERM2GENE = hallmark)

# ---- Visualise enrichment results ----
dotplot_nk <- clusterProfiler::dotplot(
  ora_nk,
  showCategory = 20
)

dotplot_nk

cnet_nk <- clusterProfiler::cnetplot(
  ora_nk,
  showCategory = 5
)

cnet_nk

heatplot_nk <- clusterProfiler::heatplot(
  ora_nk, 
  showCategory = 10
)

heatplot_nk

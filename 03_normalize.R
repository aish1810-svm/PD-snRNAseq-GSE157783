# ============================================================
# 03_normalize.R
# Purpose: Normalize counts, identify variable features, and scale
#          data in preparation for PCA/clustering.
#
# Note: the original authors already normalized this data as part of
# their own pipeline (confirmed via vst.* columns in genes.tsv), but
# only raw counts were included in the GEO supplementary files we
# have. We therefore normalize ourselves here, using standard Seurat
# defaults, to produce a self-consistent pipeline.
#
# Input:  data/processed/seurat_raw.rds
# Output: data/processed/seurat_normalized.rds
# ============================================================

library(Seurat)

seurat_obj <- readRDS("data/processed/seurat_raw.rds")

# ---- 1. Log-normalize counts ----
# Standard Seurat default: scales each cell's counts to a common
# total (10,000), then log1p-transforms. Corrects for differing total
# RNA capture per cell/nucleus, analogous to library-size
# normalization in bulk RNA-seq.
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ---- 2. Identify variable features ----
# Focuses downstream analysis (PCA, clustering) on the 2000 genes
# that vary most across cells, reducing noise from genes that are
# uniformly on/off/noisy across all cell types.
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# ---- 3. Scale data ----
# Centers each gene to mean 0, variance 1 across cells, so highly
# expressed genes don't dominate PCA purely due to scale.
seurat_obj <- ScaleData(seurat_obj)

# ---- Verification ----
head(VariableFeatures(seurat_obj))
seurat_obj

# ---- Save ----
saveRDS(seurat_obj, "data/processed/seurat_normalized.rds")
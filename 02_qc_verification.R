# ============================================================
# 02_qc_verification.R
# Purpose: Verify data quality of GSE157783 snRNA-seq data.
#
# IMPORTANT CONTEXT (confirmed via the original authors' GitHub repo,
# SpielmannLab/pd_human_midbrain_snrnaseq):
#   - CellRanger 3.0 was used for initial transcript counting
#   - Scrublet was already used to identify and remove doublets
#   - Seurat3/monocle3 were used for normalization, clustering, and
#     annotation BEFORE this UMI matrix + cell metadata were exported
#     as GEO supplementary files
#   - Confirmed independently: genes.tsv contains vst.* columns,
#     the output of Seurat's FindVariableFeatures() - i.e. this gene
#     set has already been through variable feature selection
#   - Confirmed independently: 0 of 13 canonical mitochondrial gene
#     Ensembl IDs are present in this matrix - MT genes were removed
#     upstream, so MT%-based filtering is not possible on this data
#
# CONCLUSION: this is NOT raw/unfiltered data. Doublet removal and
# low-quality-cell filtering were already performed by the original
# authors. This script therefore VERIFIES data quality rather than
# performing new filtering - re-filtering already-curated data with
# arbitrary cutoffs would risk discarding valid cells for no reason.
#
# Input:  data/processed/seurat_raw.rds
# Output: results/tables/qc_summary.csv (no cells removed)
# ============================================================

library(Seurat)

seurat_obj <- readRDS("data/processed/seurat_raw.rds")

# ---- Verify count distributions ----
count_summary   <- summary(seurat_obj$nCount_RNA)
feature_summary <- summary(seurat_obj$nFeature_RNA)

print(count_summary)
print(feature_summary)

# A solid non-near-zero floor on both nCount_RNA (min 1345) and
# nFeature_RNA (min 954) is consistent with already-filtered data.
# Raw/unfiltered droplets typically show large numbers of near-zero
# values (empty droplets); we don't see that pattern here.
stopifnot(min(seurat_obj$nCount_RNA) > 500)
stopifnot(min(seurat_obj$nFeature_RNA) > 200)

# ---- Per-sample summary, for the record ----
qc_by_sample <- data.frame(
  patient       = names(table(seurat_obj$patient)),
  n_cells       = as.integer(table(seurat_obj$patient)),
  median_nCount = tapply(seurat_obj$nCount_RNA, seurat_obj$patient, median),
  median_nFeature = tapply(seurat_obj$nFeature_RNA, seurat_obj$patient, median)
)
print(qc_by_sample)

# ---- Save QC record (documentation, not a filtering output) ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(qc_by_sample, "results/tables/qc_summary.csv", row.names = FALSE)

# No cells are removed in this script - seurat_raw.rds remains the
# object to carry forward into normalization (Script 3).
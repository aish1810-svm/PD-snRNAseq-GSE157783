# ============================================================
# 07_de_validation.R
# Purpose: Cross-validate the 24 microglia-specific DE genes from
#          06_microglia_de.R using an independent statistical
#          approach - CPM normalization + Wilcoxon rank-sum test -
#          rather than relying on DESeq2 alone.
#
# RATIONALE: DESeq2 uses model-based normalization (size factors) and
# a negative-binomial framework. Wilcoxon on simple CPM-normalized
# values makes different assumptions entirely. Agreement between the
# two methods is stronger evidence of a real signal than either
# method alone; disagreement would flag a result as method-dependent
# and worth treating with more caution.
#
# NOTE: with only 11 pseudobulk samples (5 PD, 6 Control), Wilcoxon
# can only produce a small number of distinct p-values (a known
# property of rank tests with small n) - identical p-values across
# multiple genes are expected and not an error.
#
# Input:  data/processed/microglia_pseudobulk_mat.rds
#         data/processed/microglia_specific_genes.rds
#         (both produced by 06_microglia_de.R)
# Output: results/tables/microglia_de_validation.csv
# ============================================================

pseudobulk_mat     <- readRDS("data/processed/microglia_pseudobulk_mat.rds")
microglia_specific <- readRDS("data/processed/microglia_specific_genes.rds")

# ---- 1. CPM-normalize pseudobulk counts (independent of DESeq2) ----
cpm_mat <- sweep(pseudobulk_mat, 2, colSums(pseudobulk_mat), FUN = "/") * 1e6
log_cpm_mat <- log2(cpm_mat + 1)

# ---- 2. Wilcoxon test on the 24 microglia-specific genes ----
patient_condition <- ifelse(grepl("^PD", colnames(log_cpm_mat)), "PD", "Control")

wilcox_results <- lapply(microglia_specific$gene_id, function(g) {
  vals <- log_cpm_mat[g, ]
  test <- tryCatch(
    wilcox.test(vals ~ patient_condition),
    error = function(e) NULL
  )
  data.frame(
    gene_id = g,
    wilcox_p = if (!is.null(test)) test$p.value else NA
  )
})
wilcox_df <- do.call(rbind, wilcox_results)

# ---- 3. Compare DESeq2 vs Wilcoxon results ----
comparison <- merge(
  microglia_specific[, c("gene_id", "gene_symbol", "log2FoldChange", "padj")],
  wilcox_df,
  by = "gene_id"
)
comparison <- comparison[order(comparison$padj), ]

print(comparison)

# ---- Save ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, "results/tables/microglia_de_validation.csv", row.names = FALSE)
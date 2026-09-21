# ============================================================
# 06_microglia_de.R
# Purpose: Pseudobulk differential expression analysis comparing PD vs
#          Control specifically within Microglia, using DESeq2.
#
# METHOD NOTE: cells are aggregated (summed) per patient into
# pseudobulk profiles before testing, rather than treating individual
# cells as independent replicates. Testing at the cell level would
# massively inflate statistical power by ignoring that cells from the
# same patient are not independent observations (pseudoreplication) -
# a well-known pitfall in single-cell DE analysis. Pseudobulk treats
# each of the 11 patients as one true biological replicate, matching
# standard bulk RNA-seq practice.
#
# IMPORTANT CAVEAT - AMBIENT RNA CONTAMINATION:
# Of 621 genes significant at padj < 0.1, only 24 (3.9%) have their
# highest expression in Microglia when checked against all cell types.
# The remaining 597 genes - including most top hits by p-value - peak
# in other cell types, overwhelmingly Oligodendrocytes (the most
# abundant cell type in this tissue). This is consistent with ambient
# RNA contamination: background RNA from abundant cell types leaking
# into other cells' droplets/nuclei during dissociation, a known
# technical artifact in snRNA-seq. Because oligodendrocyte proportion
# itself trends differently between PD and Control (see
# 05_celltype_proportions.R), this contamination signature can
# masquerade as condition-associated "differential expression" in
# microglia without reflecting real microglial biology.
#
# We therefore report two result sets:
#   1. Full DE results (all 621 significant genes) - for transparency
#   2. Filtered "microglia-specific" results (24 genes, where Microglia
#      is each gene's top-expressing cell type) - the higher-confidence
#      candidate list
# Even the filtered list should be interpreted cautiously given the
# small sample size (11 pseudobulk samples total).
#
# Input:  data/processed/seurat_clustered.rds
# Output: results/tables/microglia_de_full.csv
#         results/tables/microglia_de_specific_filtered.csv
# ============================================================

library(Seurat)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

seurat_obj <- readRDS("data/processed/seurat_clustered.rds")

# ---- 1. Subset to Microglia ----
microglia <- subset(seurat_obj, subset = cell_ontology == "Microglia")

# ---- 2. Build pseudobulk matrix (sum counts per patient) ----
patients <- unique(microglia$patient)

pseudobulk_list <- lapply(patients, function(p) {
  cells_p <- colnames(microglia)[microglia$patient == p]
  rowSums(microglia@assays$RNA$counts[, cells_p, drop = FALSE])
})

pseudobulk_mat <- do.call(cbind, pseudobulk_list)
colnames(pseudobulk_mat) <- patients

# ---- 3. Run DESeq2 ----
coldata <- data.frame(
  patient   = colnames(pseudobulk_mat),
  condition = ifelse(grepl("^PD", colnames(pseudobulk_mat)), "PD", "Control")
)
rownames(coldata) <- coldata$patient

stopifnot(all(colnames(pseudobulk_mat) == rownames(coldata)))
coldata$condition <- factor(coldata$condition, levels = c("Control", "PD"))

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_mat,
  colData   = coldata,
  design    = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "PD", "Control"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# ---- 4. Map gene symbols ----
res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# ---- 5. Flag likely ambient RNA contamination ----
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.1, ]

avg_exp_all <- AggregateExpression(
  seurat_obj,
  features = sig_genes$gene_id,
  group.by = "cell_ontology"
)$RNA
avg_exp_all <- as.data.frame(avg_exp_all)
avg_exp_all$gene_id <- rownames(avg_exp_all)
avg_exp_all$top_celltype <- colnames(avg_exp_all)[apply(avg_exp_all[, 1:12], 1, which.max)]

sig_genes <- merge(sig_genes, avg_exp_all[, c("gene_id", "top_celltype")], by = "gene_id")
sig_genes$likely_microglia_specific <- sig_genes$top_celltype == "Microglia"
sig_genes <- sig_genes[order(sig_genes$padj), ]

microglia_specific <- sig_genes[sig_genes$likely_microglia_specific, ]

print(table(sig_genes$likely_microglia_specific))
print(microglia_specific[, c("gene_id", "gene_symbol", "log2FoldChange", "padj", "top_celltype")])

# ---- 6. Save ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(sig_genes, "results/tables/microglia_de_full.csv", row.names = FALSE)
write.csv(microglia_specific, "results/tables/microglia_de_specific_filtered.csv", row.names = FALSE)

# Save intermediate objects so 07_de_validation.R can run standalone,
# independent of this script's in-memory state
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(pseudobulk_mat, "data/processed/microglia_pseudobulk_mat.rds")
saveRDS(microglia_specific, "data/processed/microglia_specific_genes.rds")
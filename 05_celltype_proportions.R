# ============================================================
# 05_celltype_proportions.R
# Purpose: Compare cell-type proportions between PD and Control
#          samples, using author-provided cell_ontology labels.
#
# Input:  data/processed/seurat_clustered.rds
# Output: results/tables/celltype_proportions.csv
#         results/tables/proportion_test_results.csv
#         results/figures/celltype_proportions_boxplot.png
#
# NOTE ON INTERPRETATION: with only 11 samples split across 12 cell
# types, this is a low-powered comparison. After Benjamini-Hochberg
# correction for the 12 simultaneous tests, no cell type reaches the
# conventional 0.05 significance threshold. Microglia shows the
# strongest trend (uncorrected p = 0.0087, adjusted p = 0.104) and is
# directionally consistent with known PD neuroinflammation biology,
# but should be reported as a suggestive trend, not a confirmed
# finding. DaNs show p = 1 due to near-zero proportions in most
# samples of both conditions - a power/sparsity limitation, not
# evidence of "no difference."
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

seurat_obj <- readRDS("data/processed/seurat_clustered.rds")

# ---- 1. Build per-sample cell-type proportion table ----
cell_counts <- as.data.frame(table(seurat_obj$patient, seurat_obj$cell_ontology))
colnames(cell_counts) <- c("patient", "cell_type", "n_cells")

sample_totals <- as.data.frame(table(seurat_obj$patient))
colnames(sample_totals) <- c("patient", "total_cells")

cell_props <- merge(cell_counts, sample_totals, by = "patient")
cell_props$proportion <- cell_props$n_cells / cell_props$total_cells

condition_map <- unique(data.frame(patient = seurat_obj$patient, condition = seurat_obj$condition))
cell_props <- merge(cell_props, condition_map, by = "patient")

stopifnot(nrow(cell_props) == length(unique(cell_props$patient)) * length(unique(cell_props$cell_type)))

# ---- 2. Visualize ----
ggplot(cell_props, aes(x = cell_type, y = proportion, fill = condition)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell type proportions: PD vs Control", y = "Proportion of cells", x = "Cell type")

ggsave("results/figures/celltype_proportions_boxplot.png", width = 12, height = 6, dpi = 300)

# ---- 3. Statistical testing (Wilcoxon, per cell type) ----
prop_tests <- cell_props %>%
  group_by(cell_type) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(proportion ~ condition)$p.value,
      error = function(e) NA
    )
  ) %>%
  arrange(p_value)

# Multiple-testing correction across the 12 cell types tested
prop_tests$p_adj <- p.adjust(prop_tests$p_value, method = "BH")

print(prop_tests)

# ---- 4. Save ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(cell_props, "results/tables/celltype_proportions.csv", row.names = FALSE)
write.csv(prop_tests, "results/tables/proportion_test_results.csv", row.names = FALSE)
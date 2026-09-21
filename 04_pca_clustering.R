# ============================================================
# 04_pca_clustering.R
# Purpose: Run PCA, cluster cells, generate UMAP embedding, validate
#          our independent clustering against the cell_ontology labels
#          already provided by the original authors, and save a
#          comparison figure.
#
# Input:  data/processed/seurat_normalized.rds
# Output: data/processed/seurat_clustered.rds
#         results/tables/cluster_vs_celltype.csv
#         results/figures/umap_clusters_vs_celltype.png
# ============================================================

library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("data/processed/seurat_normalized.rds")

# ---- 1. PCA ----
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# ElbowPlot(seurat_obj, ndims = 50) was inspected interactively:
# curve flattens around PC10. We use the first 15 PCs (slightly past
# the visual elbow) for downstream steps, since missing real signal
# is worse than including a few marginal PCs.

# ---- 2. Clustering ----
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# ---- 3. UMAP (visualization only - does not affect cluster assignment) ----
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# ---- 4. Validate clusters against author-provided cell type labels ----
# NOTE: both panels below use OUR OWN UMAP embedding (computed entirely
# from our own PCA/clustering pipeline in this script). No coordinates,
# figures, or analysis from the original paper are used anywhere. The
# right-hand panel colors our UMAP using the cell_ontology column that
# was included as author-provided metadata in the GEO supplementary
# files (IPDCO_hg_midbrain_cell.tsv) - i.e. only the cell type NAMES
# came from the authors, not the embedding, clustering, or any figure.
# This is a standard validation check: does our independent clustering
# land in the same neighborhoods as the authors' expert-assigned labels?
cluster_vs_celltype <- table(seurat_obj$seurat_clusters, seurat_obj$cell_ontology)
print(cluster_vs_celltype)

# ---- 5. Comparison figure ----
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Our UMAP, colored by our own clusters")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_ontology", label = TRUE) +
  ggtitle("Our UMAP, colored by author-provided cell type labels")

combined_plot <- p1 + p2

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/umap_clusters_vs_celltype.png", combined_plot,
       width = 16, height = 6, dpi = 300)

# ---- 6. Save ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame.matrix(cluster_vs_celltype),
          "results/tables/cluster_vs_celltype.csv")

saveRDS(seurat_obj, "data/processed/seurat_clustered.rds")
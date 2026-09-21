# ============================================================
# 01_load_data.R
# Purpose: Load GSE157783 snRNA-seq data (UMI matrix, cell metadata,
#          gene metadata) and assemble into a Seurat object.
# Input:   data/raw/IPDCO_hg_midbrain_UMI.tsv
#          data/raw/IPDCO_hg_midbrain_cell.tsv
#          data/raw/IPDCO_hg_midbrain_genes.tsv
# Output:  data/processed/seurat_raw.rds
# ============================================================

library(data.table)
library(Matrix)
library(Seurat)

# ---- 1. Load metadata files ----
cell_df  <- read.delim("data/raw/IPDCO_hg_midbrain_cell.tsv",  stringsAsFactors = FALSE)
genes_df <- read.delim("data/raw/IPDCO_hg_midbrain_genes.tsv", stringsAsFactors = FALSE)

# genes.tsv's "row" column must be sequential 1..N for row-position matching to be valid
stopifnot(all(genes_df$row == seq_len(nrow(genes_df))))

# ---- 2. Load UMI count matrix ----
umi_dt <- fread("data/raw/IPDCO_hg_midbrain_UMI.tsv", header = TRUE)

# Structural checks confirmed interactively - re-verified here so the script
# fails loudly if run on different/corrupted files rather than silently continuing
stopifnot(nrow(umi_dt) == nrow(genes_df))
stopifnot(ncol(umi_dt) == nrow(cell_df))
stopifnot(all(colnames(umi_dt) == cell_df$barcode))

# ---- 3. Convert to sparse matrix in chunks ----
# Converting the full ~2GB table to a dense matrix in one shot requires holding
# both the data.table AND an ~8GB dense copy in memory simultaneously, which
# exceeds available RAM on a typical machine. Instead, convert small column
# (cell) chunks to sparse format one at a time, so memory use stays low.
n_cells    <- ncol(umi_dt)
chunk_size <- 5000
n_chunks   <- ceiling(n_cells / chunk_size)
sparse_list <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  start_col <- (i - 1) * chunk_size + 1
  end_col   <- min(i * chunk_size, n_cells)
  
  chunk_mat <- as.matrix(umi_dt[, start_col:end_col, with = FALSE])
  sparse_list[[i]] <- as(chunk_mat, "CsparseMatrix")
  rm(chunk_mat)
  gc()
  
  cat("Processed columns", start_col, "to", end_col, "of", n_cells, "\n")
}

umi_mat <- do.call(cbind, sparse_list)
rm(sparse_list, umi_dt)
gc()

rownames(umi_mat) <- genes_df$gene
colnames(umi_mat) <- cell_df$barcode

# Re-verify final dimensions after chunked reassembly
stopifnot(dim(umi_mat)[1] == nrow(genes_df))
stopifnot(dim(umi_mat)[2] == nrow(cell_df))

# ---- 4. Build per-cell metadata ----
meta <- data.frame(
  patient       = cell_df$patient,
  cell_ontology = cell_df$cell_ontology,
  row.names     = cell_df$barcode
)

# Check actual patient labels before assuming any naming pattern
print(table(meta$patient))

meta$condition <- ifelse(grepl("^PD", meta$patient), "PD", "Control")

# ---- 5. Create the Seurat object ----
seurat_obj <- CreateSeuratObject(
  counts    = umi_mat,
  project   = "GSE157783_PD_snRNAseq",
  meta.data = meta
)

# ---- 6. Sanity checks ----
seurat_obj
table(seurat_obj$condition)
table(seurat_obj$patient)

# ---- 7. Save ----
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(seurat_obj, "data/processed/seurat_raw.rds")
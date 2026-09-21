# Single-Nucleus RNA-seq Analysis of Parkinson's Disease Substantia Nigra (GSE157783)

## Overview

This project analyzes single-nucleus RNA-seq (snRNA-seq) data from human
midbrain tissue (6 control, 5 idiopathic Parkinson's disease samples),
originally published by Smajić et al., *Brain* 2022 ("Single-cell sequencing
of human midbrain reveals glial activation and a Parkinson-specific neuronal
state"). Raw data source: [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783).

This is an independent re-analysis built from the GEO supplementary files
(pooled UMI count matrix + cell/gene metadata), not a reproduction of the
original authors' code. The goal was to build a full snRNA-seq pipeline from
raw files through to a defensible biological finding, using Seurat and
DESeq2, and to validate results rigorously rather than take them at face
value.

## Data

- **Source**: GEO accession GSE157783, 3 supplementary files (UMI matrix,
  cell metadata with author-assigned cell type labels, gene metadata)
- **Scale**: 26,737 genes × 41,435 cells (nuclei), across 11 samples
  (5 PD, 6 Control)
- **Note**: this data was already processed upstream by the original authors
  (doublet removal via Scrublet, clustering/annotation via Seurat3 +
  monocle3, confirmed via the authors' GitHub repository) before being
  deposited as GEO supplementary files. Only raw counts were included in
  the supplementary files, so normalization/clustering were re-run
  independently here.

## Pipeline

| Script | Purpose |
|---|---|
| `01_load_data.R` | Load raw GEO files, build Seurat object |
| `02_qc_verification.R` | Verify (not re-filter) data quality; document upstream QC |
| `03_normalize.R` | Log-normalization, variable feature selection, scaling |
| `04_pca_clustering.R` | PCA, Louvain clustering, UMAP, validation against author cell-type labels |
| `05_celltype_proportions.R` | PD vs Control cell-type composition, Wilcoxon tests, BH correction |
| `06_microglia_de.R` | Pseudobulk DESeq2: PD vs Control differential expression in Microglia, ambient RNA contamination check |
| `07_de_validation.R` | Cross-validation of DE results using an independent method (Wilcoxon on CPM) |

Run in order. Requires R 4.6.1, Rtools45, and the packages: Seurat,
data.table, Matrix, DESeq2, org.Hs.eg.db, AnnotationDbi, dplyr, ggplot2,
patchwork.

## Key Findings

**1. Clustering validation.** Independent clustering (18 clusters, modularity
0.94) closely recovered the cell-type structure the original authors
identified with their own pipeline — most clusters mapped cleanly to a
single dominant cell type (e.g., Oligodendrocytes, Microglia, Astrocytes),
confirmed both numerically (cross-tabulation) and visually (UMAP).

**2. Cell-type proportions.** Microglia showed the strongest PD vs Control
proportion trend (uncorrected p = 0.0087), directionally consistent with
known PD neuroinflammation biology. After Benjamini-Hochberg correction for
multiple testing (12 cell types), this did not reach the conventional 0.05
threshold (adjusted p = 0.104) — reported as a suggestive trend, not a
confirmed finding.

**3. Microglia differential expression (PD vs Control).** Pseudobulk DESeq2
(patients as biological replicates, avoiding cell-level pseudoreplication)
identified 621 genes significant at padj < 0.1.

**4. Ambient RNA contamination.** Cross-checking each significant gene
against its expression across all cell types revealed that 597 of 621
genes (96%) actually peak in other cell types — predominantly
Oligodendrocytes, the most abundant cell type in this tissue. This is
consistent with ambient RNA contamination, a known technical artifact in
snRNA-seq. Filtering to the 24 genes where Microglia is the top-expressing
cell type produced a higher-confidence candidate list.

**5. GPNMB.** Among the 24 microglia-specific candidates, GPNMB stood out:
upregulated in PD (log2FC = +1.40, DESeq2 padj = 0.010), and independently
confirmed via a second, unrelated statistical method (CPM + Wilcoxon,
p = 0.0043). GPNMB is a well-established Parkinson's disease genetic risk
gene, previously linked to activated/disease-associated microglia in
neurodegeneration. Recovering this signal independently, in a different
dataset than where it was originally discovered, is a genuine (if modest,
given n=11) validation result.

**6. External validation against the original publication.** The original
Smajić et al. (2022, *Brain*) paper independently reports both an increase
in nigral microglia proportion in PD and GPNMB specifically as a defining
marker of a pro-inflammatory, disease-associated microglia subpopulation
in their own (differently-implemented) analysis pipeline. This external,
independently-published confirmation - using different code on the same
underlying tissue - substantially strengthens confidence that the
microglia proportion trend and GPNMB findings reflect real biology rather
than a pipeline artifact.

## Limitations

- **Small sample size**: 11 total samples (5 PD, 6 Control) throughout.
  All statistical tests, especially the proportion and DE analyses, are
  underpowered by conventional standards; results should be read as
  hypothesis-generating, not confirmatory.
- **Sparse dopaminergic neurons (DaNs)**: despite being the cell type most
  directly relevant to PD pathology, DaNs were too rare (74 cells total)
  to form a distinct cluster or support meaningful statistical testing
  (proportion test: p = 1, driven by near-zero values in most samples of
  both conditions). This reflects known DaN depletion in PD tissue and
  is a real limitation of this dataset/analysis, not an error.
- **Ambient RNA contamination**: affects the whole dataset, not just the
  genes explicitly filtered in Script 6. Any single-cell/nucleus finding
  from this data should be interpreted with this in mind, especially for
  comparisons involving Oligodendrocytes (the dominant, most contamination-
  prone cell type here).
- **Reliance on pre-processed input data**: raw FASTQ files were not used;
  this analysis builds on the authors' already-curated supplementary
  matrix, inheriting any upstream processing decisions (e.g. gene set
  restriction — only 26,737 of the full transcriptome, no mitochondrial
  genes present) that could not be independently verified or altered.
- **Condition labels (PD vs Control)**: derived from GEO sample naming
  conventions (IPD1-IPD5, C1-C6) matched against the Series-level
  description ("...Idiopathic Parkinson's Disease patients and matched
  controls"), not confirmed against each individual sample's own GEO
  metadata page. This is well-grounded but not independently verified
  at the per-sample level.

## Reproducibility

All processed intermediate objects (`data/processed/*.rds`) and result
tables (`results/tables/*.csv`) are saved at each pipeline stage. Scripts
were run interactively during development; a fresh-session, start-to-finish
rerun of all 7 scripts in order is recommended before treating this as a
final, polished record.

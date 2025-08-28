# CITEscope

An integrated and user, friendly R package + Shiny app for comprehensive **CITEseq** analysis.

## Prerequisites
**Requires R >= 4.1**
(The upper limit for large file uploads has been set inside the package. No additional operations are required.)

---

## Installation

```r
# If devtools is not installed
install.packages("devtools")

# Install CITEscope from GitHub
devtools::install_github("Weichen-S/CITEscope")

```

## Bioconductor dependencies(if missing)
```r
install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex", "org.Mm.eg.db"))
```
## For Seurat → Monocle3 conversion (optional)
```r
BiocManager::install("SeuratWrappers")
```
## Quick start
```r
library(CITEscope)
CITEscope::run_app()
```
# About
---
CITEscope is an integrated R package and visualization platform designed for comprehensive CITE-seq data analysis. It streamlines preprocessing, multimodal integration, annotation, Infection Status Classification,functional module scoring, trajectory inference, and downstream interpretation within a unified workflow. With its interactive visualization interface, CITEscope enables users to intuitively explore complex single-cell multi-omics datasets.
<img width="838" height="472" alt="new" src="https://github.com/user-attachments/assets/0dc7d4f4-a528-4931-8598-750c477fa1ac" />

# Workflow
---
## Data Preprocessing & Quality Control
**Input:** 

10x Genomics outputs (ZIP of matrix.mtx, features.tsv, barcodes.tsv) or uploaded TSV.

ADT Requirements: The recommended data should include ADT (antibody tag) for the use of multimodal analysis.

**Processing steps:**

1.Split RNA and ADT modalities (based on feature type).

2.RNA normalization: SCTransform (method = "glmGamPoi", min.cells = 3, min.features = 200).

3.ADT normalization: CLR + ScaleData.

4.Identify variable features for RNA and ADT.

**Output:** 

A Seurat object containing RNA + ADT assays; metadata preview table (first 10 rows).


## WNN Clustering
**Functions:** 

Run PCA on RNA (SCT) and ADT, build WNN graph, UMAP reduction, Louvain clustering.

**Output:**

WNN UMAP + clusters

SCT.weight and ADT.weight FeaturePlots


## Automatic Cell Type Annotation
Automated annotation via SingleR (Aran et al., 2019).

**Reference datasets** (from celldex):

1. Human: HumanPrimaryCellAtlasData

2. Mouse: MouseRNAseqData


Marker support: Predefined immune marker panel (e.g., CD3, CD4, CD8, CD19, CD14, NCAM1, EPCAM).

**provided：** 

Core marker presence/absence table + full searchable gene list.

**Output:**
SingleR-labeled UMAP.

ADT marker UMAP 

SCT top variable feature UMAP.

Custom marker UMAP

violin plots.


## Infection Status Classification
Based on NP-ADT signal.

Thresholding NP expression to assign infected vs bystander cells.

**Output :**

Violin plot: NP-ADT distribution.

Bar plot: Proportions of infected vs bystander.

UMAP colored by infection status
## Functional Module Scoring
**Built-in gene sets:**

1. Hallmark IFN-α + IFN-γ (MSigDB H).

2. MSigDB C7 Interferon sets.

3. Cross-species Core ISGs (human, mouse, cow, sheep).

4. PRR signaling (DDX58, IFIH1, MAVS, TBK1, IRF3, IRF7).

**Additional modules:**

1. Exhaustion (MSigDB C7 “exhaust” sets).

2. Cytokine/Inflammatory Response (MSigDB H “inflammatory response”).

**Custom:** User-uploaded CSV (first column = gene list).

**Output:**

UMAP plots for ISG, Exhaustion, Cytokine modules.

Combined UMAP (three modules merged).

DotPlot (Group × Module scores).

Heatmap (group mean scores).

RidgePlot (distribution across groups).

Pseudobulk Bar plot (mean ± SEM).

Scatter correlation: ISG_Score vs NP-ADT 

## Trajectory Inference & Graph Learning
Root selection:

Manual cluster, or

ADT-based (select cluster with highest quantile of NP-ADT).

ADT integration:

ADT markers projected on trajectory UMAP.

ADT vs pseudotime smooth curve fitting.

**Output:**

UMAP + learned trajectory graph (pseudotime).

UMAP colored by ADT marker.

ADT vs pseudotime smoothed curve.


## Differential Expression & Pathway Enrichment
**Output:**

MA plot

Volcano plot

Heatmap of top 10 DEGs

GO/KEGG enrichment plots (dot/bar)

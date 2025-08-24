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


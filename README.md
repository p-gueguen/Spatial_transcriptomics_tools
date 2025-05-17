# 🌌 Spatial Analyses

**Other resource** → [awesome_spatial_omics](https://github.com/crazyhottommy/awesome_spatial_omics)

- [Spatial-omic SIB tutorial](https://elixir-europe-training.github.io/ELIXIR-SCO-spatial-omics/presentations.html)

## Table of Contents
- [General Tools](#general-tools)
- [Analysis Pipeline Steps](#analysis-pipeline-steps)
  - [QC](#qc)
  - [Normalization](#normalization)
  - [Denoising](#denoising)
  - [Bias Correction](#bias-correction)
  - [Cell Segmentation](#cell-segmentation)
  - [Cell Annotation](#cell-annotation)
  - [Cell Deconvolution](#cell-deconvolution)
  - [Differential Expression](#differential-expression)
  - [Spatially Variable Genes](#spatially-variable-genes)
  - [Integration](#integration)
  - [Cell Niches & Tissue Domains](#cell-niches--tissue-domains)
  - [Cell Distances & Neighborhood](#cell-distances--neighborhood)
  - [Spatial Trajectories](#spatial-trajectories)
  - [Cell-Cell Communication](#cell-cell-communication)
  - [Metacells & Scalability](#metacells--scalability)
  - [Subcellular Analysis](#subcellular-analysis)
  - [Copy Number Variations](#copy-number-variations)
  - [Transcription Factors & Gene Regulatory Networks](#transcription-factors--gene-regulatory-networks)
- [Technical Enhancements](#technical-enhancements)
  - [Slide Alignment](#slide-alignment)
  - [Super Resolution](#super-resolution)
  - [Transcripts + Histology](#transcripts--histology)
- [Benchmarks](#benchmarks)
- [Datasets & Foundation Models](#datasets--foundation-models)

## General Tools

- [Best practices Bioconductor](https://lmweber.org/OSTA/) - Principles for statistical analysis of spatial transcriptomics data
- [squidpy](https://squidpy.readthedocs.io/en/stable/) - Spatial single cell analysis toolkit from scverse
- [Giotto](https://giottosuite.readthedocs.io/en/latest/) - Comprehensive spatial data analysis suite
- [Vitessce](http://vitessce.io/) - Visual integration tool for exploration of spatial single cell experiments
- [Voyager](https://github.com/pachterlab/voyager) - Spatial transcriptomics visualization from Pachter lab
- [BASS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02734-7) - Multiple sample analysis
- [SpaVAE](https://github.com/ttgump/spaVAE) - All-purpose tool for dimension reduction, visualization, clustering, batch integration, denoising, differential expression, spatial interpolation, and resolution enhancement
- [sopa](https://github.com/gustaveroussy/sopa) - Spatial omics processing and analysis
- [SpatialAgent](https://www.biorxiv.org/content/10.1101/2025.04.03.646459v1.full.pdf) - An autonomous AI agent for spatial biology
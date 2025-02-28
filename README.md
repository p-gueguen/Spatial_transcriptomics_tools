# 🌌 Spatial Analyses

**Other resource** → [awesome_spatial_omics](https://github.com/crazyhottommy/awesome_spatial_omics)

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
  - [Spatially Variable Genes](#spatially-variable-genes)
  - [Cell Niches & Tissue Domains](#cell-niches--tissue-domains)
  - [Cell Distances & Neighborhood](#cell-distances--neighborhood)
  - [Spatial Trajectories](#spatial-trajectories)
  - [Cell-Cell Communication](#cell-cell-communication)
  - [Metacells](#metacells)
  - [Subcellular Analysis](#subcellular-analysis)
  - [Copy Number Variations](#copy-number-variations)
  - [Transcription Factors & Gene Regulatory Networks](#transcription-factors--gene-regulatory-networks)
- [Technical Enhancements](#technical-enhancements)
  - [Slide Alignment](#slide-alignment)
  - [Super Resolution](#super-resolution)
- [Benchmarks](#benchmarks)
- [Datasets & Foundation Models](#datasets--foundation-models)

## General Tools

- [Best practices Bioconductor](https://lmweber.org/PrinciplesSTA/devel/) - Principles for statistical analysis of spatial transcriptomics data
- [squidpy](https://squidpy.readthedocs.io/en/stable/) - Spatial single cell analysis toolkit from scverse
- [Giotto](https://giottosuite.readthedocs.io/en/latest/) - Comprehensive spatial data analysis suite
- [Vitessce](http://vitessce.io/) - Visual integration tool for exploration of spatial single cell experiments
- [Voyager](https://github.com/pachterlab/voyager) - Spatial transcriptomics visualization from Pachter lab
- [BASS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03361-0) - Multiple sample analysis
- [SpaVAE](https://github.com/ttgump/spaVAE) - All-purpose tool for dimension reduction, visualization, clustering, batch integration, denoising, differential expression, spatial interpolation, and resolution enhancement
- [sopa](https://github.com/gustaveroussy/sopa?tab=readme-ov-file) - Spatial omics processing and analysis

## Analysis Pipeline Steps

### QC

- [SpaceTrooper](https://htmlpreview.github.io/?https://github.com/drighelli/SpaceTrooper/blob/main/vignette/introduction.html) - Quality control for spatial transcriptomics

### Normalization

- [Cell volume normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w#Fig7) - Recommended for imaging-based techniques, especially with small probe lists

### Denoising

- **Note:** Gene imputation is [not recommended](https://github.com/BayraktarLab/cell2location/issues/379)
- [Sprod](https://www.nature.com/articles/s41592-022-01560-w#Fig2) - Spatial denoising method

### Bias Correction

- [ResolVI](https://www.biorxiv.org/content/biorxiv/early/2025/01/24/2025.01.20.634005.full.pdf) - Bias correction method
- [Statial](https://sydneybiox.github.io/Statial/) - Correction of spill-over effects

### Cell Segmentation

- [Baysor](https://github.com/kharchenkolab/Baysor) - Bayesian segmentation of spatial transcriptomics data
- [Cellpose](https://github.com/MouseLand/cellpose) - Generalist algorithm for cellular segmentation
  - [Cellpose 3](https://www.biorxiv.org/content/10.1101/2024.02.10.579780v2) - With supersampling/restoration capabilities
- [DeepCell](https://github.com/vanvalenlab/deepcell-tf) - Deep learning library for single cell analysis
- [Bo Wang's method](https://www.nature.com/articles/s41592-024-02233-6#Fig3) - Better than SOTA segmentation (Nature Methods 2024)
- [Bin2Cell](https://www.biorxiv.org/content/10.1101/2024.06.19.599766v1) - Segmentation for VisiumHD data
- [ENACT](https://www.biorxiv.org/content/10.1101/2024.10.17.618905v1.full.pdf) - Segmentation for VisiumHD data
- [Proseg](https://www.biorxiv.org/content/10.1101/2024.04.25.591218v1.full.pdf) - Probabilistic segmentation method
- [ComSeg](https://www.nature.com/articles/s42003-024-06480-3) - Transcript-based segmentation
- [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) - Feature-based image segmentation
- [Xenium cell boundary](https://kb.10xgenomics.com/hc/en-us/articles/30205122555917-Boundary-stain-shows-cell-morphology-better-than-the-interior-stain) - Alternative when interior staining fails
- [Bioimage.io](https://bioimage.io/#/) - Repository of AI models for segmentation
- [ST-cellseg](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012254) - Segmentation for spatial transcriptomics
- [CelloType](https://github.com/maxpmx/CelloType) - Cell type detection and segmentation

### Spatially Variable Genes

- [PROST](https://www.nature.com/articles/s41467-024-44835-w) - Detection of spatially variable genes
- [SPARK-X](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0#Fig1) - Detection of spatially variable genes, best performing

### Cell Annotation

- [STEM](https://www.nature.com/articles/s42003-023-05640-1) - Cell type annotation method
- [TACIT](https://www.biorxiv.org/content/10.1101/2024.05.31.596861v1) - Automated cell type identification
- [moscot](https://github.com/theislab/moscot) - Optimal transport-based cell mapping
- [RCTD](https://github.com/dmcable/spacexr) - Cell type annotation from reference data
- [Annotability](https://www.nature.com/articles/s43588-024-00721-5#Fig2) - Assessment of cell annotation quality
- [CELLama](https://www.biorxiv.org/content/10.1101/2024.05.08.593094v1.full) - Cell annotation model
- [TACCO](https://www.nature.com/articles/s41587-023-01657-3) - Transfer of annotations between single-cell datasets
- [TANGRAM](https://github.com/broadinstitute/Tangram) - Mapping single-cell to spatial data

### Cell Deconvolution

- [RCTD](https://github.com/dmcable/spacexr) - Robust cell type decomposition
- [Cell2location](https://github.com/BayraktarLab/cell2location) - Mapping scRNA-seq to spatial data

### Cell Niches & Tissue Domains

> **(Smaller)** Cell types → Cell modules/neighborhoods → Niches/tissue domains **(Larger)**

- [BANKSY](https://www.nature.com/articles/s41588-024-01664-3) - Unified cell typing and tissue domain segmentation
- [staryfish](https://twitter.com/Cancer_dynamics/status/1770811468578971905) - Tissue domain detection
- [TISSUE](https://buff.ly/49Oc0M2) - Tissue domain identification
- [CellCharter](https://www.nature.com/articles/s41588-023-01588-4) - Hierarchical niche detection
- [SpatialGLUE](https://www.nature.com/articles/s41592-024-02316-4) - Multi-omics cell niche identification
- [smoothclust](https://github.com/lmweber/smoothclust) - Spatial clustering
- [SpaTopic](https://www.science.org/doi/10.1126/sciadv.adp4942) - Spatial topic modeling
- [hdWGCNA](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) - Weighted gene correlation network analysis
- [GASTON](https://github.com/raphael-group/GASTON) - Graph-based spatial domain detection
- [SpatialMNN](https://github.com/Pixel-Dream/spatialMNN) - Identification of shared niches between slides

### Cell Distances & Neighborhood

- [CRAWDAD](https://github.com/JEFworks-Lab/CRAWDAD) - Cell relationship analysis with directional adjacency distributions
- [HoodscanR](https://www.biorxiv.org/content/10.1101/2024.03.26.586902v1) - Neighborhood analysis
- [SpicyR](https://sydneybiox.github.io/spicyR/) - Spatial analysis in R

### Spatial Trajectories

- [spaTrack](https://www.biorxiv.org/content/10.1101/2023.09.04.556175v2) - Spatial trajectory analysis

### Cell-Cell Communication

- [Spatia](https://www.nature.com/articles/s41592-024-02408-1#Fig2) - Spatial cell-cell interaction analysis
- [CellAgentChat](https://github.com/mcgilldinglab/CellAgentChat) - Agent-based cell communication modeling

### Metacells

- [MetaSpot](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btae734/7919601) - Metacell analysis for spatial data

### Subcellular Analysis

- [Sprawl](https://elifesciences.org/reviewed-preprints/87517) - Subcellular transcript localization
- [Bento](https://www.notion.so/Spatial-analyses-0d451532f4c64fc599cb6ceb469ab523?pvs=21) - Subcellular analysis
- [FISHfactor](https://academic.oup.com/bioinformatics/article/39/5/btad183/7114027) - Analysis of subcellular transcript patterns
- [InSTAnT](https://www.nature.com/articles/s41467-024-49457-w) - Intracellular spatial transcript analysis

### Copy Number Variations

- [CalicoST](https://www.nature.com/articles/s41592-024-02438-9#Fig1) - CNV detection in spatial data

### Transcription Factors & Gene Regulatory Networks

- [STAN](https://www.biorxiv.org/content/10.1101/2024.06.26.600782v2.full.pdf) - Spatial transcription factor analysis

## Technical Enhancements

### Slide Alignment

- [PASTE](https://www.nature.com/articles/s41592-022-01459-6)/[PASTE2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9881963/) - Probabilistic alignment of spatial transcriptomics experiments

### Super Resolution

- [TESLA](https://www.cell.com/cell-systems/pdf/S2405-4712(23)00084-4.pdf) - Super resolution for 10X Visium
- [istar](https://www.notion.so/Spatial-analyses-0d451532f4c64fc599cb6ceb469ab523?pvs=21) - Super resolution for Visium
- [BayesSPACE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763026/) - Subspot resolution

## Benchmarks

- [Deconvolution benchmark](https://www.biorxiv.org/content/10.1101/2024.04.03.586404v1.full) - Comprehensive comparison
- [RCTD and Cell2location benchmark](https://www.biorxiv.org/content/10.1101/2023.03.22.533802v3.full.pdf) - Claims these are the best methods
- [Spatial clustering benchmark](https://www.nature.com/articles/s41592-024-02215-8) - Comparison of clustering methods
- [Nature Communications review](https://www.nature.com/articles/s41467-023-37168-7#Fig2) - Confirms Cell2location performance
- [Open problems benchmark](https://openproblems.bio/results/spatial_decomposition/) - Cell2location is top performer

## Datasets & Foundation Models

### Foundation Models

- [scGPT-spatial](https://www.biorxiv.org/content/10.1101/2025.02.05.636714v1.full.pdf) - Language model for spatial transcriptomics
- [Phikon-v2](https://huggingface.co/owkin/phikon-v2) - Spatial biology foundation model
- [Bioptimus H-optimus-0](https://huggingface.co/bioptimus/H-optimus-0) - Biology-focused foundation model
- [DeepCell dataset](https://exploredata.deepcell.com/cell-visualizations/9/versions/1) - CNN + human features embeddings

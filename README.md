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

## Analysis Pipeline Steps

### QC

- [SpaceTrooper](https://htmlpreview.github.io/?https://github.com/drighelli/SpaceTrooper/blob/main/vignette/introduction.html) - Quality control for spatial transcriptomics
- [GrandQC](https://github.com/cpath-ukk/grandqc) - Comprehensive solution for quality control in digital pathology
- [SpotSweeper](https://www.nature.com/articles/s41592-025-02713-3#Fig1) - Spatially aware quality control for spatial transcriptomics
- [KRONOS](https://github.com/mahmoodlab/kronos) - Foundation Model for Multiplex Spatial Proteomic Images

### Normalization

- [Cell volume normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w#Fig7) - Recommended for imaging-based techniques, especially with small probe lists

### Denoising

- **Note:** Gene imputation is [not recommended](https://github.com/BayraktarLab/cell2location/issues/379)
- [Sprod](https://www.nature.com/articles/s41592-022-01560-w#Fig2) - Spatial denoising method

### Bias Correction

- [ResolVI](https://www.biorxiv.org/content/biorxiv/early/2025/01/24/2025.01.20.634005.full.pdf) - Bias correction method
- [Statial](https://sydneybiox.github.io/Statial/) - Correction of spill-over effects
- [ovrl.py](https://github.com/HiDiHlabs/ovrl.py) - A python tool to investigate vertical signal properties of imaging-based spatial transcriptomics data

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
- [SAINSC](https://onlinelibrary.wiley.com/doi/10.1002/smtd.202401123) - Segmentation for sequencing-based spatial data
- [BIDCell](https://github.com/SydneyBioX/BIDCell) - Biologically-informed deep learning for subcellular spatial transcriptomics segmentation
- [FastReseg](https://github.com/Nanostring-Biostats/FastReseg) - Using transcript locations to refine image-based cell segmentation results
- [Segger](https://www.biorxiv.org/content/10.1101/2025.03.14.643160v1) - Fast and accurate cell segmentation of imaging-based spatial transcriptomics data

### Cell Annotation

- [STEM](https://www.nature.com/articles/s42003-023-05640-1) - Cell type annotation method
- [TACIT](https://www.biorxiv.org/content/10.1101/2024.05.31.596861v1) - Automated cell type identification
- [moscot](https://github.com/theislab/moscot) - Optimal transport-based cell mapping
- [RCTD](https://github.com/dmcable/spacexr) - Cell type annotation from reference data
- [Annotability](https://www.nature.com/articles/s43588-024-00721-5#Fig2) - Assessment of cell annotation quality
- [CELLama](https://www.biorxiv.org/content/10.1101/2024.05.08.593094v1.full) - Cell annotation model
- [TACCO](https://www.nature.com/articles/s41587-023-01657-3) - Transfer of annotations between single-cell datasets
- [TANGRAM](https://github.com/broadinstitute/Tangram) - Mapping single-cell to spatial data
- [MMoCHi](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(24)00328-X) - Cell annotation method
- [CytoSPACE](https://www.nature.com/articles/s41587-023-01697-9) - High-resolution alignment of single-cell and spatial transcriptomes
- [ABCT](https://github.com/ercsb-sp/ABCT) - Cell type annotation method
- [STHD](https://www.biorxiv.org/content/10.1101/2024.06.20.599803v2) - Cell annotation for VisiumHD
- [STELLAR](https://www.nature.com/articles/s41592-022-01651-8) - Annotation of spatially resolved single-cell data with STELLAR
- [SPLIT](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full.pdf) - SPLIT effectively resolves mixed signals and enhances cell-type purity

### Cell Deconvolution

- [RCTD](https://github.com/dmcable/spacexr) - Robust cell type decomposition
- [Cell2location](https://github.com/BayraktarLab/cell2location) - Mapping scRNA-seq to spatial data
- [SPOTlight](https://github.com/MarcElosua/SPOTlight) - Seeded NMF regression to deconvolute spatial spots
- [CARD](https://github.com/YingMa0107/CARD) - Spatially informed cell-type deconvolution

### Differential Expression

- [C-SIDE](https://github.com/dmcable/spacexr) - Cell type-Specific Inference of Differential Expression in spatial transcriptomics
- [Niche-DE](https://github.com/kaishumason/NicheDE) - Niche-differential gene expression analysis identifying context-dependent cell-cell interactions
- [Vespucci](https://github.com/neurorestore/Vespucci) - Prioritize spatial regions involved in the response to an experimental perturbation in spatial transcriptomics

### Spatially Variable Genes

- [PROST](https://www.nature.com/articles/s41467-024-44835-w) - Detection of spatially variable genes
- [SpatialDE](https://github.com/Teichlab/SpatialDE) - Spatial differential expression analysis
- [SPARK-X](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0#Fig1) - Detection of spatially variable genes, best performing
- [trendsceek](https://github.com/edsgard/trendsceek) - Identification of spatial expression trends
- [nnSVG](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1) - Scalable identification of spatially variable genes using nearest-neighbor Gaussian processes

### Integration

- [PRECAST](https://www.nature.com/articles/s41467-023-35947-w) - Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data

### Cell Niches & Tissue Domains

> **(Smaller)** Cell types → Cell modules/neighborhoods → Niches/tissue domains **(Larger)**

- [BANKSY](https://www.nature.com/articles/s41588-024-01664-3) - Unified cell typing and tissue domain segmentation
- [TISSUE](https://buff.ly/49Oc0M2) - Tissue domain identification
- [CellCharter](https://www.nature.com/articles/s41588-023-01588-4) - Hierarchical niche detection
- [SpatialGLUE](https://www.nature.com/articles/s41592-024-02316-4) - Multi-omics cell niche identification
- [smoothclust](https://github.com/lmweber/smoothclust) - Spatial clustering
- [SpaTopic](https://www.science.org/doi/10.1126/sciadv.adp4942) - Spatial topic modeling
- [hdWGCNA](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) - Weighted gene correlation network analysis
- [GASTON](https://github.com/raphael-group/GASTON) - Graph-based spatial domain detection
- [SpatialMNN](https://github.com/Pixel-Dream/spatialMNN) - Identification of shared niches between slides
- [NicheCompass](https://github.com/Lotfollahi-lab/nichecompass) - End-to-end analysis of spatial multi-omics data

### Cell Distances & Neighborhood

- [CRAWDAD](https://github.com/JEFworks-Lab/CRAWDAD) - Cell relationship analysis with directional adjacency distributions
- [HoodscanR](https://www.biorxiv.org/content/10.1101/2024.03.26.586902v1) - Neighborhood analysis
- [SpicyR](https://sydneybiox.github.io/spicyR/) - Spatial analysis in R
- [MISTy](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02663-5) - Explainable multiview framework for dissecting spatial relationships from highly multiplexed data
- [SpatialCorr](https://www.biorxiv.org/content/10.1101/2022.02.04.479191v1.full) - Identifying gene sets with spatially varying correlation structure

### Spatial Trajectories

- [spaTrack](https://www.biorxiv.org/content/10.1101/2023.09.04.556175v2) - Spatial trajectory analysis
- [scSpace](https://www.biorxiv.org/content/10.1101/2022.05.07.491043v1) - Reconstruction of cell pseudo-space from single-cell RNA sequencing data
- [SOCS](https://www.biorxiv.org/content/10.1101/2025.03.19.644194v1.full.pdf) - Accurate trajectory inference in time-series spatial transcriptomics with structurally-constrained optimal transport

### Cell-Cell Communication

- [Spatia](https://www.nature.com/articles/s41592-024-02408-1#Fig2) - Spatial cell-cell interaction analysis
- [CellAgentChat](https://github.com/mcgilldinglab/CellAgentChat) - Agent-based cell communication modeling
- [SpaTalk](https://github.com/ZJUFanLab/SpaTalk) - Knowledge-graph-based cell-cell communication inference
- [SpaOTsc](https://github.com/zcang/SpaOTsc) - Inferring spatial and signaling relationships between cells
- [MISTy](https://github.com/saezlab/mistyR) - Explainable multi-view framework for dissecting intercellular signaling
- [DeepLinc](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02692-0) - De novo reconstruction of cell interaction landscapes
- [CellChat](https://github.com/jinworks/CellChat?tab=readme-ov-file) - Inferrence of cell-cell communication from multiple spatially resolved transcriptomics datasets
- [COMMOT](https://github.com/zcang/COMMOT) - Screening cell-cell communication in spatial transcriptomics via collective optimal transport
- [FlowSig](https://www.nature.com/articles/s41592-024-02380-w) - Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics
- [CellNEST](https://www.nature.com/articles/s41592-025-02721-3) - Cell–cell relay networks using attention mechanisms on spatial transcriptomics

### Metacells & Scalability

- [MetaSpot](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btae734/7919601) - Metacell analysis for spatial data
- [SEraster](https://jef.works/blog/2024/07/23/spatial-bootstrapping-with-seraster/) - Rasterization method for spatial data processing

### Subcellular Analysis

- [Sprawl](https://elifesciences.org/articles/87517) - Subcellular transcript localization
- [Bento](https://github.com/YeoLab/bento-tools) - Python toolkit for subcellular analysis of spatial transcriptomics data
- [FISHfactor](https://academic.oup.com/bioinformatics/article/39/5/btad183/7114027) - Analysis of subcellular transcript patterns
- [InSTAnT](https://www.nature.com/articles/s41467-024-49457-w) - Intracellular spatial transcript analysis

### Copy Number Variations

- [CalicoST](https://www.nature.com/articles/s41592-024-02438-9#Fig1) - CNV detection in spatial data

### Transcription Factors & Gene Regulatory Networks

- [STAN](https://www.biorxiv.org/content/10.1101/2024.06.26.600782v2.full.pdf) - Spatial transcription factor analysis

## Technical Enhancements

### Slide Alignment

- [PASTE](https://www.nature.com/articles/s41592-022-01459-6)/[PASTE2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9881963/) - Probabilistic alignment of spatial transcriptomics experiments
- [SPIRAL](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03078-6#Fig2) - Integrating and aligning spatially resolved transcriptomics data across different experiments, conditions, and technologies
- [TOAST](https://www.biorxiv.org/content/10.1101/2025.04.15.648894v1.full.pdf) - Topography Aware Optimal Transport for Alignment of Spatial Omics Data
- [SANTO](https://pmc.ncbi.nlm.nih.gov/articles/PMC11258319/pdf/41467_2024_Article_50308.pdf) - A coarse-to-fine alignment and stitching method for spatial omics

### Super Resolution

- [TESLA](https://www.cell.com/cell-systems/pdf/S2405-4712(23)00084-4.pdf) - Super resolution for 10X Visium
- [istar](https://github.com/daviddaiweizhang/istar) - Super resolution for Visium
- [BayesSPACE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763026/) - Subspot resolution
- [Spotiphy](https://www.nature.com/articles/s41592-025-02622-5) - Super resolution tool for spatial data

### Transcripts + Histology

- [ST-Net](https://github.com/bryanhe/ST-Net) - Integrating spatial gene expression and tumor morphology via deep learning
- [SpaceDIVA](https://www.biorxiv.org/content/10.1101/2025.02.19.638201v1.full.pdf) - Integration of transcript data with histological images
- [HEST](https://github.com/mahmoodlab/HEST) - Dataset for Spatial Transcriptomics and Histology Image Analysis
- [CellLENS](https://github.com/sggao/celllens/) - Cell Local Environment Neighborhood Scan
- [DeepSpot](https://www.medrxiv.org/content/10.1101/2025.02.09.25321567v2.full.pdf) - Leveraging Spatial Context for Enhanced Spatial Transcriptomics Prediction from H&E Images

## Benchmarks

- [Deconvolution benchmark](https://www.biorxiv.org/content/10.1101/2024.04.03.586404v1.full) - Comprehensive comparison
- [RCTD and Cell2location benchmark](https://www.biorxiv.org/content/10.1101/2023.03.22.533802v3.full.pdf) - Claims these are the best methods
- [Spatial clustering benchmark](https://www.nature.com/articles/s41592-024-02215-8) - Comparison of clustering methods
- [Nature Communications review](https://www.nature.com/articles/s41467-023-37168-7#Fig2) - Confirms Cell2location performance
- [Open problems benchmark](https://openproblems.bio/results/spatial_decomposition/) - Cell2location is top performer
- [Neighborhod benchmark](https://www.biorxiv.org/content/10.1101/2025.03.31.646289v2.full.pdf) - New COZI method top performer
- [Kaiko.ai FM benchmark](https://kaiko-ai.github.io/eva/main/leaderboards/) - WSI benchmark
- [Benchmarking of spatial transcriptomics platforms across six cancer types](https://www.biorxiv.org/content/10.1101/2024.05.21.593407v2) - Comprehensive platform comparison

## Datasets & Foundation Models

### Datasets
- [HISSTA](https://github.com/ercsb-sp/HISSTA/tree/v1.0)

### Foundation Models

- [scGPT-spatial](https://www.biorxiv.org/content/10.1101/2025.02.05.636714v1.full.pdf) - Language model for spatial transcriptomics
- [Phikon-v2](https://huggingface.co/owkin/phikon-v2) - Spatial biology foundation model
- [Bioptimus H-optimus-0](https://huggingface.co/bioptimus/H-optimus-0) - Biology-focused foundation model 
- [Bioptimus H-optimus-1](https://www.bioptimus.com/news/bioptimus-launches-h-optimus-1) - Latest biology-focused foundation model from Bioptimus
- [DeepCell dataset](https://exploredata.deepcell.com/cell-visualizations/9/versions/1) - CNN + human features embeddings
- [Virchow](https://huggingface.co/paige-ai/Virchow) - Foundation model for computational pathology
- [UNI](https://github.com/mahmoodlab/UNI) and [UNI2](https://huggingface.co/MahmoodLab/UNI2-h) - Universal pathology foundation models
- [CONCH](https://github.com/mahmoodlab/CONCH) - Contrastive learning for histopathology
- [GIGApath](https://github.com/prov-gigapath/prov-gigapath) - Large-scale pathology foundation model
- [OmiCLIP](https://www.nature.com/articles/s41592-025-02707-1) - A visual–omics foundation model to bridge histopathology with spatial transcriptomics

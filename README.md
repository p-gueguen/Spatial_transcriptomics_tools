# 🌌 Spatial Transcriptomics Tools

**Other resource** → [awesome_spatial_omics](https://github.com/crazyhottommy/awesome_spatial_omics)
- [Spatial-omic SIB tutorial](https://elixir-europe-training.github.io/ELIXIR-SCO-spatial-omics/presentations.html)

## Table of Contents

- [General Tools](#general-tools)
  - [Nextflow Pipelines](#nextflow--pipelines)
- [Analysis Pipeline Steps](#analysis-pipeline-steps)
  - [ROI Selection](#roi-selection)
  - [QC](#qc)
  - [Normalization](#normalization)
  - [Gene Imputation & Denoising](#gene-imputation--denoising)
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
- [LazySlide](https://github.com/rendeirolab/LazySlide) - Framework for whole slide image (WSI) analysis
- [pasta](https://robinsonlabuzh.github.io/pasta/00-home.html) - Point pattern and lattice data analysis from Robinson lab
- [rakaia](https://github.com/camlab-bioml/rakaia) - Scalable interactive visualization and analysis of spatial omics including spatial transcriptomics, in the browser ([Website](https://rakaia.io/))
- [semla](https://ludvigla.github.io/semla/index.html) - Useful tools for Spatially Resolved Transcriptomics data analysis and visualization
- [sosta](https://github.com/sgunz/sosta) - Spatial Omic Structure Analysis
- [SPATA2](https://themilolab.github.io/SPATA2/index.html) - Spatial transcriptomics analysis toolkit
- [Thor](https://doi.org/10.1038/s41467-025-62593-1) - Comprehensive platform for cell-level analysis with anti-shrinking Markov diffusion and 10 modular tools paired with Mjolnir web interface
- [VR-Omics](https://doi.org/10.1186/s13059-025-03630-6) - Free platform-agnostic software with end-to-end automated processing of multi-slice spatial transcriptomics data through biologist-friendly GUI | [Windows](https://doi.org/10.6084/m9.figshare.28259834.v3) | [MacOS](https://doi.org/10.6084/m9.figshare.28340495.v4) | [GitHub](https://github.com/Ramialison-Lab/VR-Omics)
- [CosMx-Analysis-Scratch-Space](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/) - Analysis resources and tools for CosMx SMI spatial transcriptomics | [GitHub](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space)


### Nextflow / Pipelines

- [nf-core/spatialxe](https://nf-co.re/spatialxe/dev/) - Nextflow pipeline for Xenium spatial transcriptomics analysis
- [nf-core/sopa](https://nf-co.re/sopa/dev/) - Spatial Omics Pipeline Analysis (SOPA) for processing spatial transcriptomics data
- [Allen Immunology Xenium Pipeline](https://apps.allenimmunology.org/user-documentation/data-ingest/use-the-xenium-pipeline/) - HISE platform pipeline for Xenium data processing


## Analysis Pipeline Steps

### ROI Selection

- [S2Omics](https://github.com/ddb-qiwang/S2Omics) - Designing smart spatial omics experiments with S2Omics

### QC

- [SpaceTrooper](https://htmlpreview.github.io/?https://github.com/drighelli/SpaceTrooper/blob/main/vignette/introduction.html) - Quality control for spatial transcriptomics
- [GrandQC](https://github.com/cpath-ukk/grandqc) - Comprehensive solution for quality control in digital pathology
- [SpotSweeper](https://www.nature.com/articles/s41592-025-02713-3#Fig1) - Spatially aware quality control for spatial transcriptomics
- [MerQuaCo](https://elifesciences.org/reviewed-preprints/105149) - A computational tool for quality control in image-based spatial transcriptomics

### Normalization

- [Cell volume normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w#Fig7) - Recommended for imaging-based techniques, especially with small probe lists
- [SpaNorm](https://doi.org/10.1186/s13059-025-03565-y) - First spatially-aware normalization method that concurrently models library size effects and underlying biology | [Bioconductor](https://doi.org/10.18129/B9.bioc.SpaNorm)

### Gene Imputation & Denoising

- **Note:** Gene imputation is [not recommended](https://github.com/BayraktarLab/cell2location/issues/379) for deconvolution tasks
- [SpaGE](https://github.com/tabdelaal/SpaGE) - Spatial gene expression prediction with best overall performance
- [SpaGCN](https://github.com/jianhuupenn/SpaGCN) - Spatial graph convolutional network for gene correlation analysis
- [Tangram](https://github.com/broadinstitute/Tangram) - Transcript distribution prediction and spatial mapping
- [SpaOTsc](https://github.com/zcang/SpaOTsc) - Spatial imputation via optimal transport
- [Seurat integration workflow](https://satijalab.org/seurat/articles/integration_mapping.html) - Transfer gene expression from scRNA-seq reference
- [Sprod](https://www.nature.com/articles/s41592-022-01560-w#Fig2) - Spatial denoising method
- [TISSUE](https://github.com/sunericd/TISSUE) - Transcript imputation with spatial single-cell uncertainty estimation

### Bias Correction

- [ResolVI](https://www.biorxiv.org/content/biorxiv/early/2025/01/24/2025.01.20.634005.full.pdf) - Bias correction method
- [Statial](https://sydneybiox.github.io/Statial/) - Correction of spill-over effects
- [ovrl.py](https://github.com/HiDiHlabs/ovrl.py) - A python tool to investigate vertical signal properties of imaging-based spatial transcriptomics data
- [SPLIT](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full.pdf) - SPLIT effectively resolves mixed signals and enhances cell-type purity


### Cell Segmentation

#### Imaging-based Segmentation

- [Baysor](https://github.com/kharchenkolab/Baysor) - Bayesian segmentation of spatial transcriptomics data
- [Cellpose](https://github.com/MouseLand/cellpose) - Generalist algorithm for cellular segmentation
  - [Cellpose 3](https://www.biorxiv.org/content/10.1101/2024.02.10.579780v2) - With supersampling/restoration capabilities
  - [Cellpose-SAM](https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1) - Cell and nucleus segmentation with superhuman generalization, works in 3D with various image conditions
- [DeepCell](https://github.com/vanvalenlab/deepcell-tf) - Deep learning library for single cell analysis
- [Bo Wang's method](https://www.nature.com/articles/s41592-024-02233-6#Fig3) - Better than SOTA segmentation (Nature Methods 2024)
- [Proseg](https://www.biorxiv.org/content/10.1101/2024.04.25.591218v1.full.pdf) - Probabilistic segmentation method
- [ComSeg](https://github.com/fish-quant/ComSeg) - Transcript-based point cloud segmentation
- [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) - Feature-based image segmentation
- [Xenium cell boundary](https://kb.10xgenomics.com/hc/en-us/articles/30205122555917-Boundary-stain-shows-cell-morphology-better-than-the-interior-stain) - Alternative when interior staining fails
- [Bioimage.io](https://bioimage.io/#/) - Repository of AI models for segmentation
- [ST-cellseg](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012254) - Segmentation for spatial transcriptomics
- [CelloType](https://github.com/maxpmx/CelloType) - Cell type detection and segmentation
- [SAINSC](https://onlinelibrary.wiley.com/doi/10.1002/smtd.202401123) - Segmentation for sequencing-based spatial data
- [BIDCell](https://github.com/SydneyBioX/BIDCell) - Biologically-informed deep learning for subcellular spatial transcriptomics segmentation
- [FastReseg](https://github.com/Nanostring-Biostats/FastReseg) - Using transcript locations to refine image-based cell segmentation results
- [Segger](https://www.biorxiv.org/content/10.1101/2025.03.14.643160v1) - Fast and accurate cell segmentation of imaging-based spatial transcriptomics data
- [Bering](https://doi.org/10.1038/s41467-025-60898-9) - Graph deep learning for joint noise-aware cell segmentation and molecular annotation in 2D and 3D spatial transcriptomics
- [STP](https://doi.org/10.1038/s41467-025-59782-3) - Single-cell Partition for subcellular spatially-resolved transcriptomics integrating data with nuclei-stained images
- [Deep learning-based segmentation](https://doi.org/10.1126/sciadv.adw4871) - Extensively trained nuclear and membrane segmentation models for precise transcript assignment in CosMx SMI data

**Segmentation-free methods:**
- [SSAM](https://github.com/HiDiHlabs/ssam) - Subcellular segmentation-free analysis by multidimensional mRNA density
- [Points2Regions](https://github.com/wahlby-lab/Points2Regions) - Transcript-based region identification without segmentation

#### VisiumHD Segmentation

- [Bin2Cell](https://www.biorxiv.org/content/10.1101/2024.06.19.599766v1) - Segmentation for VisiumHD data
- [ENACT](https://www.biorxiv.org/content/10.1101/2024.10.17.618905v1.full.pdf) - Enhanced accuracy for VisiumHD segmentation
- [STHD](https://www.biorxiv.org/content/10.1101/2024.06.20.599803v2) - Cell annotation for VisiumHD

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
- [ABCT](https://github.com/ercsb-sp/ABCT) - Anchor-based Cell Typer
- [STHD](https://www.biorxiv.org/content/10.1101/2024.06.20.599803v2) - Cell annotation for VisiumHD
- [STELLAR](https://www.nature.com/articles/s41592-022-01651-8) - Annotation of spatially resolved single-cell data with STELLAR
- [Vesalius](https://www.nature.com/articles/s41467-025-62782-y#Fig2) - Multi-scale and multi-context interpretable mapping of cell states across heterogeneous spatial samples
- [STALocator](https://doi.org/S2405-4712(25)00028-6) - ST-Aided Locator using deep learning to localize cells from single-cell RNA-seq data onto tissue slices
- [CMAP](https://doi.org/10.1038/s41467-025-61667-4) - Cellular Mapping of Attributes with Position, maps large-scale individual cells to precise spatial locations using divide-and-conquer strategy
- [TransST](https://doi.org/10.1186/s12859-025-06099-z) - Transfer learning framework leveraging cell-labeled information from external sources for cell-level heterogeneity inference | [GitHub](https://github.com/shuoshuoliu/TransST)
- [InSituType](https://github.com/Nanostring-Biostats/InSituType) - Cell typing for CosMx SMI spatial transcriptomics
- [HieraType](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/HieraType) - Hierarchical cell typing using RNA + protein for CosMx SMI
- [CosMx-Cell-Profiles](https://github.com/Nanostring-Biostats/CosMx-Cell-Profiles) - Collection of reference datasets for CosMx SMI

### Cell Deconvolution

- [RCTD](https://github.com/dmcable/spacexr) - Robust cell type decomposition
- [Cell2location](https://github.com/BayraktarLab/cell2location) - Mapping scRNA-seq to spatial data
- [SPOTlight](https://github.com/MarcElosua/SPOTlight) - Seeded NMF regression to deconvolute spatial spots
- [CARD](https://github.com/YingMa0107/CARD) - Spatially informed cell-type deconvolution
- [FlashDeconv](https://github.com/cafferychen777/flashdeconv) - High-performance deconvolution using randomized sketching, achieves linear O(N) scaling for Visium HD and other high-resolution platforms

### Differential Expression

- [C-SIDE](https://github.com/dmcable/spacexr) - Cell type-Specific Inference of Differential Expression in spatial transcriptomics
- [Niche-DE](https://github.com/kaishumason/NicheDE) - Niche-differential gene expression analysis identifying context-dependent cell-cell interactions
- [smiDE](https://www.biorxiv.org/content/10.1101/2024.08.02.606405v2) - Spatial differential expression method
- [spatialGE](https://github.com/FridleyLab/spatialGE) - Spatial gene expression analysis
- [Vespucci](https://github.com/neurorestore/Vespucci) - Prioritize spatial regions involved in the response to an experimental perturbation in spatial transcriptomics

### Spatially Variable Genes

- [PROST](https://www.nature.com/articles/s41467-024-44835-w) - Detection of spatially variable genes
- [SpatialDE](https://github.com/Teichlab/SpatialDE) - Spatial differential expression analysis
- [SPARK-X](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0#Fig1) - Detection of spatially variable genes, best performing
- [Hotspot](https://github.com/YosefLab/Hotspot) - Identify informative gene modules with lowest false positive rate
- [SOMDE](https://github.com/XuegongLab/somde) - Self-organizing map for spatially variable gene detection with optimization
- [trendsceek](https://github.com/edsgard/trendsceek) - Identification of spatial expression trends
- [nnSVG](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v1) - Scalable identification of spatially variable genes using nearest-neighbor Gaussian processes

### Integration

- [PRECAST](https://www.nature.com/articles/s41467-023-35947-w) - Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data
- [MISO](https://doi.org/10.1038/s41592-024-02574-2) - MultI-modal Spatial Omics for versatile feature extraction and clustering integrating multiple modalities including gene expression, protein, epigenetics, metabolomics, and histology
- [SIMO](https://doi.org/10.1038/s41467-025-56523-4) - Spatial Integration of Multi-Omics through probabilistic alignment integrating spatial transcriptomics with multiple single-cell modalities
- [GSI](https://doi.org/10.1093/bioinformatics/btaf350) - Gene Spatial Integration using deep learning with representation learning to extract spatial distribution of genes | [GitHub](https://github.com/Riandanis/Spatial_Integration_GSI) | [Zenodo](https://doi.org/10.5281/zenodo.15165223)
- [SPACE-seq](https://doi.org/10.1073/pnas.2424070122) - Unified molecular approach for spatial multiomics enabling simultaneous analysis of chromatin accessibility, mitochondrial DNA mutations, and gene expression on standard 10× Genomics Visium CytAssist platform

### Cell Niches & Tissue Domains

> (**Smaller**) Cell types → Cell modules/neighborhoods → Niches/tissue domains (**Larger**)

- [BANKSY](https://www.nature.com/articles/s41588-024-01664-3) - Unified cell typing and tissue domain segmentation
- [CellCharter](https://www.nature.com/articles/s41588-023-01588-4) - Hierarchical niche detection
- [SpatialGLUE](https://www.nature.com/articles/s41592-024-02316-4) - Multi-omics cell niche identification
- [smoothclust](https://github.com/lmweber/smoothclust) - Spatial clustering
- [SpaTopic](https://www.science.org/doi/10.1126/sciadv.adp4942) - Spatial topic modeling
- [hdWGCNA](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) - Weighted gene correlation network analysis
- [GASTON](https://github.com/raphael-group/GASTON) - Graph-based spatial domain detection
- [SpatialMNN](https://github.com/Pixel-Dream/spatialMNN) - Identification of shared niches between slides
- [NicheCompass](https://github.com/Lotfollahi-lab/nichecompass) - End-to-end analysis of spatial multi-omics data
- [Proust](https://genome.cshlp.org/content/35/7/1621) - Spatial domain detection using contrastive self-supervised learning for spatial multi-omics technologies (multi-modal domains)
- [STAMP](https://doi.org/10.1038/s41592-024-02463-8) - Spatial Transcriptomics Analysis with topic Modeling, provides interpretable dimension reduction through deep generative modeling discovering tissue domains and cellular communication patterns
- [DeepGFT](https://doi.org/10.1186/s13059-025-03631-5) - Combines deep learning with graph Fourier transform for spatial domain identification | [GitHub](https://github.com/jxLiu-bio/DeepGFT) | [Zenodo](https://doi.org/10.5281/zenodo.15073243)
- [novae](https://github.com/MICS-Lab/novae) - Deep learning framework for spatial domain detection and tissue organization analysis

### Cell Distances & Neighborhood

- [CRAWDAD](https://github.com/JEFworks-Lab/CRAWDAD) - Cell relationship analysis with directional adjacency distributions
- [HoodscanR](https://www.biorxiv.org/content/10.1101/2024.03.26.586902v1) - Neighborhood analysis
- [SpicyR](https://sydneybiox.github.io/spicyR/) - Spatial analysis in R
- [MISTy](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02663-5) - Explainable multiview framework for dissecting spatial relationships from highly multiplexed data
- [SpatialCorr](https://www.biorxiv.org/content/10.1101/2022.02.04.479191v1.full) - Identifying gene sets with spatially varying correlation structure
- [CatsCradle](https://github.com/AnnaLaddach/CatsCradle) - Spatial analysis framework for tissue neighbourhoods

### Spatial Trajectories

- [spaTrack](https://www.biorxiv.org/content/10.1101/2023.09.04.556175v2) - Spatial trajectory analysis
- [scSpace](https://www.biorxiv.org/content/10.1101/2022.05.07.491043v1) - Reconstruction of cell pseudo-space from single-cell RNA sequencing data
- [SOCS](https://www.biorxiv.org/content/10.1101/2025.03.19.644194v1.full.pdf) - Accurate trajectory inference in time-series spatial transcriptomics with structurally-constrained optimal transport
- [STORIES](https://doi.org/10.1038/s41592-025-02855-4) - Spatiotemporal Reconstruction Using Optimal Transport for cell trajectory inference from spatial transcriptomics profiled at multiple time points

### Cell-Cell Communication

- [Spatia](https://www.nature.com/articles/s41592-024-02408-1#Fig2) - Spatial cell-cell interaction analysis
- [CellAgentChat](https://github.com/mcgilldinglab/CellAgentChat) - Agent-based cell communication modeling
- [SpaTalk](https://github.com/ZJUFanLab/SpaTalk) - Knowledge-graph-based cell-cell communication inference
- [SpaOTsc](https://github.com/zcang/SpaOTsc) - Inferring spatial and signaling relationships between cells
- [MISTy](https://github.com/saezlab/mistyR) - Explainable multi-view framework for dissecting intercellular signaling
- [DeepLinc](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02692-0) - De novo reconstruction of cell interaction landscapes
- [CellChat](https://github.com/jinworks/CellChat?tab=readme-ov-file) - Inferrence of cell-cell communication from multiple spatially resolved transcriptomics datasets
- [COMMOT](https://github.com/zcang/COMMOT) - Screening cell-cell communication in spatial transcriptomics via collective optimal transport
- [NicheNet](https://github.com/saeyslab/nichenetr) - Linking ligands to downstream target gene regulation
- [DeepTalk](https://github.com/JiangBioLab/DeepTalk) - Single-cell resolution cell-cell communication using deep learning
- [CellNEST](https://github.com/schwartzlab-methods/CellNEST) - Cell–cell relay networks using attention mechanisms on spatial transcriptomics
- [FlowSig](https://www.nature.com/articles/s41592-024-02380-w) - Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics

### Metacells & Scalability

- [MetaSpot](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btae734/7919601) - Metacell analysis for spatial data
- [SEraster](https://jef.works/blog/2024/07/23/spatial-bootstrapping-with-seraster/) - Rasterization method for spatial data processing

### Subcellular Analysis

- [Sprawl](https://elifesciences.org/articles/87517) - Subcellular transcript localization
- [Bento](https://github.com/YeoLab/bento-tools) - Python toolkit for subcellular analysis of spatial transcriptomics data
- [FISHfactor](https://academic.oup.com/bioinformatics/article/39/5/btad183/7114027) - Analysis of subcellular transcript patterns
- [InSTAnT](https://www.nature.com/articles/s41467-024-49457-w) - Intracellular spatial transcript analysis
- [troutpy](https://github.com/theislab/troutpy) - Analysis of transcripts outside segmented cells in spatial transcriptomics data

### Copy Number Variations

- [CalicoST](https://www.nature.com/articles/s41592-024-02438-9#Fig1) - CNV detection in spatial data
- [inSituCNV](https://github.com/Moldia/InSituCNV) - Inference of Copy Number Variations in Image-Based Spatial Transcriptomics

### Transcription Factors & Gene Regulatory Networks

- [STAN](https://www.biorxiv.org/content/10.1101/2024.06.26.600782v2.full.pdf) - Spatial transcription factor analysis

## Technical Enhancements

### Slide Alignment

- [PASTE](https://www.nature.com/articles/s41592-022-01459-6)/[PASTE2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9881963/) - Probabilistic alignment of spatial transcriptomics experiments
- [SPIRAL](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03078-6#Fig2) - Integrating and aligning spatially resolved transcriptomics data across different experiments, conditions, and technologies
- [TOAST](https://www.biorxiv.org/content/10.1101/2025.04.15.648894v1.full.pdf) - Topography Aware Optimal Transport for Alignment of Spatial Omics Data
- [STalign](https://jef.works/STalign/notebooks/xenium-xenium-alignment.html) - Alignment of spatial transcriptomics data using diffeomorphic metric mapping
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
- [SpotWhisperer](https://medical-epigenomics.org/papers/spotwhisperer/#citation) - Molecularly informed analysis of histopathology images using natural language
- [STPath](https://github.com/Graph-and-Geometric-Learning/STPath) - A Generative Foundation Model for Integrating Spatial Transcriptomics and Whole Slide Images
- [AESTETIK](https://github.com/ratschlab/aestetik) - AutoEncoder for Spatial Transcriptomics Expression with Topology and Image Knowledge

## Benchmarks

- [Deconvolution benchmark](https://www.biorxiv.org/content/10.1101/2024.04.03.586404v1.full) - Comprehensive comparison
- [RCTD and Cell2location benchmark](https://www.biorxiv.org/content/10.1101/2023.03.22.533802v3.full.pdf) - Claims these are the best methods
- [Spatial clustering benchmark](https://www.nature.com/articles/s41592-024-02215-8) - Comparison of clustering methods
- [Spatialbench](https://github.com/latchbio/spatialbench) - Benchmark for evaluating AI agents on spatial biology analysis tasks
- [Nature Communications review](https://www.nature.com/articles/s41467-023-37168-7#Fig2) - Confirms Cell2location performance
- [Open problems benchmark](https://openproblems.bio/results/spatial_decomposition/) - Cell2location is top performer
- [Neighborhod benchmark](https://www.biorxiv.org/content/10.1101/2025.03.31.646289v2.full.pdf) - New COZI method top performer
- [Kaiko.ai FM benchmark EVA](https://kaiko-ai.github.io/eva/main/leaderboards/) - WSI benchmark
- [Benchmarking of spatial transcriptomics platforms across six cancer types](https://www.biorxiv.org/content/10.1101/2024.05.21.593407v2) - Comprehensive platform comparison
- [PathBench](https://birkhoffkiki.github.io/PathBench/) - Pathology benchmark
- [SPATCH Benchmark - 2025](https://www.nature.com/articles/s41467-025-64292-3) - Showing Xenium performs best
- [Thunder](https://github.com/MICS-Lab/thunder) - Pathology benchmark

## Datasets & Foundation Models

### Datasets

- [HISSTA](https://github.com/ercsb-sp/HISSTA/tree/v1.0) - Histopathology spatial transcriptomics dataset
- [STOmicsDB](https://db.cngb.org/stomics/) - Spatial transcriptomics database

### Foundation Models

- [KRONOS](https://github.com/mahmoodlab/KRONOS) - Foundation Model for Multiplex Spatial Proteomic Images
- [scGPT-spatial](https://www.biorxiv.org/content/10.1101/2025.02.05.636714v1.full.pdf) - Language model for spatial transcriptomics
- [Phikon-v2](https://huggingface.co/owkin/phikon-v2) - Spatial biology foundation model
- [Bioptimus H-optimus-0](https://huggingface.co/bioptimus/H-optimus-0) - Biology-focused foundation model
- [Bioptimus H-optimus-1](https://www.bioptimus.com/news/bioptimus-launches-h-optimus-1) - Latest biology-focused foundation model from Bioptimus
- [DeepCell dataset](https://exploredata.deepcell.com/cell-visualizations/9/versions/1) - CNN + human features embeddings
- [TITAN](https://www.nature.com/articles/s41591-025-03982-3) - A multimodal whole-slide foundation model for pathology
- [Virchow](https://huggingface.co/paige-ai/Virchow) - Foundation model for computational pathology
- [UNI](https://github.com/mahmoodlab/UNI) and [UNI2](https://huggingface.co/MahmoodLab/UNI2-h) - Universal pathology foundation models
- [CONCH](https://github.com/mahmoodlab/CONCH) - Contrastive learning for histopathology
- [GIGApath](https://github.com/prov-gigapath/prov-gigapath) - Large-scale pathology foundation model
- [OmiCLIP](https://www.nature.com/articles/s41592-025-02707-1) - A visual–omics foundation model to bridge histopathology with spatial transcriptomics
- [Nicheformer](https://doi.org/10.1038/s41592-025-02814-z) - Transformer-based foundation model pretrained on SpatialCorpus-110M containing over 110 million cells for spatial composition and label prediction | [GitHub](https://github.com/theislab/nicheformer-data/)

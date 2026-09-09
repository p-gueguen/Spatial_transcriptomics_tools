# 🌌 Spatial Transcriptomics Tools

**Other resource** → [awesome_spatial_omics](https://github.com/crazyhottommy/awesome_spatial_omics)
- [Spatial-omic SIB tutorial](https://elixir-europe-training.github.io/ELIXIR-SCO-spatial-omics/presentations.html)

## Table of Contents

- [General Tools](#general-tools)
  - [Nextflow Pipelines](#nextflow--pipelines)
- [Viewers & Interactive Annotation](#viewers--interactive-annotation)
- [Probe & Panel Design](#probe--panel-design)
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
  - [Isoform Analysis](#isoform-analysis)
  - [Transcription Factors & Gene Regulatory Networks](#transcription-factors--gene-regulatory-networks)
- [Technical Enhancements](#technical-enhancements)
  - [Slide Alignment](#slide-alignment)
  - [Super Resolution](#super-resolution)
  - [Transcripts + Histology](#transcripts--histology)
- [Benchmarks](#benchmarks)
  - [Simulators & Ground Truth](#simulators--ground-truth)
- [Datasets & Foundation Models](#datasets--foundation-models)


## General Tools

- [Best practices Bioconductor](https://lmweber.org/OSTA/) - [R] - Principles for statistical analysis of spatial transcriptomics data
- [squidpy](https://github.com/scverse/squidpy) - [Python] - Spatial single cell analysis toolkit from scverse
- [Giotto](https://github.com/giotto-suite/Giotto) - [R/Python] - Comprehensive spatial data analysis suite
- [Vitessce](https://github.com/vitessce/vitessce) - [JavaScript] - Visual integration tool for exploration of spatial single cell experiments
- [Voyager](https://github.com/pachterlab/voyager) - [R] - Spatial transcriptomics visualization from Pachter lab
- [BASS](https://github.com/zhengli09/BASS) - [R] - Multiple sample analysis
- [SpaVAE](https://github.com/ttgump/spaVAE) - [Python] - All-purpose tool for dimension reduction, visualization, clustering, batch integration, denoising, differential expression, spatial interpolation, and resolution enhancement | [Implementations](https://github.com/hrlblab/computer_vision_spatial_omics)
- [sopa](https://github.com/prism-oncology/sopa) - [Python] - Spatial omics processing and analysis
- [SpatialAgent](https://github.com/Genentech/SpatialAgent) - [Python] - An autonomous AI agent for spatial biology
- [ChatSpatial](https://github.com/cafferychen777/ChatSpatial) - [Python] - MCP server enabling spatial transcriptomics analysis via natural language, integrating 60+ methods including SpaGCN, Cell2location, LIANA+, CellRank for Visium, Xenium, MERFISH | [Paper](https://doi.org/10.64898/2026.02.26.708361) | [Docs](https://cafferychen777.github.io/ChatSpatial/) | [PyPI](https://pypi.org/project/chatspatial/)
- [STAgent](https://github.com/LiuLab-Bioelectronics-Harvard/STAgent) - [Python] - Autonomous multimodal LLM agent (Claude/GPT/Gemini) for end-to-end ST analysis of .h5ad data, with code generation, visual reasoning over tissue images, literature retrieval, and report synthesis via a Streamlit interface | [Paper](https://doi.org/10.1101/2025.04.01.646731)
- [LazySlide](https://github.com/rendeirolab/LazySlide) - [Python] - Framework for whole slide image (WSI) analysis
- [pasta](https://robinsonlabuzh.github.io/pasta/00-home.html) - [R] - Point pattern and lattice data analysis from Robinson lab
- [rakaia](https://github.com/camlab-bioml/rakaia) - [JavaScript] - Scalable interactive visualization and analysis of spatial omics including spatial transcriptomics, in the browser ([Website](https://rakaia.io/))
- [semla](https://github.com/spatial-research/semla) - [R] - Useful tools for Spatially Resolved Transcriptomics data analysis and visualization | [Docs](https://spatial-research.github.io/semla/)
- [sosta](https://github.com/sgunz/sosta) - [Python] - Spatial Omic Structure Analysis
- [SPATA2](https://themilolab.github.io/SPATA2/index.html) - [R] - Spatial transcriptomics analysis toolkit
- [spatial-omics-tutorials](https://github.com/pnucolab/spatial-omics-tutorials) - [Python/R] - Tutorials and best-practices for spatial omics data analysis from BIML 2025
- [MGC-BioSB-Spatial-Omics-Analysis-2025](https://github.com/mahfouzlab/MGC-BioSB-Spatial-Omics-Analysis-2025) - [Python] - Workshop materials for spatial omics analysis
- [SMINT](https://github.com/JurgenKriel/SMINT) - [Python] - Spatial Multi-Omics Integration Toolkit for transcriptomics and metabolomics | [Docs](https://jurgenkriel.github.io/SMINT/)
- [Thor](https://github.com/GuangyuWangLab2021/Thor) - [Python] - Comprehensive platform for cell-level analysis with anti-shrinking Markov diffusion and 10 modular tools paired with Mjolnir web interface
- [VR-Omics](https://doi.org/10.1186/s13059-025-03630-6) - [GUI] - Free platform-agnostic software with end-to-end automated processing of multi-slice spatial transcriptomics data through biologist-friendly GUI | [Windows](https://doi.org/10.6084/m9.figshare.28259834.v3) | [MacOS](https://doi.org/10.6084/m9.figshare.28340495.v4) | [GitHub](https://github.com/Ramialison-Lab/VR-Omics)
- [CosMx-Analysis-Scratch-Space](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/) - [R/Python] - Analysis resources and tools for CosMx SMI spatial transcriptomics | [GitHub](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space)
- [SpaceSequest](https://github.com/interactivereport/SpaceSequest) - [R] - Unified pipeline for analysis, visualization, and publication of spatial transcriptomics data from Visium, Visium HD, Xenium, GeoMx, and CosMx | [Tutorial](https://interactivereport.github.io/SpaceSequest/tutorial/docs/index.html)
- [VST-DAVis](https://github.com/GudaLab/VST-DAVis) - [R Shiny] - Browser-based GUI for end-to-end Visium HD spatial transcriptomics analysis including QC, clustering, cell annotation, pathway enrichment, CellChat, and trajectory analysis
- [BrainConnect](https://github.com/CPenglab/BrainConnect) - [Python] - Integrative analysis of mouse brain connectivity and whole-brain spatial transcriptomics using LSTM networks to predict connectivity strength from regional gene expression
- [stLearn](https://github.com/BiomedicalMachineLearning/stLearn) - [Python] - Combines expression, spatial distance and tissue morphology for clustering, pseudo-time-space trajectories and cell-cell interaction
- [gsMap](https://github.com/JianYang-Lab/gsMap) - [Python] - Integrates GWAS summary statistics with spatial transcriptomics to map cells and spatial regions associated with human complex traits
- [VoltRon](https://github.com/BIMSBbioinfo/VoltRon) - [R] - Multi-resolution, multi-omic spatial toolbox with built-in image registration across assays and resolutions
- [spatialLIBD](https://github.com/LieberInstitute/spatialLIBD) - [R] - Bioconductor package and Shiny app for interactive visualization of Visium and spot-based spatial data


### Nextflow / Pipelines

- [nf-core/spatialaxe](https://nf-co.re/spatialaxe/) - [Nextflow] - Nextflow pipeline for Xenium and Artera spatial transcriptomics analysis (renamed from spatialxe)
- [nf-core/sopa](https://nf-co.re/sopa/dev/) - [Nextflow] - Spatial Omics Pipeline Analysis (SOPA) for processing spatial transcriptomics data
- [Allen Immunology Xenium Pipeline](https://apps.allenimmunology.org/user-documentation/data-ingest/use-the-xenium-pipeline/) - [Web] - HISE platform pipeline for Xenium data processing
- [SCALPEL](https://github.com/AllenInstitute/Spatial-Transcriptomics-Processing-Pipeline) - [Nextflow] - Allen Institute pipeline for large-scale ST atlas construction with 3D segmentation, doublet detection (SOLO), MapMyCells label transfer, and CCF registration
- [NNclinSSOAP](https://github.com/NovoNordisk-OpenSource/nnclinssoap) - [Nextflow/R] - GxP-ready clinical pipeline for 10x Xenium spatial transcriptomics and scRNA-seq analysis with Docker/Apptainer containerization
- [nf-core/spatialvi](https://github.com/nf-core/spatialvi) - [Nextflow] - nf-core pipeline for Visium data: spot counts, spatial coordinates and image data through QC, normalization and clustering
- [WebAtlas](https://github.com/haniffalab/webatlas-pipeline) - [Python/Nextflow] - Converts h5ad, SpaceRanger, Xenium and MERSCOPE output into Zarr + OME-TIFF for browser-based Vitessce atlases | [Paper](https://doi.org/10.1038/s41592-024-02371-x)


## Viewers & Interactive Annotation

- [QuPath](https://github.com/qupath/qupath) - [Java] - De-facto open-source whole-slide image viewer and annotator: ROI drawing, cell detection, stain deconvolution, scripting, Xenium/Visium image overlays | [Paper](https://doi.org/10.1038/s41598-017-17204-5)
- [InstanSeg](https://github.com/instanseg/InstanSeg) - [Python] - Embedding-based nucleus and cell segmentation, shipped as the official QuPath extension; faster and more accurate than Cellpose/StarDist on six public datasets | [Paper](https://arxiv.org/abs/2408.15954)
- [napari-spatialdata](https://github.com/scverse/napari-spatialdata) - [Python] - Desktop napari viewer and annotator for spatial data: multi-scale images, transcripts and polygons with a live scanpy round-trip
- [TissUUmaps](https://github.com/TissUUmaps/TissUUmaps) - [JavaScript/Python] - GPU-accelerated browser viewer for 10^7+ transcripts over multi-resolution tissue images; Jupyter-embeddable, with in-situ decoding QC modules | [Paper](https://doi.org/10.1016/j.heliyon.2023.e15306)

## Probe & Panel Design

### Panel Selection

- [Spapros](https://github.com/theislab/spapros) - [Python] - Gene panel selection for targeted spatial transcriptomics (Visium HD, Xenium, MERFISH, etc.) from an scRNA-seq reference, optimizing for cell-type ID and within-cell-type variation
- [scGPD](https://github.com/TinaGuo/scGPD) - [Python] - Deep-learning gene panel design for targeted ST; a two-stage, gene-gene correlation-aware binary gating mechanism selects compact, nonredundant marker panels from an scRNA-seq reference | [Paper](https://doi.org/10.1093/bib/bbag160)
- [ReconST](https://github.com/haoranlustat/ReconST) - [Python] - Gated-autoencoder panel design that selects the gene subset best reconstructing the whole transcriptome (rather than relying on known markers) for MERFISH/seqFISH/Xenium/MERSCOPE | [Paper](https://www.biorxiv.org/content/10.1101/2025.10.08.681071v1)

### Probe Sequence Design (Visium / Visium HD / FFPE)

- [gene2probe](https://github.com/Teichlab/gene2probe) - [Python] - Pipeline for designing custom probes against human genes, tuned to VisiumHD / Visium FFPE constraints
- [ProbeST](https://github.com/giacomellolab/ProbeST) - [Snakemake] - Custom Visium-compatible probe design pipeline for any species (including prokaryotes and host–pathogen)

### Off-target QC

- [Off-Target Probe Tracker](https://github.com/JEFworks-Lab/off-target-probe-tracker) - [Python] - Pipeline to predict off-target binding from probe sequences, useful for auditing Xenium and other targeted ST panels

## Analysis Pipeline Steps

### ROI Selection

- [S2Omics](https://github.com/ddb-qiwang/S2Omics) - [Python] - Designing smart spatial omics experiments with S2Omics
- [SOFisher](https://github.com/Heathcliff-Ng/SOFisher) - [Python] - Reinforcement-learning agent that chooses the next field-of-view position from previously sampled FOVs, cutting the acquisition needed to reach regions of interest | [Paper](https://doi.org/10.1038/s41467-026-73404-6)

### QC

- [SpaceTrooper](https://github.com/drighelli/SpaceTrooper) - [R] - Quality control for spatial transcriptomics
- [GrandQC](https://github.com/cpath-ukk/grandqc) - [Python] - Comprehensive solution for quality control in digital pathology
- [SpotSweeper](https://github.com/MicTott/SpotSweeper) - [R] - Spatially aware quality control for spatial transcriptomics
- [MerQuaCo](https://github.com/AllenInstitute/merquaco) - [Python] - A computational tool for quality control in image-based spatial transcriptomics
- [SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM) - [R] - Standardized QC-metric suite for imaging-based ST (Xenium/CosMx/MERSCOPE): transcripts/cell, global FDR, signal-to-noise, Moran's I and more; the software arm of the Spatial Touchstone reproducibility framework | [Paper](https://doi.org/10.1038/s41587-025-02811-9)
- [SpatialArtifacts](https://bioconductor.org/packages/SpatialArtifacts) - [R] - Identification and classification of spatial artifacts (edge and interior) in Visium and Visium HD data via a two-step outlier-detection plus image-processing workflow | [GitHub](https://github.com/CambridgeCat13/SpatialArtifacts)
- [destriping-GLM](https://github.com/paolamalsot/destriping-GLM) - [Python] - Removes Visium HD stripe artifacts by GLM modelling of per-bin nuclear counts | [Preprint](https://doi.org/10.64898/2026.05.04.722591)

### Normalization

- [Cell volume normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w#Fig7) - [R] - Recommended for imaging-based techniques, especially with small probe lists
- [SpaNorm](https://github.com/bhuvad/SpaNorm) - [R] - First spatially-aware normalization method that concurrently models library size effects and underlying biology | [Bioconductor](https://doi.org/10.18129/B9.bioc.SpaNorm)
- [TranspaceR](https://github.com/BOSTLAB/TranspaceR) - [R] - Statistical framework for imaging-based (Xenium/CosMx) and Visium HD ST: novel QC, scalable spatially-variable-gene selection, and data-driven optimal normalization plus low-dimensional embedding | [Paper](https://doi.org/10.1186/s44330-026-00055-w)

### Gene Imputation & Denoising

- **Note:** Gene imputation is [not recommended](https://github.com/BayraktarLab/cell2location/issues/379) for deconvolution tasks
- [SpaGE](https://github.com/tabdelaal/SpaGE) - [Python] - Spatial gene expression prediction with best overall performance
- [SpaGCN](https://github.com/jianhuupenn/SpaGCN) - [Python] - Spatial graph convolutional network for gene correlation analysis
- [Tangram](https://github.com/broadinstitute/Tangram) - [Python] - Transcript distribution prediction and spatial mapping
- [SpaOTsc](https://github.com/zcang/SpaOTsc) - [Python] - Spatial imputation via optimal transport
- [Seurat integration workflow](https://satijalab.org/seurat/articles/integration_mapping.html) - [R] - Transfer gene expression from scRNA-seq reference
- [Sprod](https://github.com/yunguan-wang/SPROD) - [Python] - Spatial denoising method
- [TISSUE](https://github.com/sunericd/TISSUE) - [Python] - Transcript imputation with spatial single-cell uncertainty estimation
- [SpaIM](https://github.com/QSong-github/SpaIM) - [Python] - Single-cell Spatial Transcriptomics Imputation via Style Transfer
- [Cellpin](https://github.com/Saur-Lab/cellpin) - [Python] - Lightweight probabilistic model for reference-based imputation and denoising of spatial transcriptomes from an scRNA-seq reference, enabling transcriptome-wide imputation and atlas-to-spatial label transfer for targeted-panel and full-transcriptome data | [Paper](https://doi.org/10.64898/2026.06.02.729566)

### Bias Correction

- [ResolVI](https://github.com/scverse/scvi-tools) - [Python] - Bias correction method
- [Statial](https://sydneybiox.github.io/Statial/) - [R] - Correction of spill-over effects
- [ovrl.py](https://github.com/HiDiHlabs/ovrl.py) - [Python] - A python tool to investigate vertical signal properties of imaging-based spatial transcriptomics data
- [SPLIT](https://github.com/bdsc-tds/SPLIT) - [R] - SPLIT effectively resolves mixed signals and enhances cell-type purity
- [cellAdmix](https://github.com/kharchenkolab/cellAdmix) - [R] - From Kharchenko lab - Evaluating and correcting cell admixtures in imaging-based spatial transcriptomics data.
- [DenoIST](https://github.com/aaronkwc/DenoIST) - [R] - Denoising Image-based Spatial Transcriptomics data
- [MisTIC](https://github.com/yunguan-wang/MisTic-Wanglab) - [Python] - A probabilistic model for correcting mis-assigned transcripts due to cell segmentation errors
- [TRACER](https://github.com/imlong4real/TRACER) - [Python] - Tissue Reconstruction via Associative Clique Extraction and Relation-mapping
- [RESCUE](https://github.com/brunoyjlee/RESCUE) - [R] - Negative-selection method that partitions ST expression into reference-explained "canonical" and sparse "idiosyncratic" components, recovering biological signal (fragile cell types, neurites, extracellular transcripts) lost by reference-based deconvolution/segmentation | [Paper](https://doi.org/10.1038/s41467-026-71720-5)
- [XeniumClean](https://github.com/rachael-z/XeniumClean) - [R] - Spatial neighbour-aware transcript cleanup for imaging-based ST (Xenium, MERSCOPE, CosMx, Atera): removes biologically implausible transcripts by combining a single-cell RNA-seq reference with spatial neighbourhood information to erase genes that cannot plausibly originate from a cell's own type
- [SPARKLE](https://github.com/WangShuai-3/SPARKLE) - [Python] - Evidence-constrained kernel estimator that removes local ambient RNA leakage in imaging-based ST | [Preprint](https://doi.org/10.64898/2026.08.12.744394)
- [CLEAR-ST](https://github.com/holab-hku/CLEAR-ST) - [Python] - Physics-informed probabilistic decontamination modelling mRNA lateral diffusion between neighbouring cells | [Preprint](https://doi.org/10.64898/2026.08.13.744615)
- [DeSpotX](https://github.com/Gentles-lab/DeSpotX) - [Python] - Identifiability-based deep generative decontamination of single-cell-resolution ST | [Preprint](https://doi.org/10.64898/2026.05.12.724704)

### Cell Segmentation

#### Imaging-based Segmentation

- [Baysor](https://github.com/kharchenkolab/Baysor) - [Julia] - Bayesian segmentation of spatial transcriptomics data
- [Cellpose](https://github.com/MouseLand/cellpose) - [Python] - Generalist algorithm for cellular segmentation
  - [Cellpose 3](https://github.com/MouseLand/cellpose) - With supersampling/restoration capabilities
  - [Cellpose-SAM](https://github.com/MouseLand/cellpose) - Cell and nucleus segmentation with superhuman generalization, works in 3D with various image conditions
- [DeepCell](https://github.com/vanvalenlab/deepcell-tf) - [Python] - Deep learning library for single cell analysis
- [Bo Wang's method](https://www.nature.com/articles/s41592-024-02233-6#Fig3) - [Python] - Better than SOTA segmentation (Nature Methods 2024)
- [Proseg](https://github.com/dcjones/proseg) - [Rust] - Probabilistic segmentation method
- [ComSeg](https://github.com/fish-quant/ComSeg) - [Python] - Transcript-based point cloud segmentation
- [FICTURE](https://github.com/seqscope/ficture) - [Python] - Feature-based image segmentation
- [Xenium cell boundary](https://kb.10xgenomics.com/hc/en-us/articles/30205122555917-Boundary-stain-shows-cell-morphology-better-than-the-interior-stain) - [Web] - Alternative when interior staining fails
- [Bioimage.io](https://bioimage.io/#/) - [Web] - Repository of AI models for segmentation
- [ST-cellseg](https://github.com/xjtulyc/ST-Cellseg) - [Python] - Segmentation for spatial transcriptomics
- [CelloType](https://github.com/maxpmx/CelloType) - [Python] - Cell type detection and segmentation
- [SAINSC](https://github.com/HiDiHlabs/sainsc) - [Python] - Segmentation for sequencing-based spatial data
- [BIDCell](https://github.com/SydneyBioX/BIDCell) - [Python] - Biologically-informed deep learning for subcellular spatial transcriptomics segmentation
- [FastReseg](https://github.com/Nanostring-Biostats/FastReseg) - [R] - Using transcript locations to refine image-based cell segmentation results
- [Segger](https://github.com/dpeerlab/segger) - [Python] - Fast and accurate cell segmentation of imaging-based spatial transcriptomics data
- [Bering](https://github.com/jian-shu-lab/Bering) - [Python] - Graph deep learning for joint noise-aware cell segmentation and molecular annotation in 2D and 3D spatial transcriptomics
- [STP](https://github.com/leihouyeung/STP) - [Python] - Single-cell Partition for subcellular spatially-resolved transcriptomics integrating data with nuclei-stained images
- [Deep learning-based segmentation](https://doi.org/10.1126/sciadv.adw4871) - [Python] - Extensively trained nuclear and membrane segmentation models for precise transcript assignment in CosMx SMI data
- [CellSAM](https://github.com/vanvalenlab/cellSAM) - [Python] - Foundation model for cell segmentation achieving state-of-the-art performance across cellular targets (bacteria, tissue, yeast, cell culture) and imaging modalities (brightfield, fluorescence, phase, multiplexed) | [Paper](https://www.nature.com/articles/s41592-025-02879-w) | [Web App](https://cellsam.deepcell.org)
- [DISSECT](https://github.com/zenglab-pku/DISSECT) - [Python] - Diffusion-based cell segmentation combining cytological image segmentation with transcriptome-guided boundary refinement for Xenium, CosMx, and Stereo-seq | [PyPI](https://pypi.org/project/dissect-st/)

**Segmentation-free methods:**
- [SSAM](https://github.com/HiDiHlabs/ssam) - [Python] - Subcellular segmentation-free analysis by multidimensional mRNA density (repo archived 2024-10, reference implementation)
- [Points2Regions](https://github.com/wahlby-lab/Points2Regions) - [Python] - Transcript-based region identification without segmentation
- [Cellist](https://github.com/wanglabtongji/Cellist) - [Python] - Multi-modal segmentation combining image and expression signal, benchmarked across Stereo-seq, Seq-Scope, seqFISH+, STARmap and Xenium | [Paper](https://doi.org/10.1038/s41588-026-02610-1)
- [RNA2seg](https://github.com/fish-quant/rna2seg) - [Python] - Generalist segmentation model trained on 4M MERFISH/CosMx cells, fusing RNA point clouds with membrane and nuclear stains; zero-shot capable | [Paper](https://doi.org/10.1186/s13059-025-03908-9)

#### VisiumHD Segmentation

- [Bin2Cell](https://github.com/Teichlab/bin2cell) - [Python] - Segmentation for VisiumHD data
- [ENACT](https://github.com/Sanofi-Public/enact-pipeline) - [Python] - Enhanced accuracy for VisiumHD segmentation
- [STHD](https://github.com/yi-zhang/STHD) - [Python] - Cell annotation for VisiumHD
- [Segmentation and annotation pipeline for VisiumHD with Proseg and Novae](https://rafael-silva-oliveira.github.io/blog/2026/segmentation-and-annotation/) - [Python] - Tutorial with a proposed VisiumHD pipeline using Proseg and Novae
- [SMURF](https://github.com/The-Mitra-Lab/SMURF) - [Python] - Soft-segmentation and manifold-unrolling framework that assigns 2 µm Visium HD bins to StarDist/Cellpose nuclei (fractionally distributing transcripts to optimize per-cell cluster similarity) and "unrolls" cells onto Cartesian coordinates for single-cell reconstruction | [Paper](https://doi.org/10.1101/2025.05.28.656357) | [PyPI](https://pypi.org/project/pysmurf/)

### Cell Annotation

- [STEM](https://github.com/WhirlFirst/STEM) - [Python] - Cell type annotation method
- [TACIT](https://github.com/huynhkl953/TACIT) - [Python] - Automated cell type identification
- [moscot](https://github.com/theislab/moscot) - [Python] - Optimal transport-based cell mapping
- [CELLama](https://github.com/portrai-io/CELLama) - [Python] - Cell annotation model
- [TACCO](https://github.com/simonwm/tacco) - [Python] - Transfer of annotations between single-cell datasets
- [Tangram](https://github.com/broadinstitute/Tangram) - [Python] - Label transfer by mapping an scRNA-seq reference onto spatial coordinates (also listed under Gene Imputation for its expression-prediction use)
- [MMoCHi](https://github.com/donnafarberlab/MMoCHi) - [Python] - Cell annotation method
- [CytoSPACE](https://github.com/digitalcytometry/cytospace) - [Python] - High-resolution alignment of single-cell and spatial transcriptomes
- [ABCT](https://github.com/ercsb-sp/ABCT) - [R] - Anchor-based Cell Typer
- [STELLAR](https://github.com/snap-stanford/stellar) - [Python] - Annotation of spatially resolved single-cell data with STELLAR
- [Vesalius](https://github.com/WonLab-CS/Vesalius) - [R] - Multi-scale and multi-context interpretable mapping of cell states across heterogeneous spatial samples
- [STALocator](https://github.com/zhanglabtools/STALocator) - [Python] - ST-Aided Locator using deep learning to localize cells from single-cell RNA-seq data onto tissue slices
- [CMAP](https://github.com/SuoLab-GZLab/CMAP) - [Python] - Cellular Mapping of Attributes with Position, maps large-scale individual cells to precise spatial locations using divide-and-conquer strategy
- [TransST](https://doi.org/10.1186/s12859-025-06099-z) - [Python] - Transfer learning framework leveraging cell-labeled information from external sources for cell-level heterogeneity inference | [GitHub](https://github.com/shuoshuoliu/TransST)
- [InSituType](https://github.com/Nanostring-Biostats/InSituType) - [R] - Cell typing for CosMx SMI spatial transcriptomics
- [HieraType](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/HieraType) - [R] - Hierarchical cell typing using RNA + protein for CosMx SMI
- [CosMx-Cell-Profiles](https://github.com/Nanostring-Biostats/CosMx-Cell-Profiles) - [R] - Collection of reference datasets for CosMx SMI
- [GARDEN](https://github.com/Briskzxm/GARDEN) - [Python] - Graph-based dynamic attention framework for identifying rare pathogenic cell populations (disease-driving cells often missed by standard methods), enables 3D tissue reconstruction
- [Spatial-ID](https://github.com/TencentAILabHealthcare/spatialID) - [Python] - Supervision-based cell typing for high-throughput cell-level SRT via transfer learning from scRNA-seq + spatial embedding | [Paper](https://doi.org/10.1038/s41467-022-35288-0)
- [HiCAT](https://github.com/jinghuang-stats/HiCAT) - [Python] - Hierarchical atlas-guided annotation transfer for cohort-scale spatial omics | [Preprint](https://doi.org/10.64898/2026.05.27.728266)

### Cell Deconvolution

- [RCTD](https://github.com/dmcable/spacexr) - [R] - Robust cell type decomposition
- [rctd-py](https://github.com/p-gueguen/rctd-py) - [Python] - Python reimplementation of the RCTD algorithm with GPU acceleration 
- [Cell2location](https://github.com/BayraktarLab/cell2location) - [Python] - Mapping scRNA-seq to spatial data
- [SPOTlight](https://github.com/MarcElosua/SPOTlight) - [R] - Seeded NMF regression to deconvolute spatial spots
- [CARD](https://github.com/YMa-lab/CARD) - [R] - Spatially informed cell-type deconvolution
- [FlashDeconv](https://github.com/cafferychen777/flashdeconv) - [Python] - Atlas-scale spatial deconvolution via structure-preserving sketching with linear O(N) scalability | [Paper](https://doi.org/10.64898/2025.12.22.696108) | [PyPI](https://pypi.org/project/flashdeconv/)
- [UCASpatial](https://github.com/BIGHanLab/UCASpatial) - [R] - Ultra-precision spatial deconvolution using entropy-based weighting of cell-identity genes to robustly map low-abundance and transcriptionally heterogeneous cell subpopulations | [Paper](https://doi.org/10.1038/s41467-026-70645-3) | [Docs](https://bighanlab.github.io/UCASpatial/)
- [RETROFIT](https://github.com/qunhualilab/retrofit) - [R] - Bayesian reference-free deconvolution requiring no single-cell reference or marker genes; effective down to Visium HD near-single-cell resolution | [Paper](https://doi.org/10.1038/s41467-026-74928-7) | [Bioconductor](https://bioconductor.org/packages/retrofit)
- [SpaCET](https://github.com/data2intelligence/SpaCET) - [R] - Reference-free tumour-aware deconvolution that infers malignant cell fractions from CNV signal | [Paper](https://doi.org/10.1038/s41467-023-36062-6)
- [spacedeconv](https://github.com/omnideconv/spacedeconv) - [R] - Unified interface to 30+ spatial deconvolution tools behind one input/output contract

### Differential Expression

- [C-SIDE](https://github.com/dmcable/spacexr) - [R] - Cell type-Specific Inference of Differential Expression in spatial transcriptomics
- [Niche-DE](https://github.com/kaishumason/NicheDE) - [R] - Niche-differential gene expression analysis identifying context-dependent cell-cell interactions
- [smiDE](https://www.biorxiv.org/content/10.1101/2024.08.02.606405v2) - [R] - Spatial differential expression method | [GitHub](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/smiDE)
- [spatialGE](https://github.com/FridleyLab/spatialGE) - [R] - Spatial gene expression analysis
- [Vespucci](https://github.com/neurorestore/Vespucci) - [R] - Prioritize spatial regions involved in the response to an experimental perturbation in spatial transcriptomics
- [CSDE](https://github.com/YosefLab/CSDE) - [Python] - Corrected Spatial Differential Expression using Prediction-Powered Inference to account for preprocessing uncertainties (segmentation, quantification, cell typing)
- [SpNeigh](https://github.com/jinming-cheng/SpNeigh) - [R] - Boundary- and gradient-aware spatial differential expression for high-resolution ST (Xenium, MERFISH, Visium HD): neighborhood extraction, distance-weighted/spline DE, and spatial enrichment scoring | [Paper](https://doi.org/10.1093/nargab/lqag039)
- [SpaceMarkers](https://github.com/DeshpandeLab/SpaceMarkers) - [R] - Infers interaction-driven molecular changes from overlapping latent-space patterns: genes that change because two cell populations are adjacent | [Paper](https://doi.org/10.1016/j.cels.2023.03.004)

### Spatially Variable Genes

- [PROST](https://github.com/Tang-Lab-super/PROST) - [Python] - Detection of spatially variable genes
- [SpatialDE](https://github.com/Teichlab/SpatialDE) - [Python] - Spatial differential expression analysis
- [SPARK-X](https://github.com/xzhoulab/SPARK) - [R] - Detection of spatially variable genes, best performing
- [Hotspot](https://github.com/YosefLab/Hotspot) - [Python] - Identify informative gene modules with lowest false positive rate
- [SOMDE](https://github.com/XuegongLab/somde) - [Python] - Self-organizing map for spatially variable gene detection with optimization
- [trendsceek](https://github.com/edsgard/trendsceek) - [R] - Identification of spatial expression trends
- [nnSVG](https://github.com/lmweber/nnSVG) - [R] - Scalable identification of spatially variable genes using nearest-neighbor Gaussian processes
- [SLOPER](https://github.com/chitra-lab/SLOPER) - [Python] - Score-based learning of Poisson-modeled expression rates for spatial gene modules and tissue organization patterns
- [FlashS](https://github.com/cafferychen777/FlashS) - [Python] - Frequency-domain Gaussian kernel testing for SVG detection, where expression sparsity accelerates rather than hinders computation | [PyPI](https://pypi.org/project/flashs/)
- [MERINGUE](https://github.com/JEFworks-Lab/MERINGUE) - [R] - Spatial autocorrelation and cross-correlation analysis robust to non-uniform cell density | [Paper](https://doi.org/10.1101/gr.271288.120)

### Integration

- [PRECAST](https://github.com/feiyoung/PRECAST) - [R] - Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data
- [MISO](https://github.com/kpcoleman/miso) - [Python] - MultI-modal Spatial Omics for versatile feature extraction and clustering (Coleman Lab)
- [SpatialMETA](https://github.com/WanluLiuLab/SpatialMETA) - [Python] - Joint analysis of spatial transcriptomics and metabolomics via CVAE | [Docs](https://spatialmeta.readthedocs.io/)
- [pyWNN](https://github.com/dylkot/pyWNN) - [Python] - Weighted Nearest Neighbors (WNN) implementation for Scanpy
- [MOFA-FLEX](https://github.com/bioFAM/mofaflex) - [Python] - Factor model framework for integrating omics data with prior knowledge | [Docs](https://mofaflex.readthedocs.io/)
- [SIMO](https://github.com/ZJUFanLab/SIMO) - [Python] - Spatial Integration of Multi-Omics through probabilistic alignment integrating spatial transcriptomics with multiple single-cell modalities
- [GSI](https://doi.org/10.1093/bioinformatics/btaf350) - [Python] - Gene Spatial Integration using deep learning with representation learning to extract spatial distribution of genes | [GitHub](https://github.com/Riandanis/Spatial_Integration_GSI) | [Zenodo](https://doi.org/10.5281/zenodo.15165223)
- [SPACE-seq](https://doi.org/10.1073/pnas.2424070122) - [Paper] - Unified molecular approach for spatial multiomics enabling simultaneous analysis of chromatin accessibility, mitochondrial DNA mutations, and gene expression on standard 10× Genomics Visium CytAssist platform
- [LLOKI](https://github.com/elliehaber07/LLOKI) - [Python] - Cross-platform spatial transcriptomics integration using optimal transport and scGPT foundation models for unified features across different gene panels (RECOMB 2025)
- [LYNX](https://github.com/azizilab/Lynx) - [Python] - Deep generative model integrating paired spatial multi-modal data (RNA, protein, metabolomics, H&E) from adjacent sections to infer microenvironmental gradients and cell-state transitions | [Docs](https://lynx-spatial.readthedocs.io)
- [INSPIRE](https://github.com/jiazhao97/INSPIRE) - [Python] - Adversarial GNN plus NMF integration of many ST datasets, yielding interpretable spatial factors and gene programs; scales to Stereo-seq | [Paper](https://doi.org/10.1038/s41588-026-02579-x)
- [SpaMosaic](https://github.com/JinmiaoChenLab/SpaMosaic) - [Python] - Mosaic integration of spatial multi-omics with only partially overlapping modalities, with missing-modality imputation | [Paper](https://doi.org/10.1038/s41588-026-02573-3)

### Cell Niches & Tissue Domains

> (**Smaller**) Cell types → Cell modules/neighborhoods → Niches/tissue domains (**Larger**)

- [BANKSY](https://github.com/prabhakarlab/Banksy) - [R/Python] - Unified cell typing and tissue domain segmentation
- [CellCharter](https://github.com/CSOgroup/cellcharter) - [Python] - Hierarchical niche detection
- [SpatialGLUE](https://github.com/JinmiaoChenLab/SpatialGlue) - [Python] - Multi-omics cell niche identification
- [smoothclust](https://github.com/lmweber/smoothclust) - [R] - Spatial clustering
- [SpaTopic](https://github.com/compbioNJU/SpaTopic) - [R] - Spatial topic modeling
- [hdWGCNA](https://smorabit.github.io/hdWGCNA/articles/ST_basics.html) - [R] - Weighted gene correlation network analysis
- [GASTON](https://github.com/raphael-group/GASTON) - [Python] - Graph-based spatial domain detection
- [PAST](https://github.com/lizhen18THU/PAST) - [Python] - Prior-based self-attention method for spatial transcriptomics tissue domain identification
- [SpatialMNN](https://github.com/Pixel-Dream/spatialMNN) - [R] - Identification of shared niches between slides
- [NicheCompass](https://github.com/Lotfollahi-lab/nichecompass) - [Python] - End-to-end analysis of spatial multi-omics data
- [Proust](https://github.com/JianingYao/proust_paper) - [Python] - Spatial domain detection using contrastive self-supervised learning for spatial multi-omics technologies (multi-modal domains)
- [STAMP](https://github.com/JinmiaoChenLab/scTM) - [Python] - Spatial Transcriptomics Analysis with topic Modeling, provides interpretable dimension reduction through deep generative modeling discovering tissue domains and cellular communication patterns
- [DeepGFT](https://doi.org/10.1186/s13059-025-03631-5) - [Python] - Combines deep learning with graph Fourier transform for spatial domain identification | [GitHub](https://github.com/jxLiu-bio/DeepGFT) | [Zenodo](https://doi.org/10.5281/zenodo.15073243)
- [SpatialFusion](https://github.com/uhlerlab/spatialfusion) - [Python] - A lightweight multimodal foundation model for pathway-informed spatial niche mapping
- [SpaHDmap](https://github.com/sldyns/SpaHDmap) - [Python] - High-definition spatial embedding integrating expression NMF with histology image encoder-decoder for spatial domain detection at enhanced resolution
- [SpatialEcoTyper](https://github.com/digitalcytometry/spatialecotyper) - [R] - Discovers and recovers spatially distinct multicellular communities from ST, scRNA-seq, and bulk data | [Paper](https://doi.org/10.1038/s41586-026-10452-4) | [Docs](https://digitalcytometry.github.io/spatialecotyper)
- [STAGATE](https://github.com/zhanglabtools/STAGATE) - [Python] - Graph-attention autoencoder for spatial domain identification; standard benchmark comparator (upstream inactive since 2022) | [Paper](https://doi.org/10.1038/s41467-022-29439-6)
- [GraphST](https://github.com/JinmiaoChenLab/GraphST) - [Python] - Graph self-supervised contrastive learning for spatial domains, vertical/horizontal integration and deconvolution | [Paper](https://doi.org/10.1038/s41467-023-36796-3)
- [lisaClust](https://github.com/SydneyBioX/lisaClust) - [R] - Clusters tissue into spatial regions using local indicators of spatial association between cell types
- [SpaNiche](https://github.com/SiyuanHuang1/SpaNiche) - [R] - Graph-regularized joint NMF over cell abundance and ligand-receptor expression to call colocalization patterns and multi-sample ecotypes | [Paper](https://doi.org/10.1186/s13059-026-04069-z)

### Cell Distances & Neighborhood

- [CRAWDAD](https://github.com/JEFworks-Lab/CRAWDAD) - [R] - Cell relationship analysis with directional adjacency distributions
- [HoodscanR](https://github.com/DavisLaboratory/hoodscanR) - [R] - Neighborhood analysis
- [SpicyR](https://sydneybiox.github.io/spicyR/) - [R] - Spatial analysis in R
- [MISTy](https://github.com/saezlab/mistyR) - [R] - Explainable multiview framework for dissecting spatial relationships from highly multiplexed data
- [SpatialCorr](https://github.com/mbernste/SpatialCorr) - [Python] - Identifying gene sets with spatially varying correlation structure
- [CatsCradle](https://github.com/AnnaLaddach/CatsCradle) - [R] - Spatial analysis framework for tissue neighbourhoods
- [RIPPLE](https://github.com/Maier-Lab/RIPPLE) - [R] - Replicate-aware detection of cell-type-anchored proximity gradients, with GPU permutation testing | [Preprint](https://doi.org/10.64898/2026.07.23.740288)

### Spatial Trajectories

- [spaTrack](https://github.com/yzf072/spaTrack) - [Python] - Spatial trajectory analysis
- [scSpace](https://github.com/ZJUFanLab/scSpace) - [Python] - Reconstruction of cell pseudo-space from single-cell RNA sequencing data
- [SOCS](https://github.com/algo-bio-lab/SOCS) - [Python] - Accurate trajectory inference in time-series spatial transcriptomics with structurally-constrained optimal transport
- [STORIES](https://github.com/cantinilab/stories) - [Python] - Spatiotemporal Reconstruction Using Optimal Transport for cell trajectory inference from spatial transcriptomics profiled at multiple time points
- [ONTrac](https://github.com/gyuanlab/ONTraC) - [Python] - Ordered Niche Trajectory Construction

### Cell-Cell Communication

- [Spatia](https://github.com/yunguan-wang/Spacia) - [Python] - Spatial cell-cell interaction analysis
- [CellAgentChat](https://github.com/mcgilldinglab/CellAgentChat) - [Python] - Agent-based cell communication modeling
- [SpaTalk](https://github.com/ZJUFanLab/SpaTalk) - [R] - Knowledge-graph-based cell-cell communication inference
- [SpaOTsc](https://github.com/zcang/SpaOTsc) - [Python] - Infers signalling relationships between cells via optimal transport (also listed under Gene Imputation for its spatial-mapping use)
- [MISTy](https://github.com/saezlab/mistyR) - [R] - Multi-view modelling of ligand-receptor and pathway views to attribute expression to intercellular signalling (also listed under Cell Distances for its neighbourhood views)
- [DeepLinc](https://github.com/xryanglab/DeepLinc) - [Python] - De novo reconstruction of cell interaction landscapes
- [CellChat](https://github.com/jinworks/CellChat?tab=readme-ov-file) - [R] - Inferrence of cell-cell communication from multiple spatially resolved transcriptomics datasets
- [COMMOT](https://github.com/zcang/COMMOT) - [Python] - Screening cell-cell communication in spatial transcriptomics via collective optimal transport
- [NicheNet](https://github.com/saeyslab/nichenetr) - [R] - Linking ligands to downstream target gene regulation
- [DeepTalk](https://github.com/JiangBioLab/DeepTalk) - [Python] - Single-cell resolution cell-cell communication using deep learning
- [CellNEST](https://github.com/schwartzlab-methods/CellNEST) - [Python] - Cell–cell relay networks using attention mechanisms on spatial transcriptomics
- [FlowSig](https://github.com/axelalmet/flowsig) - [Python] - Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics
- [CytoSignal](https://github.com/welch-lab/cytosignal) - [R] - Detects locations and dynamics of ligand-receptor signaling at single-cell resolution from spatial transcriptomic data (VeloCytoSignal captures temporal signaling velocity) | [Paper](https://doi.org/10.1038/s41588-026-02624-9)
- [LIANA+](https://github.com/scverse/liana) - [Python] - All-in-one CCC framework: 20+ ligand-receptor methods, consensus ranking, spatially-aware bivariate scores and MOFA-based multi-condition analysis | [Paper](https://doi.org/10.1038/s41556-024-01469-w)
- [SpatialDM](https://github.com/StatBiomed/SpatialDM) - [Python] - Bivariate Moran's I on ligand-receptor pairs with an analytical null, giving per-spot interaction calls | [Paper](https://doi.org/10.1038/s41467-023-39608-w)
- [FineST](https://github.com/StatBiomed/FineST) - [Python] - Fuses histology foundation-model features with ST for nuclei-resolved imputation and ligand-receptor interaction discovery | [Paper](https://doi.org/10.1038/s41467-026-70528-7)
- [MaskTalk](https://github.com/JiaPP1994/MaskTalk) - [Python] - Cell-identity-gated spatial lag model for target-aware cell-cell communication at single-cell resolution | [Preprint](https://doi.org/10.64898/2026.08.02.741948)

### Metacells & Scalability

- [SuperSpot](https://github.com/GfellerLab/SuperSpot) - [R] - Metacell analysis for spatial data
- [SEraster](https://github.com/JEFworks-Lab/SEraster) - [R] - Rasterization method for spatial data processing

### Subcellular Analysis

- [Sprawl](https://github.com/salzman-lab/SPRAWL) - [Python] - Subcellular transcript localization
- [Bento](https://github.com/YeoLab/bento-tools) - [Python] - Python toolkit for subcellular analysis of spatial transcriptomics data
- [FISHfactor](https://github.com/bioFAM/FISHFactor) - [Python] - Analysis of subcellular transcript patterns
- [InSTAnT](https://github.com/bhavaygg/InSTAnT) - [Python] - Intracellular spatial transcript analysis
- [troutpy](https://github.com/theislab/troutpy) - [Python] - Analysis of transcripts outside segmented cells in spatial transcriptomics data
- [smoppix](https://github.com/sthawinke/smoppix) - [R] - Nonparametric probabilistic-index tests for uni- and bivariate single-molecule localization patterns, avoiding segmentation and density estimation | [Paper](https://doi.org/10.1186/s13059-026-03976-5)
- [pyTrance](https://github.com/rajewsky-lab/pytrance) - [Python] - Finds co-localizing RNAs in subcellular imaging-based ST via latent embeddings | [Preprint](https://doi.org/10.64898/2026.05.07.723470)

### Copy Number Variations

- [CalicoST](https://github.com/raphael-group/CalicoST) - [Python] - CNV detection in spatial data
- [inSituCNV](https://github.com/Moldia/InSituCNV) - [Python] - Inference of Copy Number Variations in Image-Based Spatial Transcriptomics
- [SpaCNA](https://github.com/XiDsLab/SpaCNA) - [R] - Spatially aware CNA detection using morphology-matched neighbour aggregation and a hidden Markov random field; supports 3D | [Paper](https://doi.org/10.1038/s41467-026-72284-0)
- [fastCNV](https://github.com/bioinfo-MUST/fastCNV) - [R] - Fast putative CNV detection in single-cell and spatial transcriptomics data
- [SpatialInferCNV](https://github.com/aerickso/SpatialInferCNV) - [R] - Clone calling from Visium in cancer; underpins Erickson et al., Nature 608:360 (reference implementation, inactive since 2022)
- [SPICE](https://github.com/kzb193/SPICE) - [R] - Calls somatic copy-number events from spatially resolved transcriptomics | [Preprint](https://doi.org/10.64898/2026.06.30.735508)

### Isoform Analysis

- [SPLISOSM](https://github.com/JiayuSuPKU/SPLISOSM) - [Python] - Spatial isoform statistical modeling for detecting isoform-resolution patterns (alternative splicing, polyadenylation) from spatial transcriptomics data | [Paper](https://www.nature.com/articles/s41587-025-02965-6) | [Docs](https://splisosm.readthedocs.io/)
- [Sicelore-2.1](https://github.com/ucagenomix/sicelore-2.1) - [Java] - Barcode and UMI assignment plus isoform quantification for Nanopore long-read Visium; produces the isoform-level spatial counts that downstream tests consume | [Paper](https://doi.org/10.1093/nar/gkad169)

### Transcription Factors & Gene Regulatory Networks

- [STAN](https://github.com/osmanbeyoglulab/STAN) - [R] - Spatial transcription factor analysis
- [SpaGRN](https://github.com/BGI-Qingdao/SpaGRN) - [Python] - Spatially aware GRN inference using ligand-receptor-mediated spatial co-expression (bivariate Moran/Geary) on Stereo-seq and imaging ST

## Technical Enhancements

### Slide Alignment

- [PASTE](https://github.com/raphael-group/paste)/[PASTE2](https://github.com/raphael-group/paste2) - [Python] - Probabilistic alignment of spatial transcriptomics experiments
- [SPIRAL](https://github.com/guott15/SPIRAL) - [R] - Integrating and aligning spatially resolved transcriptomics data across different experiments, conditions, and technologies
- [TOAST](https://github.com/cecca46/TOAST) - [Python] - Topography Aware Optimal Transport for Alignment of Spatial Omics Data
- [STalign](https://github.com/JEFworks-Lab/STalign) - [Python] - Alignment of spatial transcriptomics data using diffeomorphic metric mapping (JEFworks Lab)
- [SPCoral](https://github.com/LiHongCSBLab/SPCoral) - [Python] - Spatial multi-modal alignment and integration
- [SPOmiAlign](https://github.com/wangyiyuyang/SPOmiAlign) - [Python] - Multi-modal spatial alignment and integration for transcriptomics and metabolomics
- [MALDI-MSI Overlay](https://github.com/M4i-Imaging-Mass-Spectrometry/MALDI-MSI---Spatial-Transcriptomics-Overlay) - [Python] - Script for co-registration of MALDI-MSI and spatial transcriptomics from One Slide Two Worlds
- [SpaMTP](https://github.com/GenomicsMachineLearning/SpaMTP) - [R] - Spatial multi-task prediction and alignment
- [SANTO](https://github.com/leihouyeung/SANTO) - [Python] - A coarse-to-fine alignment and stitching method for spatial omics
- [3d-OT](https://github.com/dbjzs/3d-OT) - [Python] - Geometry-aware heterogeneous slice alignment of spatial multi-omics via soft-correspondence optimal transport, handling nonrigid deformation | [Paper](https://doi.org/10.1038/s41592-026-03034-9)
- [VALIS](https://github.com/MathOnco/valis) - [Python] - Automated rigid and non-rigid registration of brightfield/IF whole-slide image series and 3D serial-section reconstruction; registers the images rather than the coordinates | [Paper](https://doi.org/10.1038/s41467-023-40218-9)

### Super Resolution

- [TESLA](https://github.com/jianhuupenn/TESLA) - [Python] - Super resolution for 10X Visium
- [istar](https://github.com/daviddaiweizhang/istar) - [Python] - Super resolution for Visium
- [BayesSPACE](https://github.com/edward130603/BayesSpace) - [R] - Subspot resolution
- [Spotiphy](https://github.com/jyyulab/Spotiphy) - [Python] - Super resolution tool for spatial data

### Transcripts + Histology

- [ST-Net](https://github.com/bryanhe/ST-Net) - [Python] - Integrating spatial gene expression and tumor morphology via deep learning
- [SpaceDIVA](https://github.com/hsmaan/SpatialDIVA) - [Python] - Integration of transcript data with histological images
- [HEST](https://github.com/mahmoodlab/HEST) - [Python] - Dataset for Spatial Transcriptomics and Histology Image Analysis
- [CellLENS](https://github.com/sggao/celllens/) - [Python] - Cell Local Environment Neighborhood Scan
- [DeepSpot](https://github.com/ratschlab/DeepSpot) - [Python] - Leveraging Spatial Context for Enhanced Spatial Transcriptomics Prediction from H&E Images | [Paper](https://www.medrxiv.org/content/10.1101/2025.02.09.25321567v3)
- [SpatialWhisperer](https://github.com/zinagoodlab/spatialwhisperer) - [Python] - Trimodal contrastive foundation model embedding H&E images, spatial transcriptomes, and text into a shared 512-D space for zero-shot histopathology annotation and cross-modal retrieval (ICML 2026) | [HuggingFace](https://huggingface.co/Good-Lab/spatialwhisperer)
- [DeepSpotM](https://github.com/ratschlab/DeepSpotM) - [Python] - Multimodal foundation model predicting transcriptome-wide (~19k genes) virtual ST from H&E tiles via a LoRA-adapted pathology backbone with a cross-attention gene decoder and a gene router hypernetwork that draws frozen embeddings from DNA/RNA/protein/single-cell/text FMs (Evo2, Orthrus, ProtT5, scGPT, Apertus), enabling zero-shot prediction of unseen genes | [Paper](https://doi.org/10.64898/2026.06.19.26356060) | [HuggingFace](https://huggingface.co/ratschlab/DeepSpotM)
- [DeepSpot2Cell](https://github.com/ratschlab/DeepSpot2Cell) - [Python] - Predicting virtual single-cell spatial transcriptomics from H&E images using spot-level supervision (NeurIPS 2025 Imageomics) | [Paper](https://openreview.net/forum?id=ofCkwXQKaz)
- [SpotWhisperer](https://github.com/epigen/cellxgene) - [Python] - Molecularly informed analysis of histopathology images using natural language
- [STPath](https://github.com/Graph-and-Geometric-Learning/STPath) - [Python] - A Generative Foundation Model for Integrating Spatial Transcriptomics and Whole Slide Images
- [STFlow](https://github.com/Graph-and-Geometric-Learning/STFlow) - [Python] - Scalable generation of spatial transcriptomics from histology images via whole-slide flow matching
- [AESTETIK](https://github.com/ratschlab/aestetik) - [Python] - Representation learning for multi-modal spatially resolved transcriptomics data (Bioinformatics) | [Paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag316/8692433)

## Benchmarks

- [Deconvolution benchmark](https://www.biorxiv.org/content/10.1101/2024.04.03.586404v1.full) - [Paper] - Comprehensive comparison
- [RCTD and Cell2location benchmark](https://www.biorxiv.org/content/10.1101/2023.03.22.533802v3.full.pdf) - [Paper] - Claims these are the best methods
- [Spatial clustering benchmark](https://www.nature.com/articles/s41592-024-02215-8) - [Paper] - Comparison of clustering methods
- [Spatialbench](https://github.com/latchbio/spatialbench) - [Python] - Benchmark for evaluating AI agents on spatial biology analysis tasks
- [Nature Communications review](https://www.nature.com/articles/s41467-023-37168-7#Fig2) - [Paper] - Confirms Cell2location performance
- [Open problems benchmark](https://openproblems.bio/results/spatial_decomposition/) - [Web] - Cell2location is top performer
- [Neighborhod benchmark](https://www.biorxiv.org/content/10.1101/2025.03.31.646289v2.full.pdf) - [Paper] - New COZI method top performer
- [Kaiko.ai FM benchmark EVA](https://kaiko-ai.github.io/eva/main/leaderboards/) - [Python] - WSI benchmark
- [Benchmarking of spatial transcriptomics platforms across six cancer types](https://www.biorxiv.org/content/10.1101/2024.05.21.593407v2) - [Paper] - Comprehensive platform comparison
- [PathBench](https://birkhoffkiki.github.io/PathBench/) - [Python] - Pathology benchmark
- [SPATCH Benchmark - 2025](https://www.nature.com/articles/s41467-025-64292-3) - [Paper] - Showing Xenium performs best
- [Thunder](https://github.com/MICS-Lab/thunder) - [Python] - Pathology benchmark
- [Histoboard](https://wearewaiv.github.io/histoboard/) - [Web] - Pathology leaderboards
- [Xenium_benchmarking](https://github.com/Moldia/Xenium_benchmarking) - [Python] - Independent benchmark and best-practice analysis workflows for Xenium across 25 datasets: QC, segmentation (Baysor+Cellpose), preprocessing, SVG selection, gene imputation, and domain identification (Nature Methods 2025) | [Paper](https://doi.org/10.1038/s41592-025-02617-2)
- [SACCELERATOR](https://github.com/SpatialHackathon/SACCELERATOR) - [Python] - Community benchmarking and consensus framework for spatially aware clustering across 22 methods, >170 samples and 8 platforms; argues manual anatomical labels are biased and unsuitable as ground truth | [Paper](https://doi.org/10.1038/s41592-026-03194-8)
- [spDDB](https://github.com/Zafar-Lab/spDDB) - [Python] - Benchmarks 21 deconvolution and 18 domain-detection methods across 37 datasets and 5 platforms | [Preprint](https://doi.org/10.64898/2026.05.11.724248)
- [Benchmarking-CCI](https://github.com/LitingKu/Benchmarking-CCI) - [R] - Nine spatial cell-cell interaction methods over simulations and nine datasets (Visium, Stereo-seq, Xenium); no method is optimal across resolutions, and rankings flip with the ligand-receptor database | [Paper](https://doi.org/10.1186/s13059-026-04063-5)
- [Cross-platform deconvolution benchmark](https://github.com/AGImkeller/crossplatform_deconvolution) - [Python] - Deconvolution under cross-platform technological bias; SpatialDecon and cell2location come out ahead | [Paper](https://doi.org/10.1186/s13059-026-04222-8)
- [ctSVGbench](https://github.com/FangZY-Lab/ctSVGbench) - [R] - Six cell-type-specific SVG methods over 46 real and 666 simulated datasets, with rotation and null controls; Celina best overall | [Paper](https://doi.org/10.1093/bib/bbag190)
- [BEASTsim](https://doi.org/10.1093/bib/bbag144) - [Paper] - Benchmarks four ST simulators across five datasets; SRTsim best reference-based, scDesign3 best reference-free
- [Imputation benchmark](https://doi.org/10.1093/bib/bbag027) - [Paper] - Seven imputation methods over 23 datasets and five platforms; the winner is platform-dependent (gimVI on Visium, MAGIC on Stereo-seq, Tangram on Slide-seqV2)
- [Segmentation error impact](https://doi.org/10.1038/s41588-025-02497-4) - [Paper] - Shows segmentation error confounds differential expression, neighbour influence and ligand-receptor inference, frequently dominating the result; Proseg recommended, cellAdmix corrects

### Simulators & Ground Truth

- [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3) - [R] - Unified simulator for single-cell and spatial omics with interpretable parameters; best reference-free simulator in the BEASTsim benchmark | [Paper](https://doi.org/10.1038/s41587-023-01772-1)
- [scCube](https://github.com/ZJUFanLab/scCube) - [Python] - De-novo simulation of spatial transcriptomes with controllable spatial patterns and resolutions | [Paper](https://doi.org/10.1038/s41467-024-49445-0)
- [SRTsim](https://github.com/xzhoulab/SRTsim) - [R] - Reference-based spatially resolved transcriptomics simulator; top reference-based performer in the BEASTsim benchmark

## Datasets & Foundation Models

### Datasets

- [HISSTA](https://github.com/ercsb-sp/HISSTA/tree/v1.0) - [Python/R] - Histopathology spatial transcriptomics dataset
- [STOmicsDB](https://db.cngb.org/stomics/) - [Web] - Spatial transcriptomics database
- [STHELAR](https://github.com/MICS-Lab/STHELAR) - [Python] - Multi-tissue dataset linking spatial transcriptomics (Xenium) and histology for cell-type annotation
- [DeepSpaceDB 2.0](https://deepspacedb.com/) - [Web] - Interactive web database for large-scale Xenium exploration (628 public datasets, ~1,045 samples / 129M cells) plus Visium, stored in gene- and coordinate-chunked Zarr for sub-second browser-based expression and ROI queries | [Paper](https://doi.org/10.64898/2026.01.15.699623)
- [TCGA virtual spatial transcriptomics atlas](https://huggingface.co/datasets/ratschlab/TCGA_virtual_spatial_transcriptomics_atlas) - [Python] - DeepSpot-M predicted transcriptome-wide ST for TCGA H&E (FF + FFPE); ~28.7k slides / 32 cancer types / ~296M spots (gated) | [Paper](https://www.medrxiv.org/content/10.64898/2026.06.19.26356060v1)
- [HEST Xenium virtual spatial transcriptomics](https://huggingface.co/datasets/ratschlab/HEST_Xenium_virtual_spatial_transcriptomics) - [Python] - DeepSpot-M predicted transcriptome-wide ST for 59 HEST-1k 10x Xenium H&E slides (~13.3M cells; gated) | [Paper](https://www.medrxiv.org/content/10.64898/2026.06.19.26356060v1)

### Foundation Models

#### Expression-Centric (Transcriptomics)

- [novae](https://github.com/prism-oncology/novae) - [Python] - Deep learning foundation model for spatial domain assignments and tissue organization analysis | [Paper](https://doi.org/10.1038/s41592-025-02899-6) | [Docs](https://prism-oncology.github.io/novae/)
- [scGPT-spatial](https://github.com/bowang-lab/scGPT-spatial) - [Python] - Spatial-omic foundation model pretrained on 30M spatial profiles (SpatialHuman30M) across 821 slides
- [Nicheformer](https://doi.org/10.1038/s41592-025-02814-z) - [Python] - Transformer-based foundation model pretrained on SpatialCorpus-110M containing over 110 million cells for spatial composition and label prediction | [GitHub](https://github.com/theislab/nicheformer)
- [stFormer](https://github.com/csh3/stFormer) - [Python] - Foundation model generating contextual gene representations within spatial niches with ligand–receptor aware attention
- [TERRA](https://github.com/Lotfollahi-lab/terra) - [Python] - Self-supervised foundation model producing cell- and neighborhood-level embeddings via a Graph Transformer + JEPA on masked gene tokens, reusable across downstream tasks without retraining | [Docs](https://terra-st.readthedocs.io)

#### Visual-Omics & Multimodal (H&E + ST)

- [Loki / OmiCLIP](https://github.com/GuangyuWangLab2021/Loki) - [Python] - Visual–omics foundation model and platform bridging H&E with ST | [Docs](https://guangyuwanglab2021.github.io/Loki/)
- [SpaFoundation](https://github.com/NingZhangCSUBio/SpaFoundation) - [Python] - Visual foundation model for spatial transcriptomics using histology images alone for gene expression inference and super-resolution
- [ST-Align](https://github.com/dumbgoos/ST-Align) - [Python] - Multimodal foundation model for image–gene alignment in spatial transcriptomics
- [spEMO](https://github.com/HelloWorldLTY/spEMO) - [Python] - Framework unifying embeddings from pathology foundation models and LLMs for spatial multi-omic analysis
- [FOCUS](https://github.com/LHY1007/FOCUS) - [Python] - Foundational generative model for cross-platform unified ST enhancement conditioned on H&E images, scRNA-seq references, and spatial co-expression priors (trained on 1.7M H&E-ST pairs, 10 platforms)
- [TITAN](https://github.com/mahmoodlab/TITAN) - [Python] - A multimodal whole-slide foundation model for pathology
- [SEAL](https://github.com/mahmoodlab/SEAL) - [Python] - Spatial Expression-Aligned Learning that fine-tunes pathology foundation models (UNI, CONCH) using spatial transcriptomics data for gene expression-image alignment (Mahmood Lab)
- [Phoenix](https://github.com/peng-lab/phoenix) - [Python] - Pan-cancer virtual spatial transcriptomics from routine histology via latent flow matching generative model predicting single-cell gene expression directly from H&E images, generalizing across cohorts, donors, organs, and tissues

#### Pathology & Histology

- [Virchow](https://huggingface.co/paige-ai/Virchow) - [Python] - Foundation model for computational pathology
- [UNI](https://github.com/mahmoodlab/UNI) and [UNI2](https://huggingface.co/MahmoodLab/UNI2-h) - [Python] - Universal pathology foundation models (Mahmood Lab)
- [UNI (KatherLab)](https://github.com/KatherLab/uni) - [Python] - General-purpose self-supervised pathology foundation model
- [Google Path Foundation](https://huggingface.co/google/path-foundation) - [Python] - Self-supervised embedding model for H&E histopathology
- [CONCH](https://github.com/mahmoodlab/CONCH) - [Python] - Contrastive learning for histopathology
- [GIGApath](https://github.com/prov-gigapath/prov-gigapath) - [Python] - Large-scale pathology foundation model
- [CHIEF](https://github.com/hms-dbmi/CHIEF) - [Python] - Clinical Histopathology Imaging Evaluation Foundation Model for cancer detection and prognosis
- [Phikon-v2](https://huggingface.co/owkin/phikon-v2) - [Python] - Spatial biology foundation model
- [Bioptimus H-optimus-1](https://www.bioptimus.com/news/bioptimus-launches-h-optimus-1) - [Python] - Latest biology-focused foundation model from Bioptimus
- [Atlas 2](https://github.com/cellur-m/pathology_atlas) - [Python] - Foundation models for clinical deployment
- [DeepCell dataset](https://exploredata.deepcell.com/cell-visualizations/9/versions/1) - [Web] - CNN + human features embeddings
- [NuSPIRe](https://github.com/TongjiZhanglab/NuSPIRe) - [Python] - Self-supervised nuclear-morphology model pretrained on 15.5M DAPI nuclei, for cell typing, perturbation detection and ROI/FOV selection | [Paper](https://doi.org/10.1186/s13059-026-03987-2)

#### Proteomics

- [KRONOS](https://github.com/mahmoodlab/KRONOS) - [Python] - Foundation Model for Multiplex Spatial Proteomic Images

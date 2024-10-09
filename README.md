# Spatial Transcriptomics tools

> **Smaller** | Cell types → Cell modules - cell neighborhoods → Niches/tissue domain | **Larger**

# Spatial analyses

- [Benchmark](https://www.biorxiv.org/content/10.1101/2023.03.22.533802v3.full.pdf) claims that [RCTD](https://github.com/dmcable/spacexr) and [Cell2location](https://github.com/BayraktarLab/cell2location) are the best
- [Benchmark](https://www.nature.com/articles/s41592-024-02215-8) on spatial clustering
- Other review [here](https://www.nature.com/articles/s41467-023-37168-7#Fig2), Cell2location again
- Cell2location is the best in the [open problems](https://openproblems.bio/results/spatial_decomposition/)

[A standard for sharing spatial transcriptomics data](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00171-4)

- [Detection of SVGs](https://www.nature.com/articles/s41467-024-44835-w) with PROST
- [Super resolution 10X Visium with TESLA](https://www.cell.com/cell-systems/pdf/S2405-4712(23)00084-4.pdf)
- Super resolution Visium with istar https://www.nature.com/articles/s41587-023-02019-9
- [Subspot resolution with BayesSPACE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763026/)
- https://github.com/crazyhottommy/awesome_spatial_omics

## General Tools

- Best practices

[Best Practices for Spatial Transcriptomics Analysis with Bioconductor - 2  Spatial transcriptomics](https://lmweber.org/BestPracticesST/chapters/spatial-transcriptomics.html)

- **Doing gene imputation is [not recommended](https://github.com/BayraktarLab/cell2location/issues/379)**
- scverse / [squidpy](https://squidpy.readthedocs.io/en/stable/)
- [Giotto](https://giottosuite.readthedocs.io/en/latest/) suite
- [Vitessce](http://vitessce.io/)
- [Voyager](https://github.com/pachterlab/voyager) from Pachter lab
- [scran](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03241-7) for normalization
- Multiple sample analysis with [BASS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03361-0)
- SpaVAE to do everything, reduction, visualization, clustering, batch integration, denoising, differential expression, spatial interpolation, and resolution enhancement.
- This paper recommends [cell volume normalization](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03303-w#Fig7).

### QC

- [SpaceTrooper](https://htmlpreview.github.io/?https://github.com/drighelli/SpaceTrooper/blob/main/vignette/introduction.html)

### Cell segmentation

- [Baysor](https://github.com/kharchenkolab/Baysor)
- [Cellpose](https://github.com/MouseLand/cellpose) to tune to highly specific system with few training examples
    - Able to do supersampling / restoration with [cellpose 3](https://www.biorxiv.org/content/10.1101/2024.02.10.579780v2)
- [DeepCell](https://github.com/vanvalenlab/deepcell-tf?tab=readme-ov-file)
- [Best than SOTA segmentation - Nature Methods 2024 - Bo Wang](https://www.nature.com/articles/s41592-024-02233-6#Fig3)
- [Proseg](https://www.biorxiv.org/content/10.1101/2024.04.25.591218v1.full.pdf)
- [Bin2Cell](https://www.biorxiv.org/content/10.1101/2024.06.19.599766v1) to segment Visium data
- https://bioimage.io/#/ segmentation repo
- [ComSeg](https://www.nature.com/articles/s42003-024-06480-3) uses transcripts to segment
- [FICTURE](https://www.nature.com/articles/s41592-024-02415-2)

### Cell niches

- [BANKSY](https://www.nature.com/articles/s41588-024-01664-3) **unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis**
- https://twitter.com/Cancer_dynamics/status/1770811468578971905?t=sZhz7UBO-VE7cqK97rP8sA&s=19 staryfish
- https://buff.ly/49Oc0M2 TISSUE
- Hierarchy of Cell niches with CellCharter. Looking at the most stable solutions. There is a hierarchy of niches at different layers/hierarches.
- [SpatialGLUE](https://www.nature.com/articles/s41592-024-02316-4) for cell niches **with multi-omics**
- [smoothclust](https://github.com/lmweber/smoothclust) from lukas weber
- [SpaTopic](https://www.science.org/doi/10.1126/sciadv.adp4942)

### Subcellular

- [Sprawl](https://elifesciences.org/reviewed-preprints/87517)
- [Bento](https://www.notion.so/Spatial-analyses-0d451532f4c64fc599cb6ceb469ab523?pvs=21)
- [FISHfactor](https://academic.oup.com/bioinformatics/article/39/5/btad183/7114027)
- https://www.nature.com/articles/s41467-024-49457-w
- [FICTURE](https://www.nature.com/articles/s41592-024-02415-2)

### CCC

- [Spatia](https://www.nature.com/articles/s41592-024-02408-1#Fig2)

### Cell distances

- [CRAWDAD](https://github.com/JEFworks-Lab/CRAWDAD?tab=readme-ov-file)
- [HoodscanR](https://www.biorxiv.org/content/10.1101/2024.03.26.586902v1)

### Spatial trajectories

- spaTrack

### TF/GRN

- [STAN](https://www.biorxiv.org/content/10.1101/2024.06.26.600782v2.full.pdf)

### Slide alignement

- PASTE/[PASTE2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9881963/)

### Foundation models

- https://huggingface.co/owkin/phikon-v2
- Bioptimus: https://huggingface.co/bioptimus/H-optimus-0

# Datasets

- Interesting datasets embedding CNN + human features from DeepCell https://exploredata.deepcell.com/cell-visualizations/9/versions/1

Benchmark

- https://www.biorxiv.org/content/10.1101/2024.04.03.586404v1.full

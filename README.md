# Aging-NK

Folder description：Code (Code for data analysis), Source Data (processsed and downloaded gene sets for data analysis), ExpData(Normalzied gene expression data, The raw sequencing have been deposited into GSA Human can be downloaded via accession HRA000632)

The python version used is 3.7.6, and the R version is 3.6.1

scRNA-seq analysis:
DrawFigures_inArticle.ipynb: codes used to generate data analysis graphs in the article；

FindMarkers_For_DEG.R: Calculate differential genes in the normalized expression matrix；

IntegCOVID19_Seurat3.0.R / O7Y7_Seurat3.0.R / ODYD_Seurat3.0.R: Seurat3.0 integrates single cell data from different groups for dimensionality reduction and clustering；

Monocle.PickCellByMetadata.R: Used for pseudotime and differentiation trajectory；

SCENIC_forTFs.R: Used for Transcriptional Regulons enrichment;

scRNA_calculateAUC_inDiseaseSignatiure.R: Used for Scoring disease signatures;

PAGA.ipynb: Used to generate a simpler Partition-based graph abstraction (PAGA graph).

ODYD: The scRNA-seq data of PBMCs of the young and elderly healthy individuals.

O7Y7: The scRNA-seq data of purified NK cells from PBMC of the young and elderly healthy individuals

IntegCOVID-19:  The scRNA-seq data of integrated NK.


# Aging-NK

Folder description：Code (Code for data analysis), Source Data (expression data, raw matrix and downloaded gene sets for data analysis)

The python version used is 3.7.6, and the R version is 3.6.1

scRNA-seq analysis:
DrawFigures_inArticle.ipynb: codes used to generate data analysis graphs in the article；

FindMarkers_For_DEG.R: Calculate differential genes in the normalized expression matrix；

IntegCOVID19_Seurat3.0.R / O7Y7_Seurat3.0.R / ODYD_Seurat3.0.R: Seurat3.0 integrates single cell data from different groups for dimensionality reduction and clustering；

Monocle.PickCellByMetadata.R: Used for pseudotime and differentiation trajectory；

SCENIC_forTFs.R: Used for Transcriptional Regulons enrichment;

scRNA_calculateAUC_inDiseaseSignatiure.R: Used for Scoring disease signatures;

PAGA.ipynb: Used to generate a simpler Partition-based graph abstraction (PAGA graph).




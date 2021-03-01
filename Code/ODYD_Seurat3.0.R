.libPaths("/home/huangbeibei/R/library")

Args <- commandArgs()
print(Args)

library('Seurat')
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)

varGenes <- as.integer(Args[6])
CCAdim <- as.integer(Args[7])
PCAdim <- as.integer(Args[8])
#TSNE/UMAP
DimRed <- Args[9]
#Resolution:
ResN <- Args[10]



dataDir <- '/home/huangbeibei/Aging/Seurat3.0'
PATH <- '/home/huangbeibei/Aging/Seurat3.0/AggrODYD'


outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_Res',Args[10],'_DimRed',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/IntegrateODYD')

Aggr.data <- Read10X(data.dir='/home/huangbeibei/Aging/Aggr_ODYD/outs/filtered_feature_bc_matrix/')

Aggr <- CreateSeuratObject(counts = Aggr.data, project = "DAPI", min.cells = 5, min.features = 400)
Aggr[['percent.mt']] <- PercentageFeatureSet(object=Aggr, pattern='^MT-')
Aggr <- subset(Aggr, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)


Aggr <- NormalizeData(Aggr, verbose = FALSE)
Aggr <- FindVariableFeatures(Aggr, selection.method = "vst", nfeatures = varGenes)


Aggr <- ScaleData(Aggr, verbose = FALSE)
Aggr <- RunPCA(Aggr, npcs = 40, verbose = FALSE)
print("Run PCA Done!")

pdf(file=paste0(outDir,'/ElbowPlot.PCA.pdf'))
ElbowPlot(object = Aggr, ndims = 50)
dev.off()

Aggr <- FindNeighbors(Aggr, dims = 1:PCAdim)
Aggr <- FindClusters(Aggr, resolution = as.numeric(ResN))
print("FindClusters Done!")
PCA=as.data.frame(Embeddings(object = Aggr, reduction = "pca"))
write.table(PCA,file=paste0(FilePrefix,'.pca.txt'),sep='\t',quote=F)

if (DimRed=='TSNE'){
print('TSNE')
Aggr <- RunTSNE(Aggr, dims = 1:PCAdim, reduction = "pca")
pdf(file=paste0(outDir,'/Cluster.pdf'))
DimPlot(Aggr, reduction = "tsne")
dev.off()

pdf(file=paste0(outDir,'/Batch.pdf'))
DimPlot(Aggr, reduction = "tsne", group.by = "orig.ident")
dev.off()
TSNE=as.data.frame(Embeddings(object = Aggr, reduction = "tsne"))
write.table(TSNE,file=paste0(FilePrefix,'.tsne.txt'),sep='\t',quote=F)
print("Save TSNE Done!")

}else if (DimRed=='UMAP'){
  print("UMAP")
  Aggr <- RunUMAP(Aggr, dims = 1:PCAdim, reduction = "pca")
  pdf(file=paste0(outDir,'/Cluster.pdf'))
  DimPlot(Aggr, reduction = "umap")
  dev.off()

  pdf(file=paste0(outDir,'/Batch.pdf'))
  DimPlot(Aggr, reduction = "umap", group.by = "orig.ident")
  dev.off()
  TSNE=as.data.frame(Embeddings(object = Aggr, reduction = "umap"))
  write.table(TSNE,file=paste0(FilePrefix,'.umap.txt'),sep='\t',quote=F)
  print("Save UMAP Done!")
}


#Save MetaData
Aggr@meta.data$category <- paste0(Aggr@meta.data$seurat_clusters, "_", Aggr@meta.data$orig.ident)
MetaData=as.data.frame(Aggr@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)

#Save RawDATA
RawData=as.data.frame(as.matrix(GetAssayData(object =Aggr,slot = 'counts',assay='RNA')))
write.table(RawData,file=paste0(FilePrefix,'.RawData.txt'),sep='\t',quote=F)

#SaveExpData
Data=as.data.frame(as.matrix(GetAssayData(object = Aggr)))
write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

VarGenes=VariableFeatures(object = Aggr)
VarGeneData=Data[VarGenes,]
write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.txt'),sep='\t',quote=F)


#Marker Gene
pbmc.markers=FindAllMarkers(object=Aggr,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)

MarkerGene=as.data.frame(pbmc.markers %>% group_by(cluster) %>% top_n(30,avg_logFC))
write.table(MarkerGene,file=paste(FilePrefix,'.30MarkerGene.wilcox.txt'),quote =F,sep='\t')
print("Save MarkerGene File Done!")

pdf(paste0(outDir,'/Scatter_T.pdf'),width=12, height=10)
FeaturePlot(Aggr,features = c('CCR7','CD3D','CD8A','CD79B','S100A4','CD14'))
dev.off()

pdf(paste0(outDir,'/Scatter_NK.pdf'),width=12, height=10)
FeaturePlot(Aggr,features = c('FCER1G','GZMB','KLRB1','KLRC2','PRF1'))
dev.off()

saveRDS(Aggr,paste0(outDir,'/IntegODYD.rds'))

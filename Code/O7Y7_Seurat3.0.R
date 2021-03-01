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
PATH <- '/home/huangbeibei/Aging/Seurat3.0/IntergateO7Y7'


outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_Res',Args[10],'_DimRed',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/IntegrateO7Y7')

O7.data <- Read10X(data.dir='/home/huangbeibei/Aging/Seurat3.0/O7cat_filtered_feature_bc_matrix/')
Y7.data <- Read10X(data.dir='/home/huangbeibei/Aging/Seurat3.0/Y7cat_filtered_feature_bc_matrix/')


O7 <- CreateSeuratObject(counts = O7.data, project = "O7", min.cells = 5, min.features = 400)
O7[['percent.mt']] <- PercentageFeatureSet(object=O7, pattern='^MT-')
O7 <- subset(O7, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)
Y7 <- CreateSeuratObject(counts = Y7.data, project = "Y7", min.cells = 5, min.features = 400)
Y7[['percent.mt']] <- PercentageFeatureSet(object=Y7, pattern='^MT-')
Y7 <- subset(Y7, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)

Aging.big <- merge(O7, y = c(Y7), add.cell.ids = c("O7", "Y7"), project = "Aging")

Aging_Nointeg <- NormalizeData(Aging.big, verbose = FALSE)
Aging_Nointeg <- FindVariableFeatures(Aging.big, selection.methO7 = "vst", nfeatures = varGenes,verbose = FALSE)

Aging.list <- SplitObject(Aging.big, split.by = "orig.ident")
for (i in 1:length(Aging.list)) {
  Aging.list[[i]] <- NormalizeData(Aging.list[[i]], verbose = FALSE)
  Aging.list[[i]] <- FindVariableFeatures(Aging.list[[i]], selection.methO7 = "vst", nfeatures = varGenes,
                                             verbose = FALSE)
}

reference.list <- Aging.list[c("O7", "Y7")]
Aging.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:CCAdim)
print('FindIntegrationAnchors')
Aging.integrated <- IntegrateData(anchorset = Aging.anchors, dims = 1:CCAdim)
DefaultAssay(Aging.integrated) <- "integrated"
print('IntegrateData')
#top50 <- head(VariableFeatures(Aging.integrated), 50)
#pdf(file=paste0(FilePrefix,'VarGene.pdf'),width=20, height=7)
#plot1 <- VariableFeaturePlot(Aging.integrated)
#plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))
#dev.off()

Aging.integrated <- ScaleData(Aging.integrated, verbose = FALSE)
Aging.integrated <- RunPCA(Aging.integrated, npcs = 40, verbose = FALSE)
print("Run PCA Done!")

pdf(file=paste0(outDir,'/ElbowPlot.PCA.pdf'))
ElbowPlot(object = Aging.integrated, ndims = 50)
dev.off()

Aging.integrated <- FindNeighbors(Aging.integrated, dims = 1:PCAdim)
Aging.integrated <- FindClusters(Aging.integrated, resolution = as.numeric(ResN))
print("FindClusters Done!")
PCA=as.data.frame(Embeddings(object = Aging.integrated, reduction = "pca"))
write.table(PCA,file=paste0(FilePrefix,'.pca.txt'),sep='\t',quote=F)

if (DimRed=='TSNE'){
print('TSNE')
Aging.integrated <- RunTSNE(Aging.integrated, dims = 1:PCAdim, reduction = "pca")
pdf(file=paste0(outDir,'/Cluster.pdf'))
DimPlot(Aging.integrated, reduction = "tsne")
dev.off()

pdf(file=paste0(outDir,'/Batch.pdf'))
DimPlot(Aging.integrated, reduction = "tsne", group.by = "orig.ident")
dev.off()
TSNE=as.data.frame(Embeddings(object = Aging.integrated, reduction = "tsne"))
write.table(TSNE,file=paste0(FilePrefix,'.tsne.txt'),sep='\t',quote=F)
print("Save TSNE Done!")

}else if (DimRed=='UMAP'){
  print("UMAP")
  Aging.integrated <- RunUMAP(Aging.integrated, dims = 1:PCAdim, reduction = "pca")
  pdf(file=paste0(outDir,'/Cluster.pdf'))
  DimPlot(Aging.integrated, reduction = "umap")
  dev.off()

  pdf(file=paste0(outDir,'/Batch.pdf'))
  DimPlot(Aging.integrated, reduction = "umap", group.by = "orig.ident")
  dev.off()
  TSNE=as.data.frame(Embeddings(object = Aging.integrated, reduction = "umap"))
  write.table(TSNE,file=paste0(FilePrefix,'.umap.txt'),sep='\t',quote=F)
  print("Save UMAP Done!")
}

#Save MetaData
Aging.integrated@meta.data$category <- paste0(Aging.integrated@meta.data$seurat_clusters, "_", Aging.integrated@meta.data$orig.ident)
MetaData=as.data.frame(Aging.integrated@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)

#SaveExpData
#Data=as.data.frame(as.matrix(GetAssayData(object = Aging.integrated)))
#write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

VarGenes=VariableFeatures(object = Aging.integrated)
VarGeneData=Data[VarGenes,]
write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.Integrate.txt'),sep='\t',quote=F)


###save file
saveRDS(Aging.integrated,paste0(outDir,'/IntegO7Y7.rds'))

#Aging_Nointeg
#SaveExpData
Data=as.data.frame(as.matrix(GetAssayData(object = Aging_Nointeg)))
write.table(Data,file=paste0(FilePrefix,'.dataNorm.NoIntegrate.txt'),sep='\t',quote=F)

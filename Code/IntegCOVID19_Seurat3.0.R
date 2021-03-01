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
#Resolution:
ResN <- Args[9]
#TSNE/UMAP
DimRed <- Args[10]


dataDir <- '/home/qukun/huangbb/Aging/COV19/SeuratInteg_byGroup5_O65'
PATH <- '/home/qukun/huangbb/Aging/COV19/SeuratInteg_byGroup5_O65'


outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_Res',Args[9],'_DimRed',Args[10])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/IntegNk')

O7.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/O7cat/')
Y7.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/Y7cat/')

#HC
HC_Y.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/COV19/ExtractNK/HC_Y/HC_Y/')
HC_O.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/COV19/ExtractNK/HC_O/HC_O/')

#COV
COV_Y.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/COV19/ExtractNK/COV_Y/COV_Y/')
COV_O.data <- Read10X(data.dir='/home/qukun/huangbb/Aging/COV19/ExtractNK/COV_O_year65/COV_O_year65/')

O7 <- CreateSeuratObject(counts = O7.data, project = "O7", min.cells = 5, min.features = 400)
O7[['percent.mt']] <- PercentageFeatureSet(object=O7, pattern='^MT-')
O7 <- subset(O7, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)
Y7 <- CreateSeuratObject(counts = Y7.data, project = "Y7", min.cells = 5, min.features = 400)
Y7[['percent.mt']] <- PercentageFeatureSet(object=Y7, pattern='^MT-')
Y7 <- subset(Y7, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)

#HC:
HC_Y <- CreateSeuratObject(counts = HC_Y.data, project = "HC_Y", min.cells = 5, min.features = 400)
HC_Y[['percent.mt']] <- PercentageFeatureSet(object=HC_Y, pattern='^MT-')
HC_Y <- subset(HC_Y, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)
HC_O <- CreateSeuratObject(counts = HC_O.data, project = "HC_O", min.cells = 5, min.features = 400)
HC_O[['percent.mt']] <- PercentageFeatureSet(object=HC_O, pattern='^MT-')
HC_O <- subset(HC_O, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)

#COV:
COV_Y <- CreateSeuratObject(counts = COV_Y.data, project = "COV_Y", min.cells = 5, min.features = 400)
COV_Y[['percent.mt']] <- PercentageFeatureSet(object=COV_Y, pattern='^MT-')
COV_Y <- subset(COV_Y, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)
COV_O <- CreateSeuratObject(counts = COV_O.data, project = "COV_O", min.cells = 5, min.features = 400)
COV_O[['percent.mt']] <- PercentageFeatureSet(object=COV_O, pattern='^MT-')
COV_O <- subset(COV_O, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 10)



Aging.big <- merge(O7, y = c(Y7,HC_Y,HC_O,COV_Y,COV_O), add.cell.ids = c("O7","Y7","HC_Y","HC_O","COV_Y","COV_O"), project = "Aging")

Aging.big@meta.data$Age <- c(rep("Old", nrow(O7@meta.data)),
                                 rep("Young", nrow(Y7@meta.data)+nrow(HC_Y@meta.data)),
                                 rep("Old", nrow(HC_O@meta.data)),
                                 rep("Young", nrow(COV_Y@meta.data)),
                                 rep("Old", nrow(COV_O@meta.data)))


Aging.big@meta.data$Disease <- c(rep("HC", nrow(O7@meta.data)),
                                 rep("HC", nrow(Y7@meta.data)+nrow(HC_Y@meta.data)),
                                 rep("HC", nrow(HC_O@meta.data)),
                                 rep("COV", nrow(COV_Y@meta.data)),
                                 rep("COV", nrow(COV_O@meta.data)))

Aging.big@meta.data$Group <- c(rep("HC_O", nrow(O7@meta.data)),
                                 rep("HC_Y", nrow(Y7@meta.data)+nrow(HC_Y@meta.data)),
                                 rep("HC_O", nrow(HC_O@meta.data)),
                                 rep("COV_Y", nrow(COV_Y@meta.data)),
                                 rep("COV_O", nrow(COV_O@meta.data)))

Aging.big@meta.data$Dataset <- c(rep("HC_O_1", nrow(O7@meta.data)),
                                 rep("HC_Y_1", nrow(Y7@meta.data)),
                                 rep("HC_Y_2", nrow(HC_Y@meta.data)),
                                 rep("HC_O_2", nrow(HC_O@meta.data)),
                                 rep("COV_Y", nrow(COV_Y@meta.data)),
                                 rep("COV_O", nrow(COV_O@meta.data)))


print (head(Aging.big@meta.data))

#Aging_Nointeg <- NormalizeData(Aging.big, verbose = FALSE)
#Aging_Nointeg <- FindVariableFeatures(Aging.big, selection.methO7 = "vst", nfeatures = varGenes,verbose = FALSE)

Aging.list <- SplitObject(Aging.big, split.by = "Dataset")
print (Aging.list)
for (i in 1:length(Aging.list)) {
  Aging.list[[i]] <- NormalizeData(Aging.list[[i]], verbose = FALSE)
  Aging.list[[i]] <- FindVariableFeatures(Aging.list[[i]], selection.methO7 = "vst", nfeatures = varGenes,
                                             verbose = FALSE)
}


Aging.anchors <- FindIntegrationAnchors(object.list = Aging.list, dims = 1:CCAdim)
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
Aging.integrated@meta.data$category <- paste0(Aging.integrated@meta.data$seurat_clusters, "_", Aging.integrated@meta.data$Group)
MetaData=as.data.frame(Aging.integrated@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)

#Aging.integrated@meta.data$category <- paste0(Aging.integrated@meta.data$seurat_clusters, "_", Aging.integrated@meta.data$orig.ident)
MetaData2=as.data.frame(Aging.big@meta.data)
write.table(MetaData2,file=paste0(FilePrefix,'.MetaData2.txt'),sep='\t',quote=F)

#Aging.big@meta.data

#SaveExpData
#Data=as.data.frame(as.matrix(GetAssayData(object = Aging.integrated)))
#write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

#VarGenes=VariableFeatures(object = Aging.integrated)
#VarGeneData=Data[VarGenes,]
#write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.Integrate.txt'),sep='\t',quote=F)




print("UMAP Dimplot:")
pdf(file=paste0(outDir,'/Dimplot_',DimRed,'.cluster.pdf'))
options(repr.plot.height = 4, repr.plot.width = 4)
DimPlot(Aging.integrated, reduction = "umap")
dev.off()

#saveRDS(Aging.integrated,paste0(outDir,'/IntegNK.rds'))



#Marker Gene
pbmc.markers=FindAllMarkers(object=Aging.integrated,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
AllMarkerGene=as.data.frame(pbmc.markers %>% group_by(cluster))
write.table(AllMarkerGene,file=paste(FilePrefix,'.AllMarkerGene.wilcox.Inregrate.txt'),quote =F,sep='\t')
MarkerGene=as.data.frame(pbmc.markers %>% group_by(cluster) %>% top_n(30,avg_logFC))
write.table(MarkerGene,file=paste(FilePrefix,'.30MarkerGene.wilcox.Inregrate.txt'),quote =F,sep='\t')
print("Save MarkerGene File Done!")


RawData=as.data.frame(as.matrix(GetAssayData(object = Aging.integrated,slot = 'counts')))
write.table(RawData,file=paste0(FilePrefix,'.RawData.txt'),sep='\t',quote=F)


#print("Save MarkerGene File Done!")

#Idents(Aging_Nointeg) <- "category"
#C0 <- FindMarkers(Aging_Nointeg, ident.1 = "0_O7", ident.2 = "0_Y7", verbose = FALSE,test.use = "wilcox")
#C1 <- FindMarkers(Aging_Nointeg, ident.1 = "1_O7", ident.2 = "1_Y7", verbose = FALSE,test.use = "wilcox")

#write.table(C0, file="/home/huangbeibei/Aging/Seurat3.0/IntergateO7Y7/varGenes2200_CCA25_PCA25_Res0.5_DimRedUMAP/DEgene_OvsY_NoInteg/C0_OvsY.enriched.txt",quote =F,sep='\t')
#write.table(C1, file="/home/huangbeibei/Aging/Seurat3.0/IntergateO7Y7/varGenes2200_CCA25_PCA25_Res0.5_DimRedUMAP/DEgene_OvsY_NoInteg/C1_OvsY.enriched.txt",quote =F,sep='\t')

###save file
saveRDS(Aging.integrated,paste0(outDir,'/IntegNK.rds'))

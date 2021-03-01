library(Seurat)
library(dplyr)
library(pheatmap)


Meta=read.table('/home/huangbeibei/Aging/DEgene/O7Y7/IntegrateO7Y7.MetaData.processed_forSingler.txt',sep='\t',header=TRUE,row.names = 1)


NK <- CreateSeuratObject(counts = read.table(file = '/home/huangbeibei/Aging/DEgene/O7Y7/IntegrateO7Y7_ScanpyNormLn_FilterCells_Singler.txt', row.names = 1,sep='\t',header=TRUE),meta.data=Meta,)

print ('NK object read done!')
print (NK)

print (head(NK@meta.data))

#Find Marker genes
print ('Find Marker genes...')
Idents(NK) <- "CellType"
NK.markers=FindAllMarkers(object=NK,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
MarkerGene15=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(15,avg_logFC))
MarkerGene30=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(30,avg_logFC))
MarkerGene50=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(50,avg_logFC))
write.table(MarkerGene15,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType1_MarkerGenes_Top15.txt'),quote =F,sep='\t')
write.table(MarkerGene30,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType1_MarkerGenes_Top30.txt'),quote =F,sep='\t')
write.table(MarkerGene30,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType1_MarkerGenes_Top50.txt'),quote =F,sep='\t')
print ('Find Marker done！')

#Find Marker genes
print ('Find Marker genes...')
Idents(NK) <- "CellType2"
NK.markers=FindAllMarkers(object=NK,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
MarkerGene15=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(15,avg_logFC))
MarkerGene30=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(30,avg_logFC))
MarkerGene50=as.data.frame(NK.markers %>% group_by(cluster) %>% top_n(50,avg_logFC))
write.table(MarkerGene15,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType2_MarkerGenes_Top15.txt'),quote =F,sep='\t')
write.table(MarkerGene30,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType2_MarkerGenes_Top30.txt'),quote =F,sep='\t')
write.table(MarkerGene30,file=('/home/huangbeibei/Aging/DEgene/O7Y7/CellType2_MarkerGenes_Top50.txt'),quote =F,sep='\t')
print ('Find Marker done！')

#Find Aging DEgene in CellType
print ('Find Aging DEgene genes...')
Idents(NK) <- "CellType2Group"
AdapativeNK_1 <- FindMarkers(NK, ident.1 = "O7_AdapativeNK_1", ident.2 = "Y7_AdapativeNK_1", verbose = FALSE)
AdapativeNK_2 <- FindMarkers(NK, ident.1 = "O7_AdapativeNK_2", ident.2 = "Y7_AdapativeNK_2", verbose = FALSE)
AdapativeNK_3 <- FindMarkers(NK, ident.1 = "O7_AdapativeNK_3", ident.2 = "Y7_AdapativeNK_3", verbose = FALSE)
AdapativeNK_4 <- FindMarkers(NK, ident.1 = "O7_AdapativeNK_4", ident.2 = "Y7_AdapativeNK_4", verbose = FALSE)
MatureNK <- FindMarkers(NK, ident.1 = "O7_MatureNK", ident.2 = "Y7_MatureNK", verbose = FALSE)
Nkbright <- FindMarkers(NK, ident.1 = "O7_Nkbright", ident.2 = "Y7_Nkbright", verbose = FALSE)
NKpro <- FindMarkers(NK, ident.1 = "O7_NKpro", ident.2 = "Y7_NKpro", verbose = FALSE)

write.table(AdapativeNK_1, file=("/home/huangbeibei/Aging/DEgene/O7Y7/AdapativeNK1_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(AdapativeNK_2, file=("/home/huangbeibei/Aging/DEgene/O7Y7/AdapativeNK2_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(AdapativeNK_3, file=("/home/huangbeibei/Aging/DEgene/O7Y7/AdapativeNK3_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(AdapativeNK_4, file=("/home/huangbeibei/Aging/DEgene/O7Y7/AdapativeNK4_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(MatureNK, file=("/home/huangbeibei/Aging/DEgene/O7Y7/MatureNK_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(Nkbright, file=("/home/huangbeibei/Aging/DEgene/O7Y7/Nkbright_OvsY.enriched.txt"),quote =F,sep='\t')
write.table(NKpro, file=("/home/huangbeibei/Aging/DEgene/O7Y7/NKpro_OvsY.enriched.txt"),quote =F,sep='\t')

print ('Find Aging DEgene done!')

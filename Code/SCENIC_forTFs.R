#cd /home/huangbeibei/Aging/SCENIC/ODYD

.libPaths("/home/huangbeibei/R/library")

library(SCENIC)
library(SingleCellExperiment)
library(RcisTarget)
library(rbokeh)




scenicOptions <- readRDS("./int/scenicOptions.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

cellInfo <- readRDS("int/cellInfo.Rds")
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

#CellType-Regulator
pdf(file=('output/pheatmap_regulonActivity_byCellType_Scaled.pdf'))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, color=colorRampPalette(c("#6e9ecf","white","#de593a"))(100), breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

write.table(regulonActivity_byCellType_Scaled,file=('output/pheatmap_regulonActivity_byCellType_Scaled.txt'),sep='\t',quote=F)

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))

#tSNE的坐标
#tSNE_scenic$Y

write.table(tSNE_scenic$Y,file=('output/tSNE_scenic.txt'),sep='\t',quote=F)
#NKbright_1':'Y7_NKbright_1', 'Y7_NKbright_2':'Y7_NKbright_2', 'Y7_NKpro':'Y7_NKpro'

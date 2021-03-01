
Args <- commandArgs()
print(Args)


library(AUCell)
library(GSEABase)
library(NMF)


ExpF=Args[6]
GenesetF=Args[7]
OutDir=Args[8]


#1. Load Expression matrix
exprMatrix <- read.table(ExpF)
exprMatrix <- as.matrix(exprMatrix)
print (dim(exprMatrix))


#2. Build gene-expression rankings for each cells
#For each cell, the genes are ranked from highest to lowest value
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file=paste0(OutDir,"cells_rankings.RData"))

#Name <- c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
#Name <- c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16')
#Name <- c('0','1','2','3','4','5','6','7','8')

#Name <- c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')

Name <- c('GSE10325_DEseq_SLEvsHCup_genelist_112','T2diabetes_PB_PMID25956355_Disease_up_OverlapAllGene','Arthritis_PB_GSE93272_GEO2R_Disease_up_OverlapAllGene','Atherosclerosis_PB_PMID24975946_TableS2_Disease_up_OverlapAllGene','COVup_AllCell_as_bulk','COVup_eachClusterMerge','KEGG_DILATED_CARDIOMYOPATHY','KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS','KEGG_TYPE_II_DIABETES_MELLITUS','KEGG_ALZHEIMERS_DISEASE')

for(x in Name){
    #2. Load Gene set
    #读入geneset文件
    GFr <- read.table(GenesetF)
    GF <- GFr[GFr$cluster == x,]
    rownames(GF)<-GF[,'gene']
    genes <- rownames(GF)

    geneSets <- GeneSet(genes, setName=x)
    #geneSets
    geneSets <- GeneSet(genes, setName=x)

    #cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=ceiling(0.15 * nrow(cells_rankings)))

    save(cells_AUC, file=paste0(OutDir,x,'_cellsAUC.RData'))

    CellValueDF1 <- getAUC(cells_AUC)
    CellValueDF1_T <- t(CellValueDF1)

    oldCol=rownames(CellValueDF1_T)
    newCol=gsub("\\.",'-',oldCol)
    rownames(CellValueDF1_T)=newCol
    write.table(CellValueDF1_T,file=paste0(OutDir,x,'_cellsAUC_value.txt'),sep='\t',quote=FALSE)

    pdf(file=paste0(OutDir,x,'_cellsAUC_RawCurve.pdf'))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
    dev.off()

    threshold <- cells_assignment[[x]]$aucThr$thresholds
    DFthreshold = as.data.frame(t(threshold))
    write.table(DFthreshold,file=paste0(OutDir,x,'_cellsAUC_threshold.txt'),sep='\t',quote=FALSE)
}

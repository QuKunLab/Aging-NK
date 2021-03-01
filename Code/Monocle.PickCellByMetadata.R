.libPaths("/home/huangbeibei/R/library")
Args <- commandArgs()
print(Args)

library(monocle)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)

#Read parameters
indir=Args[6]
metaData=Args[7]
num_dim=as.integer(Args[8])
varGene=as.integer(Args[9])
OutDir0=Args[10]

Path=paste0("MonocleOutput","_PC",Args[8],"_","VarGene",Args[9])
print(Path)
OutDir=paste0(OutDir0,"/",Path)
print(OutDir)
if (!dir.exists(OutDir)){
    dir.create(OutDir)}
FilePrefix=paste0(OutDir,'/Monocle')
print('OutDir')

pbmc.data=Read10X(data.dir=indir)
dense.size=object.size(x=as.matrix(x=pbmc.data))
sparse.size=object.size(x=pbmc.data)
pbmc=CreateSeuratObject(counts = pbmc.data,min.cells=10,min.features=400,project='Monocle')
mito.genes=grep(pattern='^MT-',x=rownames(x=GetAssayData(object = pbmc)),value=TRUE)
percent.mito=Matrix::colSums(GetAssayData(object = pbmc,slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(object = pbmc,slot = "counts"))
pbmc=AddMetaData(object=pbmc,metadata=percent.mito,col.name='percent.mt')
print('read10X done')

####Get cells you need
MetaData=read.table(metaData,row.names=1)
print("Filter cells")
FinalUseCells=rownames(MetaData)
####Creat Seurat object
NewPbmc=subset(pbmc,subset=nFeature_RNA>400 & nFeature_RNA < 3000 & percent.mt < 10,cells=FinalUseCells)
print("Subset Pbmc Down!")

exprs=as.matrix(GetAssayData(object =NewPbmc))
fd=as.data.frame(row.names(exprs))
colnames(fd)[1] <- "gene_short_name"
rownames(fd)<-row.names(exprs)
HSMM=newCellDataSet(exprs,featureData = new("AnnotatedDataFrame", data = fd),expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.02) #
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
print(head(pData(HSMM)))

pData(HSMM)$Total_mRNAs <- Matrix::colSums(HSMM@assayData$exprs)
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > 2000 &pData(HSMM)$Total_mRNAs < 10000]
HSMM <- detectGenes(HSMM, min_expr =0.02)

######Select DiffGene
NewPbmc=NormalizeData(object=NewPbmc,normalization.method = 'LogNormalize',scale.factor=10000)
NewPbmc=FindVariableFeatures(NewPbmc,nfeatures = varGene,mean.cutoff = c(0.05, 6),dispersion.cutoff = c(0.5, Inf))


YGene=c('DDX3Y','TTTY15','EIF1AY','RPS4Y1','XGY2','SRY','RPS4Y1','ZFY','TGIF2LY','PCDH11Y','TTTY23B','TTTY23','TSPY2','FAM197Y9','LINC00280','TTTY1B','TTTY1','TTTY2B','TTTY2','TTTY21B','TTTY21','TTTY7B','TTTY7','TTTY8B','TTTY8','AMELY','TBL1Y','TTTY16','TTTY12','LINC00279','TTTY18','TTTY19','TTTY11','RBMY1A3P','TTTY20','TSPY10','TSPY3','TSPY4','TSPY8','FAM197Y2','FAM197Y4','FAM197Y7','FAM197Y8','FAM197Y6','TSPY1','FAM197Y3','RBMY3AP','TTTY22','GYG2P1','TTTY15','USP9Y','UTY','TMSB4Y','VCY1B','VCY','NLGN4Y','NLGN4Y-AS1','FAM41AY1','FAM41AY2','FAM224B','FAM224A','XKRY2','XKRY','CDY2B','CDY2A','HSFY1','HSFY2','TTTY9A','TTTY9B','TTTY14','CD24P4','LOC107987345','TXLNGY','KDM5D','TTTY10','EIF1AY','RPS4Y2','ERVH-6','PRORY','RBMY2EP','RBMY1A1','RBMY1B','RBMY1D','RBMY1E','TTTY13','PRY2','PRY','LOC101929148','TTTY6','TTTY6B','RBMY1F','RBMY1J','TTTY5','RBMY2FP','LOC100652931','TTTY17A','TTTY17B','TTTY17C','TTTY4B','TTTY4C','TTTY4','BPY2B','BPY2C','BPY2','DAZ1','DAZ2','DAZ3','DAZ4','TTTY3B','TTTY3','CDY1','CDY1B','CSPG4P1Y','GOLGA2P2Y','GOLGA2P3Y','PRYP4')
ordering_genes=VariableFeatures(object = NewPbmc)
print('Length of ordering_genes:')
length(ordering_genes)
print('Remove ChrY gene!')
for (g in VariableFeatures(object = NewPbmc)){
    if (g %in% YGene){
        print(g)
        ordering_genes=ordering_genes[ordering_genes!=g]
    }
}
VariableFeatures(object = NewPbmc)=ordering_genes
print('Length of final VarGene:')
print(length(VariableFeatures(object = pbmc)))
print("FindVariableGenes Done!")



HSMM=setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM<- reduceDimension(HSMM,num_dim =num_dim,max_components = 2,method = 'DDRTree') #降维
HSMM <- orderCells(HSMM)


pdf(file=paste0(FilePrefix,".Pseudotime.pdf"))
print("Save Pseudotime.pdf")
plot_cell_trajectory(HSMM,color_by="Pseudotime",cell_size=0.3)
dev.off()

pdf(file=paste0(FilePrefix,".State.pdf"))
print("Save State.pdf")
plot_cell_trajectory(HSMM,cell_size=0.3)
dev.off()

PD=pData(HSMM)
write.table(PD,file=paste0(FilePrefix,".pData.txt"),quote =F,sep='\t')
print("Write table pData.txt Done!")

PseudoTime=HSMM@reducedDimS
write.table(PseudoTime,file=paste0(FilePrefix,".reducedDims.txt"),quote =F,sep='\t')
print("Write table reducedDims.txt Done!")

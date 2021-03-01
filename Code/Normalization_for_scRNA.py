
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import sys
import math
import scipy.stats
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import scipy.stats as stats

import matplotlib
colorslist = ['#002299','white','#CC3300']
cmaps = matplotlib.colors.LinearSegmentedColormap.from_list('mylist',colorslist,N=800)

import matplotlib
#colorslist = ['#6e9ecf','white','#de593a']
###fffcf0
#colorslist = ['#6a60a9','#dedcee','#fbd14b']
#colorslist = ['#6a60a9','white','#fbd14b']
colorslist = ['#6e9ecf','white','#de593a']
cc = matplotlib.colors.LinearSegmentedColormap.from_list('mylist',colorslist,N=800)


#Step0: 根据meta的细胞来调整表达矩阵的细胞
MetaF='../IntegNk.MetaData.Singler.info.txt'
MetaDF=pd.read_csv(MetaF,sep='\t',index_col=0)
CellFilter=list(Meta.index)

# Step1: Pre-process matrix(将seurat原始count矩阵行基因列细胞转换)
DataF='../IntegNk.RawData.txt'
Data=pd.read_csv(DataF,sep='\t',index_col=0)
DataT=Data.T
#DataT=DataT.loc[CellFilter]
DataT.to_csv('IntegNk_RawDataT.txt',sep='\t')

# Step2:读入列基因行细胞的原始矩阵，adata，并标准化，取对数log
#gene expression
#adata: The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.

adata = sc.read('IntegNk_RawDataT.txt', cache=True)
#Normlize
sc.pp.normalize_total(adata)

adata.T.to_df().to_csv('IntegNk_ScanpyNorm.txt',sep='\t')
print (adata.X)

#NormData=adata.T.to_df().copy()
#NormData=NormData[CellFilter]
#print (NormData.shape[1])
#NormData.to_csv('IntegNk_ScanpyNorm_FilterCellsSingler.txt',sep='\t')

#Step3: log(自然对数，𝑋=log(𝑋+1)， 输出文件
sc.pp.log1p(adata)
adata.T.to_df().to_csv('IntegNk_ScanpyNorm_Ln.txt',sep='\t')

print (adata.X)

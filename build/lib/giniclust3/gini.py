#!/usr/bin/env python
##################################################
# File Name: gini.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Mon 08 Jul 2019 03:40:50 PM EDT
################################################

import numpy as np
import scanpy as sc
import anndata
import math
import umap
from sklearn import preprocessing
from scipy.interpolate import interp1d
import statsmodels.api as sm
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def giniIndexCalculation(array):
    plusArray=array+0.0000001
    plusArray=np.sort(plusArray)
    index=np.arange(1,plusArray.shape[0]+1)
    n=plusArray.shape[0]
    indexValue=(np.sum((2*index-n-1)*plusArray))/(n*np.sum(plusArray))
    return (indexValue)

def giniIndex(filteredArray):
    funcGeneGini=[]
    funcGeneMax=[]
    for i in range(0,len(filteredArray)):
        exp=filteredArray[i]
        maxExp=np.max(exp)
        giniIndexValue=giniIndexCalculation(exp)
        funcGeneGini.append(giniIndexValue)
        funcGeneMax.append(float(maxExp))
    return(funcGeneGini,funcGeneMax)

def arctanTransform(matrix):
    transMatrix=matrix.T
    binCellNum=int(len(matrix)/1000)
    expMList=[]
    arctanExp=[]
    for i in range(0,len(transMatrix)):
        expSort=np.sort(transMatrix[i])[::-1]
        expSum=np.sum(expSort)
        expM=[]
        if (binCellNum<=9):
            for j in range(0,len(expSort)):
                binSum=np.sum(expSort[0:j])
                div=binSum/expSum
                if div>0.8:
                    expM=np.mean(expSort[0:j])
                    break
        else:
            loopNum=int((len(matrix)-binCellNum)/binCellNum)
            for j in range(loopNum):
                end=(j+1)*binCellNum
                binSum=np.sum(expSort[0:end])
                div=binSum/expSum
                if div>0.8:
                    expM=np.mean(expSort[0:end])
                    break
        arctan=np.arctan((transMatrix[i]-expM))*10+np.arctan(expM)*10
        arctanExp.append(arctan)
    arctanExp=np.array(arctanExp)
    return(arctanExp.T)

def loessRegression(funcGeneGini,funcLogGeneMax,funcGene,funcGiniCutoff,funcPvalue):
    dictGini={}
    for i in range(len(funcGene)):
        dictGini[funcGene[i]]=funcGeneGini[i]
    fit=sm.nonparametric.lowess(funcGeneGini,funcLogGeneMax,frac=0.9)
    f=interp1d(list(zip(*fit))[0],list(zip(*fit))[1], bounds_error=False)
    giniFit=f(funcLogGeneMax)
    residue=funcGeneGini-giniFit
    posRes=[]
    for i in range(len(residue)):
        posRes.append(residue[i])
    residueSort=sorted(posRes)
    quarter=int(len(residueSort)*3/4)
    cutoff=residueSort[quarter]
    quantileGini=[]
    quantileLogMax=[]
    quantileGene=[]
    outlierGini=[]
    outlierLogMax=[]
    outlierGene=[]
    for i in range(0,len(residue)):
        if (residue[i]<=cutoff):
            quantileGene.append(funcGene[i])
            quantileGini.append(funcGeneGini[i])
            quantileLogMax.append(funcLogGeneMax[i])
        else:
            outlierGini.append(funcGeneGini[i])
            outlierGene.append(funcGene[i])
            outlierLogMax.append(funcLogGeneMax[i])
    fit2=sm.nonparametric.lowess(np.array(quantileGini),np.array(quantileLogMax),frac=0.9)
    reFit=interp1d(list(zip(*fit2))[0],list(zip(*fit2))[1], bounds_error=False)
    quantileGiniReFitPredict = reFit(np.array(quantileLogMax))
    dictFunction={}
    for i in range(len(quantileLogMax)):
        dictFunction[quantileLogMax[i]]=quantileGiniReFitPredict[i]
    uniqQuantileLogMax=list(dictFunction.keys())
    sortQuantileLogMax=sorted(uniqQuantileLogMax,reverse=True)
    k=(dictFunction[sortQuantileLogMax[1]]-dictFunction[sortQuantileLogMax[2]])/(sortQuantileLogMax[1]-sortQuantileLogMax[2])
    b=dictFunction[sortQuantileLogMax[1]]-k*sortQuantileLogMax[1]
    outlierGiniReFitPredict = reFit(outlierLogMax)
    ###remove value NaN###
    for i in range(0,len(outlierGiniReFitPredict)):
        if math.isnan(outlierGiniReFitPredict[i]):
            outlierGiniReFitPredict[i]=outlierLogMax[i]*k+b
    reFitResidueQuantile=quantileGini-quantileGiniReFitPredict
    reFitResidueOutlier=outlierGini-outlierGiniReFitPredict
    reFitResidue=np.array(list(reFitResidueQuantile)+list(reFitResidueOutlier))
    newGene=quantileGene+outlierGene
    newGini=np.array(list(quantileGini)+list(outlierGini))
    newFitGini=np.array(list(quantileGiniReFitPredict)+list(outlierGiniReFitPredict))
    newLogMax=quantileLogMax+outlierLogMax
    newResidualGini=newGini-newFitGini
    pvalue=stats.norm.cdf(-abs(preprocessing.scale(newResidualGini)))
    funcSigGiniGeneGini={}
    funcSigGiniGenePvalue={}
    for i in range(0,len(pvalue)):
        if (float(pvalue[i])<funcPvalue and dictGini[newGene[i]]>=funcGiniCutoff):
            funcSigGiniGeneGini[newGene[i]]=dictGini[newGene[i]]
            funcSigGiniGenePvalue[newGene[i]]=pvalue[i]
    return(funcSigGiniGeneGini,funcSigGiniGenePvalue)

def giniValueSelectionM(funcGeneGini,funcGene,funcGiniCutoff):
    funcSigGiniGeneGini={}
    for i in range(0,len(funcGene)):
        if (float(funcGeneGini[i])>=funcGiniCutoff):
            funcSigGiniGeneGini[funcGene[i]]=funcGeneGini[i]
    return(funcSigGiniGeneGini)

def calGini(adataSC,**kwargs):
    """
    Calculate Gini Index value for each gene.
    Params
    ------
    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    selection: string, optional (Default: 'p_value')
        selection='p_value' or selection='gini_value'. 'p_value' mode indicates
        Gini genes are selected by p value of Loess regression. 'gini_value' mode
        indicates Gini genes are selected by Gini Index value. Recommend 'p_value' 
        mode.
    p_value: float, optional (Default=0.0001)
        If `selection='p_value'`, assign a p value cutoff for Gini gene selection.
    min_gini_value: float, optional (Default=0.6)
        Assign min Gini Index value cutoff in Gini gene selection. The Gini
        Index value is range from 0 to 1. Thus the min_gini_value cutoff is 
        better in range from 0 to 0.8. Larger value indicates a more stringent
        in select Gini genes, smaller values indicates more Gini genes are 
        selected

    A. Gini gene selection in 'p_value' mode: 1. p value < p_value 2. Gini Index
    value >= min_gini_value. (Recommended)
    B. In 'gini_value' mode: Gini Index value > min_gini_value.

    Returns
    -------
    Returns dictionary with Gini genes. adata.var['gini']
    """

    ginicutoff=kwargs.get('min_gini_value', 0.6)
    p_value=kwargs.get('p_value', 0.0001)
    selection=kwargs.get('selection', 'p_value')
    geneLabel=adataSC.var.index.tolist()
    cellLabel=adataSC.obs.index.tolist()
    print("Gene number is "+str(len(geneLabel)))
    print("Cell number is "+str(len(cellLabel)))

    ###calculate gini index###
    geneGini,geneMax=giniIndex(adataSC.X.T)
    allGini=np.array(list(zip(geneLabel,geneGini)))
    logGeneMax=list(map(lambda x:math.log(x+0.1,2),geneMax))
    geneGini=np.array(geneGini)
    geneMax=np.array(geneMax)
    logGeneMax=np.array(logGeneMax)
    sigGiniGene={}
    if (selection=='p_value'):
        sigGiniGene,sigGiniPvalue=loessRegression(geneGini,logGeneMax,geneLabel,ginicutoff,p_value)
    elif (selection=='gini_value'):
         sigGiniGene=giniValueSelectionM(geneGini,geneLabel,ginicutoff)
    else:
        raise SystemExit("Only p_value or gini_value mode is allowed in this step.")
    if (len(sigGiniGene)<2):
        raise SystemExit("Only "+str(len(sigGiniGene))+" Gini gene passed the cutoff, set a lower value of --min_gini_value.")

    ###Select and save high Gini genes to Anndata###
    geneGiniBool=[]
    for i in range(len(geneLabel)):
        value=0
        if geneLabel[i] in sigGiniGene.keys():
            value=1
        geneGiniBool.append(value)
    adataSC.var['gini']=np.array(geneGiniBool,dtype=bool)

def clusterGini(adataSC,**kwargs):
    """
    Cluster cell based on Gini Index value.
    Params
    ------
    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    neighbors: int, optional (Default=5)
        The size of local neighborhood used for manifold approximation. Larger
        values result in more global views of the manifold, while smaller values
        result in more local data being preserved. For rare cell identification
        this values should be in the range 2 to 15. Recommended neighbors = 5.
    resolution: float, optional (Default=0.1)
        A parameter value controlling the coarseness of the clustering. Higher 
        values lead to more clusters.
    method: string, optional (Default: 'leiden')
        method='louvain' or method='leiden'.

    Returns
    -------
    Returns dictionary with gini cluster result. adata.var['rare']
    """
    
    cluster_neighbors=kwargs.get('neighbors', 5)
    cluster_resolution=kwargs.get('resolution', 0.1)
    cluster_method=kwargs.get('method', "leiden")
    if (cluster_method!="louvain" and cluster_method!="leiden"):
        raise SystemExit("Only leiden or louvain cluster method is allowed in this step.")
    adataGini=adataSC[:,adataSC.var['gini']]
    scaleMatrix=arctanTransform((adataGini.X))
    adataScaleGini=anndata.AnnData(X=scaleMatrix)

    ###calculate neighbor and clustering###
    sc.pp.neighbors(adataScaleGini,use_rep='X',n_neighbors=cluster_neighbors)
    giniClust=[]
    if (cluster_method=="louvain"):
        sc.tl.louvain(adataScaleGini,resolution=cluster_resolution)
        giniClust=adataScaleGini.obs['louvain'].values.tolist()
    else:
        sc.tl.leiden(adataScaleGini,resolution=cluster_resolution)
        giniClust=adataScaleGini.obs['leiden'].values.tolist()
    adataSC.obs['rare']=giniClust
    return(adataScaleGini)

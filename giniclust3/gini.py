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

def rareClusterId(funcClust,markGene):
    hashCluster={}
    for element in funcClust:
        hashCluster[element]=hashCluster.get(element,0)+1
    clusterID=hashCluster.keys()
    clusterID=sorted(map(int,clusterID))
    CandidateGene=[]
    CandidateClust=[]
    newClustID=0
    hashProjectID={}
    for i in range(len(clusterID)):
        countCandidateGene=0
        subCandidateGene=[]
        #####remove large clusters with component >= 5%
        if (hashCluster[str(clusterID[i])]/len(funcClust)>=0.01):
            continue
        for j in range(len(markGene['names'])):
            #####count marker genes with more than 5 fold enriched in select cluster####
            if (markGene['logfoldchanges'][j][i]>=np.log2(5) and markGene['pvals_adj'][j][i]<0.01):
                countCandidateGene+=1
                subCandidateGene.append(markGene['names'][j][i])
        if (countCandidateGene>=3):
            newClustID+=1
            CandidateGene.extend(subCandidateGene)
            hashProjectID[str(clusterID[i])]=str(newClustID)
    newCluster=[]
    for i in range(len(funcClust)):
        if (funcClust[i] in hashProjectID.keys()):
            newCluster.append(hashProjectID[funcClust[i]])
        else:
            newCluster.append(str(0))
    setCandidateGene=list(set(CandidateGene))
    return(newCluster,setCandidateGene)
def loessRegression(funcGeneGini,funcLogGeneMax,funcGene,funcGiniCutoff):
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
        if (float(pvalue[i])<0.0001 and dictGini[newGene[i]]>=funcGiniCutoff):
            funcSigGiniGeneGini[newGene[i]]=dictGini[newGene[i]]
            funcSigGiniGenePvalue[newGene[i]]=pvalue[i]
    return(funcSigGiniGeneGini,funcSigGiniGenePvalue)

def calGini(adataSC,**kwargs):
    ginicutoff=kwargs.get('ginicutoff', 0.6)
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
    sigGiniGene,sigGiniPvalue=loessRegression(geneGini,logGeneMax,geneLabel,ginicutoff)
    sigGiniListK=sigGiniPvalue.keys()
    sigGiniListV=sigGiniPvalue.values()
    sigGiniList=np.array([list(a) for a in zip(sigGiniListK,sigGiniListV)])
    if (len(sigGiniGene)<2):
        raise SystemExit(str(len(sigGiniGene))+" giniGene found, set a lower value of --giniCutoff.")
    ###select and save high-Gini gene List###
    geneGiniBool=[]
    for i in range(len(geneLabel)):
        value=0
        if geneLabel[i] in sigGiniGene.keys():
            value=1
        geneGiniBool.append(value)
    
    adataSC.var['gini']=np.array(geneGiniBool,dtype=bool)

def clusterGini(adataSC,**kwargs):
    cluster_neighbors=kwargs.get('neighbors', 5)
    cluster_resolution=kwargs.get('resolution', 0.1)
    cluster_method=kwargs.get('method', "leiden")
    adataGini=adataSC[:,adataSC.var['gini']]
    scaleMatrix=arctanTransform((adataGini.X))
    adataScaleGini=anndata.AnnData(X=scaleMatrix)
#    adataGini.X=arctanTransform((adataGini.X))
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

#!/usr/bin/env python
##################################################
# File Name: fano.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Mon 08 Jul 2019 04:08:00 PM EDT
################################################

import numpy as np
import scanpy as sc
import anndata

def calFanoFactor(expArray,geneName):
    fanoDict={}
    exp=np.array(expArray,dtype=float)
    expT=exp.T
    fano=[]
    for i in range(0,len(expT)):
        if (np.mean(expT[i])==0):
            fano.append(0)
        else:
            genefano=np.var(expT[i])/np.mean(expT[i])
            fano.append(genefano)
    for i in range(0,len(geneName)):
        fanoDict[geneName[i]]=fano[i]
    fanoSort=sorted(fanoDict.items(),key=lambda item:item[1],reverse=True)
    highFano={}
    for num in range(0,1000):
        highFano[fanoSort[num][0]]=fanoDict[fanoSort[num][0]]
    return(highFano)

def calFano(adataSC):
    #####louvin cluster by using high Fano genes######
    cellLabel=adataSC.obs.index.tolist()
    geneLabel=adataSC.var.index.tolist()
    sc.pp.log1p(adataSC)
#    sigFanoGene=calFanoFactor(adataSC.X,geneLabel)
#    geneFanoBool=[]
#    for i in range(len(geneLabel)):
#        value=0
#        if geneLabel[i] in sigFanoGene.keys():
#            value=1
#        geneFanoBool.append(value)
#    adataSC.var['fano']=np.array(geneFanoBool,dtype=bool)
    sc.pp.highly_variable_genes(adataSC)

def clusterFano(adataSC,**kwargs):
    cluster_neighbors=kwargs.get('neighbors', 15)
    cluster_resolution=kwargs.get('resolution', 0.1)
    cluster_method=kwargs.get('method', "leiden")
    adataFano = adataSC[:,adataSC.var['highly_variable']]
    fanoGeneLabel=adataFano.var.index.tolist()
    geneDispersions=adataSC.var['highly_variable'].reset_index().values
    sc.pp.scale(adataFano)
    sc.pp.pca(adataFano)
    sc.pp.neighbors(adataFano,n_neighbors=cluster_neighbors)
    if (cluster_method=="louvain"):
        sc.tl.louvain(adataFano,resolution=cluster_resolution)
        adataSC.obs['fano']=adataFano.obs['louvain']
    else:
        sc.tl.leiden(adataFano,resolution=cluster_resolution)
        adataSC.obs['fano']=adataFano.obs['leiden']
    return(adataFano)

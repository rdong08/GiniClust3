#!/usr/bin/env python
##################################################
# File Name: test.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Thu 11 Jul 2019 11:49:02 AM EDT
################################################

import scanpy as sc
import numpy as np
import giniclust3 as gc
import anndata

####Load and filter dataset####
adataRaw=sc.read_csv("./data/GSM1599495_ES_d0_biorep_techrep1.csv",first_column_names=True)
sc.pp.filter_cells(adataRaw,min_genes=3)#####remover gene expressed less than N cell
sc.pp.filter_genes(adataRaw,min_cells=200)#####remove cell express less than M gene
adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)

####GiniIndexClust and FanoFactorClust####
gc.gini.calGini(adataSC)
adataGini=gc.gini.clusterGini(adataSC,neighbors=3)

gc.fano.calFano(adataSC)
adataFano=gc.fano.clusterFano(adataSC)

####ConsensusClust####
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster)
gc.consensus.clusterMtilde(consensusCluster)
gc.consensus.projectFinalCluster(consensusCluster)
np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')

####UMAP visualization####
adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']
gc.plot.umapGini(adataGini)
gc.plot.umapFano(adataFano)


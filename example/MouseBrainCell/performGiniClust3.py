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
import datetime

adataSC=sc.read_h5ad("adata_rmDoublet.h5ad")
starttime=datetime.datetime.now()
gc.gini.calGini(adataSC)
adataGini=gc.gini.clusterGini(adataSC)
gc.fano.calFano(adataSC)
adataFano=gc.fano.clusterFano(adataSC)
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
np.savetxt("rare.txt",consensusCluster['giniCluster'], delimiter="\t",fmt='%s')
np.savetxt("fano.txt",consensusCluster['fanoCluster'], delimiter="\t",fmt='%s')

gc.consensus.generateMtilde(consensusCluster)
gc.consensus.clusterMtilde(consensusCluster)
gc.consensus.projectFinalCluster(consensusCluster)

adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']
fanotime=datetime.datetime.now()
print ("fano time is:"+str((fanotime-starttime).seconds))
np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')
#gc.plot.umapGini(adataGini)
#gc.plot.umapFano(adataFano)


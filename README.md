# GiniClust3
GiniClust3 is a a fast and memory-saving rare cell cluster identification method for large scale single-cell gene expression data.

A schematic flow shows the pipeline
-----------------------------------
![pipeline](https://github.com/rdong08/GiniClust3/blob/master/pipeline.png)

## Prerequisites
* numpy (https://numpy.org/)
* scikit-learn (https://scikit-learn.org/stable/)
* scipy (https://www.scipy.org/)
* statsmodels (https://www.statsmodels.org/stable/index.html)
* anndata (https://anndata.readthedocs.io/en/stable/)
* scanpy (https://scanpy.readthedocs.io/en/stable/)

## Installation
```bash
python setup.py install

## Usage and example:
-----
Import associated packages
```bash
import scanpy as sc
import numpy as np
import giniclust3 as gc
import anndata
```
Read single cell file
```bash
adataRaw=sc.read_csv("./data/GSM1599495_ES_d0_biorep_techrep1.csv",first_column_names=True)
```
Filter and normalization
```bash
sc.pp.filter_cells(adataRaw,min_genes=3)
sc.pp.filter_genes(adataRaw,min_cells=200)
adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)
```

Perform GiniIndexClust
```bash
gc.gini.calGini(adataSC)
adataGini=gc.gini.clusterGini(adataSC,neighbors=3)
```
Perform FanoFactorClust
```bash
gc.fano.calFano(adataSC)
adataFano=gc.fano.clusterFano(adataSC)
```
ConsensusClust
```bash
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster)
gc.consensus.clusterMtilde(consensusCluster)
gc.consensus.projectFinalCluster(consensusCluster)
np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')
```
UMAP visualization
```bash
adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']
gc.plot.umapGini(adataGini)
gc.plot.umapFano(adataFano)
```

Citation
--------

**NA**

License
-------

Copyright (C) 2019 YuanLab.
See the [LICENSE](https://github.com/rdong08/GiniClust3/blob/master/LICENSE)
file for license rights and limitations (MIT).

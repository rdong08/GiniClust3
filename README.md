# GiniClust3
GiniClust3: a fast and memory-efficient tool for rare cell type identification

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

Install by using anaconda (recommend)
```bash
conda install -c rdong giniclust3
```
Download from Github and install
```bash
python setup.py install
```
Install by using pip
```bash
pip install giniclust3
```

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
Filter expression matrix
```bash
sc.pp.filter_cells(adataRaw,min_genes=3)
sc.pp.filter_genes(adataRaw,min_cells=200)
```
Format expression matrix
```bash
###example csv file is col:cells X row:genes. Skip this step if the input matrix is col:genes X row:cells
adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)
```
Normalization
```bash
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)
```

Perform GiniIndexClust
```bash
gc.gini.calGini(adataSC) ###Calculate Gini Index
adataGini=gc.gini.clusterGini(adataSC,neighbors=3) ###Cluster based on Gini Index
```
Perform FanoFactorClust
```bash
gc.fano.calFano(adataSC) ###Calculate Fano factor
adataFano=gc.fano.clusterFano(adataSC) ###Cluster based on Fano factor
```
ConsensusClust
```bash
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster) ###Generate consensus matrix
gc.consensus.clusterMtilde(consensusCluster) ###Cluster consensus matrix
gc.consensus.projectFinalCluster(consensusCluster) ###Projection to each cell
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

# GiniClust3
## GiniClust3: a fast and memory-efficient tool for rare cell type identification
GiniClust is a clustering method specifically designed for rare cell type detection. It uses the Gini index to identify genes that are associated with rare cell types without prior knowledge. This differs from traditional clustering methods using highly variable genes. Using a cluster-aware, weighted consensus clustering approach, we can combine the outcomes from Gini index and Fano factor-based clustering and identify both common and rare cell types. In this new version (GiniClust3), we have substantially increased the speed and reduced memory usage in order to meet the need for large data size. It can now be used to identify rare cell types from over a million cells. Previous versions of GiniClust can be found below: GiniClust (https://github.com/lanjiangboston/GiniClust). GiniClust2 (https://github.com/dtsoucas/GiniClust2).

GiniClust3 documentation is available through https://giniclust3.readthedocs.io/en/latest/, including installation instructions and tutorial.

If you use GiniClust v1.0-v3.0, please consider cite one or more of the following papers:
- Jiang L, Chen H, Pinello L, Yuan GC. GiniClust: detecting rare cell types from single-cell gene expression data with Gini index. Genome Biol. 2016 Jul 1;17(1):144. PMCID:PMC4930624
- Tsoucas D, Yuan GC. GiniClust2: a cluster-aware, weighted ensemble clustering method for cell-type detection. Genome Biol. 2018 May 10;19(1):58. PMCID:PMC5946416
- Dong R, Yuan GC. GiniClust3: a fast and memory-efficient tool for rare cell type identification. BMC Bioinformatics. 2020 Apr 25;21(1):158. PMCID:PMC7183612


A schematic overview of the GiniClust3 pipeline
-----------------------------------
<img src="https://github.com/rdong08/GiniClust3/blob/master/docs/images/pipeline.png" width="500">

## Prerequisites
* numpy (https://numpy.org/)
* scikit-learn (https://scikit-learn.org/stable/)
* scipy (https://www.scipy.org/)
* statsmodels (https://www.statsmodels.org/stable/index.html)
* scanpy (https://scanpy.readthedocs.io/en/stable/)

## Installation
Scanpy is needed to be installed first from "https://scanpy.readthedocs.io/en/stable/installation.html".

### Install by using anaconda (recommend)
```bash
conda install -c rdong giniclust3
```
### OR download from Github and install
```bash
python setup.py install
```
### OR install by using pip
```bash
pip install giniclust3
```

 Usage and example:
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
np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')
```
UMAP visualization
```bash
adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']
gc.plot.plotGini(adataGini)
gc.plot.plotFano(adataFano)
```

Citation
--------

Dong R, Yuan GC. GiniClust3: a fast and memory-efficient tool for rare cell type identification. BMC Bioinformatics. 2020 Apr 25;21(1):158. PMCID:PMC7183612

License
-------

Copyright (C) 2019 YuanLab.
See the [LICENSE](https://github.com/rdong08/GiniClust3/blob/master/LICENSE)
file for license rights and limitations (MIT).

#Usage and example:

-----
Import packages

    import scanpy as sc
    import numpy as np
    import giniclust3 as gc
    import anndata

Read single cell file

    adataRaw=sc.read_csv("./data/GSM1599495_ES_d0_biorep_techrep1.csv",first_column_names=True)

Filter expression matrix (optional)

    sc.pp.filter_cells(adataRaw,min_genes=3)
    sc.pp.filter_genes(adataRaw,min_cells=200)

Transform expression matrix (Skip this step if the input matrix is col:genes X row:cells)

    adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)

Normalization

    sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)

Perform GiniIndexClust

    gc.gini.calGini(adataSC) ###Calculate Gini Index
    adataGini=gc.gini.clusterGini(adataSC,neighbors=3) ###Use higher value of neighbor in larger dataset. Recommend (5:15)

Perform FanoFactorClust

    gc.fano.calFano(adataSC) ###Calculate Fano factor
    adataFano=gc.fano.clusterFano(adataSC) ###Cluster based on Fano factor

ConsensusClust

    consensusCluster={}
    consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
    consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
    gc.consensus.generateMtilde(consensusCluster) ###Generate consensus matrix
    gc.consensus.clusterMtilde(consensusCluster) ###Cluster consensus matrix
    np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')

UMAP visualization

    adataGini.obs['final']=consensusCluster['finalCluster']
    adataFano.obs['final']=consensusCluster['finalCluster']
    gc.plot.umapGini(adataGini)
    gc.plot.umapFano(adataFano)

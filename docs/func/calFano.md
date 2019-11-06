gc.fano.calFano
===============

Highly variable gene selection by calculate Fano factor for each gene.

    gc.fano.calFano(adata,method='scanpy')

**Parameters**

    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    method: string, optional (Default: 'scanpy')
        method='gini2' or method='scanpy'. 'gini2' mode indicates Fano factor is calculated based on var/mean as indicated in GiniClust2. 'scanpy' mode is hvg selection by using scanpy implemented function. Recommend: 'scanpy' mode.

**Returns**

Dictionary with highly variable genes. Store in adata.var['highly_variable'].

**Example**

    gc.fano.calFano(adataSC,,method='scanpy')

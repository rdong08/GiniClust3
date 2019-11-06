gc.fano.clusterFano
==================
Cluster cell based on selected highly variable genes.
    adataFano=gc.fano.clusterFano(adataSC,neighbors=15,resolution=0.1)

**Parameters**

    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    neighbors: int, optional (Default=15)
        The size of local neighborhood used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. This values should be in the range 2 to 100. Recommended neighbors = 15.
    resolution: float, optional (Default=0.1)
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.
    method: string, optional (Default: 'leiden')
        Using Louvain or Leiden to perform cluster on single cell sequencing data. `method='louvain'` or `method='leiden'`.

**Returns**

Returns dictionary with highly variable genes. adata.var['highly_variable']

**Example**

    adataFano=gc.fano.clusterFano(adataSC,neighbors=15)

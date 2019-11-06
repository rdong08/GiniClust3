gc.gini.clusterGini
==================
Cluster cell based on Gini Index value.

    adataGini=gc.gini.clusterGini(adataSC,neighbors=5,resolution=0.1)

**Parameters**

    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    neighbors: int, optional (Default=5)
        The size of local neighborhood used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. For rare cell identification this values should be in the range 2 to 15. Recommended neighbors = 5.
    resolution: float, optional (Default=0.1)
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.
    method: string, optional (Default: 'leiden')
        `method='louvain'` or `method='leiden'`.

**Returns**

Dictionary with gini cluster result. adata.var['rare']

**Example**

    adataGini=gc.gini.clusterGini(adataSC,neighbors=5)

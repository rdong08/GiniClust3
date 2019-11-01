gc.plot.plotFano：
===============

Umap or TSNE visualization based on highly variable genes.

    gc.plot.plotFano(adataFano)

**Parameters**

    adata: Anndata
        Anndata with highly variable gene information.
    method: string, optional (Default: 'umap')
        Choose a method to do the visualization. 'method' should be 'umap' or 'tsne'.

**Returns**

UMAP or TSNE plot in "figure/finalCluster_fano.pdf".

**Example**


    gc.plot.plotFano(adataFano,method='umap')

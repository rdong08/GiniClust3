gc.plot.plotGini：
===============

Umap or TSNE visualization based on Gini genes.

    gc.plot.plotGini(adataGini)

**Parameters**

    adata: Anndata
        Anndata with Gini gene information.
    method: string, optional (Default: 'umap')
        Choose a method to do the visualization. 'method' should be 'umap' or 'tsne'.

**Returns**

UMAP or TSNE plot in "figure/finalCluster_gini.pdf".

**Example**


    gc.plot.plotGini(adataGini,method='umap')

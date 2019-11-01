gc.gini.calGini：
===============
Calculate Gini Index value for each gene and then select highly variable genes.

    gc.gini.calGini(adata,selection='p_value',p_value='0.0001',min_gini_value='0.6')

**Parameters**

    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    selection: string, optional (Default: 'p_value')
        Method for Gini gene selection. `selection='p_value' or selection='gini_value'`. 'p_value' mode indicates Gini genes are selected by p value of Loess regression. 'gini_value' mode indicates Gini genes are selected by Gini Index value. Recommend 'p_value' mode.
    p_value: float, optional (Default=0.0001)
        If `selection='p_value'`, assign a p value cutoff for Gini gene selection.
    min_gini_value: float, optional (Default=0.6)
        Assign min Gini Index value cutoff in Gini gene selection. The Gini Index value is range from 0 to 1. Thus the min_gini_value cutoff is better range from 0 to 0.8. Larger value indicates a more stringent in select Gini genes, smaller values indicates more Gini genes are selected.


**Returns**

Dictionary with Gini genes. Store in adata.var['gini'].

**Example**

    gc.gini.calGini(adataSC,selection='p_value',p_value='0.0001',min_gini_value='0.6')

gc.gini.calGini
===============

```
gc.gini.calGini(adata,selection='p_value',p_value='0.0001',min_gini_value='0.6')
```

    Calculate Gini Index value for each gene.
    Params
    ------
    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    selection: string, optional (Default: 'p_value')
        selection='p_value' or selection='gini_value'. 'p_value' mode indicates
        Gini genes are selected by p value of Loess regression. 'gini_value' mode
        indicates Gini genes are selected by Gini Index value. Recommend 'p_value'
        mode.
    p_value: float, optional (Default=0.0001)
        If `selection='p_value'`, assign a p value cutoff for Gini gene selection.
    min_gini_value: float, optional (Default=0.6)
        Assign min Gini Index value cutoff in Gini gene selection. The Gini
        Index value is range from 0 to 1. Thus the min_gini_value cutoff is
        better in range from 0 to 0.8. Larger value indicates a more stringent
        in select Gini genes, smaller values indicates more Gini genes are
        selected

    A. Gini gene selection in 'p_value' mode: 1. p value < p_value 2. Gini Index
    value >= min_gini_value. (Recommended)
    B. In 'gini_value' mode: Gini Index value > min_gini_value.

    Returns
    -------
    Returns dictionary with Gini genes. adata.var['gini']

#!/usr/bin/env python
##################################################
# File Name: gini.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Mon 08 Jul 2019 03:40:50 PM EDT
################################################

import scanpy as sc

def plotGini(adataGini,**kwargs):
    """
    Umap visualization of rare clusters.
    Params
    ------
    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    method: string, optional (Default: 'umap')
        'method' should be 'umap' or 'tsne'.

    Return
    -------
    UMAP or T-SNE plot in "figure/finalCluster_gini.pdf"
    """

    method=kwargs.get('method','umap')
    if (method=='umap'):
        sc.tl.umap(adataGini)
        sc.pl.umap(adataGini,color='final',save="_finalCluster_umap_gini.pdf")
    elif(method=='tsne'):
        sc.tl.tsne(adataGini)
        sc.pl.tsne(adataGini,color='final',save="_finalCluster_tsne_gini.pdf")
    else:
        raise SystemExit("'method' should be 'umap' or 'tsne'.")   

def plotFano(adataFano,**kwargs):
    """
    Umap visualization of rare clusters.
    Params
    ------
    adata: Anndata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    method: string, optional (Default: 'umap')
        'method' should be 'umap' or 'tsne'.

    Return
    -------
    UMAP or T-SNE plot in "figure/umap_finalCluster_based_on_fano.pdf"
    """
    method=kwargs.get('method','umap')
    if (method=='umap'):
        sc.tl.umap(adataFano)
        sc.pl.umap(adataFano,color='final',save="_finalCluster_umap_fano.pdf")
    elif(method=='tsne'):
        sc.tl.tsne(adataFano)
        sc.pl.tsne(adataFano,color='final',save="_finalCluster_tsne_fano.pdf")
    else:
        raise SystemExit("'method' should be 'umap' or 'tsne'.")

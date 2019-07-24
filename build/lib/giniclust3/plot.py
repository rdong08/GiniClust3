#!/usr/bin/env python
##################################################
# File Name: gini.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Mon 08 Jul 2019 03:40:50 PM EDT
################################################

import scanpy as sc

def umapGini(adataGini,**kwargs):
    sc.tl.umap(adataGini)
    sc.pl.umap(adataGini,color='final',save="_finalCluster_based_on_gini.pdf")

def umapFano(adataFano,**kwargs):
    sc.tl.umap(adataFano)
    sc.pl.umap(adataFano,color='final',save="_finalCluster_based_on_fano.pdf")

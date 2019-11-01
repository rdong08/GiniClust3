gc.consensus.clusterMtilde：
===============

Cluster generated consensus matrix 

    gc.consensus.clusterMtilde(GCdict)

**Parameters**

    GCdict
        Dictionary with Gini, Fano cluster results and consensus matrix.
    k: int, optional
        Number of K in KMeans clustering. Default value based on Gini and Fano cluster results.

**Returns**

Returns dictionary 'GCdict['finalCluster']' with consensus matrix Mtilde and associated information.

**Example**

    gc.consensus.clusterMtilde(GCdict)

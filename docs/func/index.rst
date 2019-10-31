.. module:: giniclust3
.. automodule:: giniclust3
   :noindex:

Functions
===


Import GiniClust3 as::

   import giniclust3 as gc

GiniIndexClust: gc.gini
-------------------

Gini Index calculation and high Gini gene selection: `calGini`
~~~~~~~~~~~~~~~~~~~

.. module:: gc.gini.calGini
.. currentmodule:: giniclust3

High Gini genes identification by using loess regression. 

Rare cell identification based on Gini Index: `clusterGini`
~~~~~~~~~~~~~~~~~~~

.. module:: gc.gini.clusterGini
.. currentmodule:: giniclust3

Rare cell cluster identification based on Gini gene expression matrix. 

FanoIndexClust: gc.fano
-------------------

Highly variable gene selection: `calFano`
~~~~~~~~~~~~~~~~~~~

.. module:: gc.fano.calFano
.. currentmodule:: giniclust3

Highly variable gene identification. 

Rare cell identification based on Gini Index: `clusterGini`
~~~~~~~~~~~~~~~~~~~

.. module:: gc.fano.clusterFano
.. currentmodule:: giniclust3

Cluster cells based on selection highly variable gene expression matrix. 


ConsensusClust
-------------------

Generate consensus matrix (Mtilde): `gc.consensus.generateMtilde`
~~~~~~~~~~~~~~~~~~~

.. module:: gc.consensus.generateMtilde
.. currentmodule:: giniclust3

Generate consensus matrix based on Gini and Fano cluster results. Same as GiniClust2.

Final cluster: `gc.consensus.clusterMtilde`
~~~~~~~~~~~~~~~~~~~

Finally cluster consensus matrix by using K-Means.

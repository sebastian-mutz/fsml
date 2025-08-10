---
title: NLP: Nonlinear Procedures
---

# Overview

This is the API documentation for all nonlinear procedures.

[TOC]


<br>
# Hierarchical Clustering

## `fsml_hcluster`

### Description
The procedure is an implementation of the agglomerative hierarchical clustering method
that groups data points into clusters by iteratively merging the most similar clusters.
The procedure uses centroid linkage and the Mahalanobis distance as a measure of similarity.

The input matrix (`x`) holds observations in rows (`nd`) and variables in columns (`nv`).
The target number of clusters (`nc`) must be at least *1* and not greater than the number
of data points.

The variables are standardised before computing the covariance matrix on the transformed
data. The matrix is used for calculating the Mahalanobis distance.

Clusters are merged iteratively until the target number of clusters is reached.
The global mean (`gm`), cluster centroids (`cm`), membership assignments (`cl`),
and cluster sizes (`cc`), the covariance matrix (`cov`) and standard deviations
(`sigma`) used in the distance calculations are returned.

### Syntax
`call` [[fsml(module):s_nlp_cluster_h(interface)]]`(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`nc`: A scalar of type `integer`.

### Returns
`gm`: A rank-1 array of type `real` with dimension `nv`.

`cm`: A rank-2 array of type `real`  with dimensions `nv`, `nv`.

`cl`: A rank-1 array of type `integer` with dimension `nd`.

`cc`: A rank-1 array of type `integer` with dimension `nc`.

`cov`: A rank-2 array of type `real` with dimensions `nv`, `nv`.

`sigma`: A rank-1 array of type `real` with dimension `nv`.


<br>
# K-Means Clustering

## `fsml_kmeans`

### Description
The procedure implements the K-means clustering algorithm using the Mahalanobis
distance as the similarity measure. It accepts initial centroids (`cm_in`), refines
them iteratively, and returns the final centroids (`cm`).

The input matrix (`x`) holds observations in rows (`nd`) and variables in columns (`nv`).
The number of clusters (`nc`) must be at least *1* and not greater than the number of
data points. The procedure assigns each observation to the nearest centroid using the
Mahalanobis distance, recomputes centroids from cluster memberships, and iterates until
convergence or the iteration limit is reached. Final centroids are sorted by the first
variable, and assignments are updated accordingly.

If the covariance matrix (`cov_in`) is passed, it will be used to calculate the
Mahalanobis distance. If it is not passed, the variables are standardised before
computing the covariance matrix on the transformed data.

The global mean (`gm`), cluster centroids (`cm`), membership assignments (`cl`),
and cluster sizes (`cc`), the covariance matrix (`cov` - either `cov_in` or internally
calculated) and standard deviations (`sigma`) used in the distance calculations are returned.

### Syntax
`call` [[fsml(module):fsml_kmeans(interface)]]`(x, nd, nv, nc, cm_in, gm, cm, cl, cc, cov, sigma, [cov_in])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`nc`: A scalar of type `integer`.

`cm_in`: A rank-2 array of type `real` with dimensions `nv`, `nc`.

`cov_in` (optional): A rank-2 array of type `real` with dimensions `nv`, `nv`.

### Returns
`gm`: A rank-1 array of type `real` with dimension `nv`.

`cm`: A rank-2 array of type `real` with dimensions `nv`, `nc`.

`cl`: A rank-1 array of type `integer` with dimension `nd`.

`cc`: A rank-1 array of type `integer` with dimension `nc`.

`cov`: A rank-2 array of type `real` with dimensions `nv`, `nv`.

`sigma`: A rank-1 array of type `real` with dimension `nv`.


<br>
# Examples

```fortran
{!example/example_nlp.f90!}
```


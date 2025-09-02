---
title: NLP: Non-Linear Procedures
---

# Overview

This is the API documentation for all nonlinear procedures.

[TOC]


<br>
# Hierarchical Clustering

## `fsml_hclust`

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

**Note:** This procedure uses the pure procedure for calculating the Mahalanobis distance
`f_lin_mahalanobis_core`, which uses`chol` from the `stdlib_linalg` module.

@note
The current implementation of the procedure is exact, but very slow for large `nd`.
@endnote

### Syntax
`call` [[fsml(module):fsml_hclust(interface)]]`(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`nc`: A scalar of type `integer`.

Invalid argument values will result in the return of a sentinel value.

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
The procedure implements the k-means clustering algorithm using the Mahalanobis
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
and cluster sizes (`cc`), the covariance matrix (`cov` - either `cov_in` if passed,
or internally calculated if `cov_in` is not passed) and standard deviations (`sigma`)
used in the distance calculations are returned.

**Note:** This procedure uses the pure procedure for calculating the Mahalanobis distance
`f_lin_mahalanobis_core`, which uses`chol` from the `stdlib_linalg` module.

### Syntax
`call` [[fsml(module):fsml_kmeans(interface)]]`(x, nd, nv, nc, cm_in, gm, cm, cl, cc, cov, sigma, [cov_in])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`nc`: A scalar of type `integer`.

`cm_in`: A rank-2 array of type `real` with dimensions `nv`, `nc`.

`cov_in` (optional): A rank-2 array of type `real` with dimensions `nv`, `nv`.

Invalid argument values will result in the return of a sentinel value.

### Returns
`gm`: A rank-1 array of type `real` with dimension `nv`.

`cm`: A rank-2 array of type `real` with dimensions `nv`, `nc`.

`cl`: A rank-1 array of type `integer` with dimension `nd`.

`cc`: A rank-1 array of type `integer` with dimension `nc`.

`cov`: A rank-2 array of type `real` with dimensions `nv`, `nv`.

`sigma`: A rank-1 array of type `real` with dimension `nv`.


<br>
# Hybrid Hierarchical/K-Means Clustering

## `fsml_hkmeans`

### Description
The procedure implements a hybrid clustering approach combining agglomerative hierarchical
clustering and k-means clustering, both using the Mahalanobis distance as the similarity measure.
The hierarchical step first partitions the data into `nc` clusters by iteratively merging the most
similar clusters. The resulting centroids from are then used as initial centroids (`cm_in`)
for the k-means procedure, which refines them iteratively.

The input matrix (`x`) holds observations in rows (`nd`) and variables in columns (`nv`).
The number of clusters (`nc`) must be at least *1* and not greater than the number of data points.
In the hierarchical clustering step, variables are standardised before computing the covariance matrix
on the transformed data. The covariance matrix is passed to the k-means clustering procedure along
with the initial cluster centroids. The k-means clustering step then assigns each observation to the
nearest centroid, recomputes centroids from cluster memberships, and iterates until convergence or
the iteration limit is reached. Final centroids are sorted by the first variable, and assignments
are updated accordingly.

The global mean (`gm`), final cluster centroids (`cm`), membership assignments (`cl`), and cluster
sizes (`cc`), the covariance matrix (`cov`) and standard deviations (`sigma`) used in the distance
calculations are returned.

**Note:** This procedure uses the pure procedure for calculating the Mahalanobis distance
`f_lin_mahalanobis_core`, which uses`chol` from the `stdlib_linalg` module.

### Syntax
`call` [[fsml(module):fsml_hkmeans(interface)]]`(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`. It must be at least *1*.

`nv`: A scalar of type `integer`. It must be at least *1*.

`nc`: A scalar of type `integer`. It must be at least *1* and not greater than `nd`.

Invalid argument values will result in the return of a sentinel value.

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
{!example/modules/example_nlp.f90!}
```


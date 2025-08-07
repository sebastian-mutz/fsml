---
title: LIN: Linear Procedures
---

# Overview

This is the API documentation for all linear procedures (relying heavily on linear algebra).

[TOC]


<br>
# Principal Component Analysis

## `fsml_pca`

### Description
Principal Component Analysis (PCA) is a procedure to reduce the dimensionality
of multivariate data by identifying a set of orthogonal vectors (eigenvectores)
that represent directions of maximum variance in the dataset.

The procedure `fsml_pca` is a wrapper for `fsml_eof` and offers a simpler,
more familiar interface for non-geoscientists. The EOF interface allows for
more options to be passed that are irrelevant to standard applications of PCA.
The PCA procedure calls the EOF procedures with weights (`wt`) set to *1.0*,
and matrix options set to `opt = 0` to force the use of the covariance matrix
to be comparable to other common implementations of a PCA (e.g., sklearn).

The covariance matrix \( \mathbf{C} \) is computed as:
$$
\mathbf{C} = \frac{1}{m - 1} \mathbf{X}^\top \mathbf{X}
$$
where \( \mathbf{X} \) is the preprocessed (centred and optionally standardised) data matrix,
and \( m \) is the number of observations (rows in `x`).

A symmetric eigen-decomposition is then performed:
$$
\mathbf{C} \mathbf{E} = \mathbf{E} \Lambda
$$
where \( \mathbf{E} \) contains the EOFs (`ev`), and \( \Lambda \) is a diagonal matrix
of eigenvalues (`ew`).

The principal components or scores (PCs, `pc`) are given by:
$$
\mathbf{PC} = \mathbf{X} \mathbf{E}
$$
The number of valid PC modes is determined by the number of non-zero eigenvalues.
Arrays are initialised to zero and populated only where eigenvalues are strictly positive.

The explained variance (`r2`) for each component is computed as a fraction:
$$
r^2_j = \frac{\lambda_j}{\sum_k \lambda_k}
$$
where \( j \) is the PC index, and \( k \) spans all retained eigenvalues,
representing all principal components that explain variability in the data.

**Note:** This subroutine uses `eigh` from the `stdlib_linalg` module to compute
eigenvalues and eigenvectors of the symmetric covariance matrix.

### Syntax
`call` [[fsml(module):fsml_pca(interface)]]`(x, m, n, pc, ev, ew [, r2])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `m`, `n`.

`m`: A scalar of type `integer`.

`n`: A scalar of type `integer`.

Invalid argument values will result in the return of a sentinel value.

### Returns
`pc`: A rank-2 array of type `real` (same type as `x`) with dimensions `m`, `n`.

`ev`: A rank-2 array of type `real` (same type as `x`) with dimensions `n`, `n`.

`ew`: An rank-1 array of type `real`. It stores the eigenvalues.

`r2`: An optional return and rank-1 array of type `real`.


<br>
# Empirical Orthogonal Functions

## `fsml_eof`

### Description
Empirical Orthogonal Function (EOF) analysis is a procedure to reduce the dimensionality
of multivariate data by identifying a set of orthogonal vectors (EOFs or eigenvectores)
that represent directions of maximum variance in the dataset.
The term *EOF analysis* is often used interchangably with the geographically weighted
principal component analysis (PCA). The procedures are mathematically equivalent, but
procedures for EOF analysis offer some additional options that are mostly relevant for
geoscience. The procedure `fsml_pca` is a wrapper for `fsml_eof` that offers a simpler,
more familiar interface for non-geoscientists.

For a classic EOF analysis, the input matrix `x` holds data or observations that have been
discretised in time and space. Rows (`m`) and columns (`n`) can therefore be interpreted
as time and space dimensions, respectively. EOF analysis allows for geographical weighting,
which translates to column-wise weighting prior to analysis in the procedure.
Weights can be set by bassing the rank-1 array `wt` of dimension `n`. If this optional
argument is not passed, the procedure will default to equal weights of value \( wt=1/n \).
It is numerically more stable than *1.0*, which is the default for many implementations of a PCA.

After the weighting is applied, the covariance or correlation matrix \( \mathbf{C} \) is computed:
$$
\mathbf{C} = \frac{1}{m - 1} \mathbf{X}^\top \mathbf{X}
$$
where \( \mathbf{X} \) is the preprocessed (centred and optionally standardised) data matrix,
and \( m \) is the number of observations (rows in `x`).
The value of the optional argument `opt` determines if the covariance matrix (`opt = 0`) or
correlation matrix (`opt = 1`) is constructed. If the argument is not passed, the procedure will
default to the use of the covariance matrix, as is the standard for a regular PCA.

A symmetric eigen-decomposition is then performed:
$$
\mathbf{C} \mathbf{E} = \mathbf{E} \Lambda
$$
where \( \mathbf{E} \) contains the EOFs (`eof`), and \( \Lambda \) is a diagonal matrix
of eigenvalues (`ew`).

The principal components or scores (PCs, `pc`) are given by:
$$
\mathbf{PC} = \mathbf{X} \mathbf{E}
$$
The number of valid EOF/PC modes is determined by the number of non-zero eigenvalues.
Arrays are initialised to zero and populated only where eigenvalues are strictly positive.

The explained variance (`r2`) for each component is computed as a fraction:
$$
r^2_j = \frac{\lambda_j}{\sum_k \lambda_k}
$$
where \( j \) is the PC index, and \( k \) spans all retained eigenvalues,
representing all principal components that explain variability in the data.

EOFs may optionally be scaled (`eof_scaled`) for more convenient plotting:
$$
\text{EOF}_{\text{scaled}} = \text{EOF} \cdot \sqrt{\lambda_j}
$$

**Note:** This subroutine uses `eigh` from the `stdlib_linalg` module to compute
eigenvalues and eigenvectors of the symmetric covariance or correlation matrix.

### Syntax
`call` [[fsml(module):fsml_eof(interface)]]`(x, m, n, pc, eof, ew [, opt, wt, eof_scaled, r2])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `m`, `n`.

`m`: A scalar of type `integer`.

`n`: A scalar of type `integer`.

`opt`: An optional argument and scalar of type `integer`. It must be *0* or *1*. If not passed, it will default to *0*.

`wt`: An optional argument and rank-1 array of type `real`. If not passed, it will default to equal weights of value *1/n*.

Invalid argument values will result in the return of a sentinel value.

### Returns
`pc`: A rank-2 array of type `real` (same type as `x`) with dimensions `m`, `n`.

`eof`: A rank-2 array of type `real` (same type as `x`) with dimensions `n`, `n`.

`ew`: An rank-1 array of type `real`. It stores the eigenvalues.

`r2`: An optional return and rank-1 array of type `real`.

`eof_scaled`: An optional return and rank-2 array of type `real` (same type as `x`) with dimensions `n`, `n`.


<br>
# Linear Discriminant Analysis (2-Class)

## `fsml_lda_2class`

### Description
The 2-class multivariate Linear Discriminant Analysis (LDA) is a statistical
procedure for classification and the investigation and explanation of differences
between two groups (or classes) with regard to their attribute variables.
It quantifies the discriminability of the groups and the contribution of each of
the attribute variables to this discriminability.

The procedure finds a discriminant function that best separates the two groups.
The function can be expressed as a linear combination of the attribute variables:

$$
Y = \nu_0 + \nu_1 X_1 + \nu_2 X_2 + \dots + \nu_m X_m + \dots + \nu_M X_M
$$

where \( Y \) is the discriminant function, \( X_m (m=1...M) \) are the attribute
variables used in evaluating the differences between the groups, \( nu_m (m=1...M) \)
are the discriminant coefficients associated with each variable, \( M \) (`nv`) is
the number of variables, and \( \nu_0 \) is the y-intercept.
(**Note**: Mathematically, it is analogous to a multivariate linear regression function.)

Each attribute variable \( X_m \) contains elements \( x_{mn} (n=1…N) \) (`x`), where
\( N \) (`nd`) is the number of elements in each group. Each element is associated with a
discriminant value \( y_n \) described by:

$$
y_n = \nu_1 x_{1n} + \nu_2 x_{2n} + \dots + \nu_m x_{mn} + \dots + \nu_M x_{Mn}
$$

Geometrically, this can be visualised as elements \( y_n \) being projected on the
discriminant axis \( Y \). The optimal discriminant function is then determined by
finding an axis, on which the projected elements for the two groups are best separated.
The best separation is given by maximising the discriminant criterion \( \Gamma \) (`g`),
a signal to noise ratio, so that:

$$
\Gamma = \frac{\text{scatter between groups }}{\text{scatter within groups }}
= \frac{(\bar{y}_{G1} - \bar{y}_{G2})^2}
{\sum_{j=1}^{n_1} (y_{G1j} - \bar{y}_{G1})^2 + \sum_{j=1}^{n_2} (y_{G2j} - \bar{y}_{G2})^2}
\rightarrow \max
$$

where \( n_1 \) and \( n_2 \) are the number of elements in groups \( G1 \) and \( G2 \),
respectively. The procedure assumes that these are the same (`nd`) and only accepts 2 groups (`nc = 2`).

The discriminant coefficients are then standardised (`sa`) using the standard deviations
of respective variables. The discriminant function represents a model that best seperates
the groups and can be used as a classification model. The skill of that model is determined
by forgetting the association of each element with the groups and using the model to reclassify
the elements. The score (`score`) is the fraction of correct classifications and can be
interpreted as a measure of how well the function works as a classification model.

The procedure optionally returns the Mahalanobis distance (`mh`) as a measure of distance
between the groups.

**Note:** This subroutine uses `eigh` from the `stdlib_linalg` module to compute
eigenvalues and eigenvectors of the symmetric covariance or correlation matrix.

### Syntax

`call` [[fsml(module):fsml_lda_2class(interface)]]`(x, nd, nv, nc, sa, g, score [, mh])`

### Parameters
`x`: A rank-3 array of type `real` with dimensions `nd`, `nv`, `nc`.

`nd`: A scalar of type `integer`.

`nc`: A scalar of type `integer`. It must be *2*.

`nv`: A scalar of type `integer`.

Invalid argument values will result in the return of a sentinel value.

### Returns
`sa`: A rank-1 array of type `real` and dimension `nv`.

`g`: A scalar of type `real`.

`score`: A rank-3 array of type `real`.

`mh`: An optional return and scalar of type `real`


<br>
# Examples

```fortran
{!example/example_lin.f90!}
```


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

For a classic PCA, the input matrix (`x`) holds data or observations in rows (`nd`)
and different variables in columns (`nv`).

The procedure `fsml_pca` is a wrapper for `fsml_eof` and offers a simpler,
more familiar interface for non-geoscientists. The EOF interface allows for
more options to be passed that are irrelevant to standard applications of PCA.
The PCA procedure calls the EOF procedures with weights (`wt`) set to *1.0*,
and matrix options set to `opt = 0` to force the use of the covariance matrix
to be comparable to other common implementations of a PCA (e.g., Python sklearn).

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

**Note:** This procedure uses `eigh` from the `stdlib_linalg` module to compute
eigenvalues and eigenvectors of the symmetric covariance matrix.

@note
The procesure has no pure equivalent, because it relies on an
impure procedure for eigen-decomposition (`eigh`).
@endnote

### Syntax
`call` [[fsml(module):fsml_pca(interface)]]`(x, nd, nv, pc, ev, ew [, r2])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

Invalid argument values will result in the return of a sentinel value.

### Returns
`pc`: A rank-2 array of type `real` (same type as `x`) with dimensions `nd`, `nv`.

`ev`: A rank-2 array of type `real` (same type as `x`) with dimensions `nv`, `nv`.

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
discretised in time and space. Rows (`nd`) and columns (`nv`) can therefore be interpreted
as time and space dimensions, respectively. EOF analysis allows for geographical weighting,
which translates to column-wise weighting prior to analysis in the procedure.
Weights can be set by bassing the rank-1 array `wt` of dimension `nv`. If this optional
argument is not passed, the procedure will default to equal weights of value \( wt=1/n \).
It is numerically more stable than *1.0*, which is the default for many implementations of a PCA.

After the weighting is applied, the covariance or correlation matrix \( \mathbf{C} \) is computed:
$$
\mathbf{C} = \frac{1}{m - 1} \mathbf{X}^\top \mathbf{X}
$$
where \( \mathbf{X} \) is the preprocessed (centred and optionally standardised) data matrix,
and \( m \) is the number of observations (rows `nd` in `x`).
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

**Note:** This procedure uses `eigh` from the `stdlib_linalg` module to compute
eigenvalues and eigenvectors of the symmetric covariance or correlation matrix.

@note
The procesure has no pure equivalent, because it relies on an
impure procedure for eigen-decomposition (`eigh`).
@endnote

### Syntax
`call` [[fsml(module):fsml_eof(interface)]]`(x, nd, nv, pc, eof, ew [, opt, wt, r2, eof_scaled])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`opt`: An optional argument and scalar of type `integer`. It must be *0* or *1*. If not passed, it will default to *0*.

`wt`: An optional argument and rank-1 array of type `real`. If not passed, it will default to equal weights of value *1/n*.

Invalid argument values will result in the return of a sentinel value.

### Returns
`pc`: A rank-2 array of type `real` (same type as `x`) with dimensions `nd`, `nv`.

`eof`: A rank-2 array of type `real` (same type as `x`) with dimensions `nv`, `nv`.

`ew`: An rank-1 array of type `real`. It stores the eigenvalues.

`r2`: An optional return and rank-1 array of type `real`.

`eof_scaled`: An optional return and rank-2 array of type `real` (same type as `x`)
with dimensions `nv`, `nv`.


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

**Note:** This procedure uses `eigh` from the `stdlib_linalg` module.

@note
The procesure has no pure equivalent, because it relies on an
impure procedure for eigen-decomposition (`eigh`).
@endnote

### Syntax
`call` [[fsml(module):fsml_lda_2class(interface)]]`(x, nd, nv, nc, sa, g, score [, mh])`

### Parameters
`x`: A rank-3 array of type `real` with dimensions `nd`, `nv`, `nc`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`nc`: A scalar of type `integer`. It must be *2*.

Invalid argument values will result in the return of a sentinel value.

### Returns
`sa`: A rank-1 array of type `real` and dimension `nv`.

`g`: A scalar of type `real`.

`score`: A rank-3 array of type `real`.

`mh`: An optional return and scalar of type `real`


<br>
# Ordinary Least Squares Regression

## `fsml_ols`

### Description
The multiple linear Ordinary Least Squares (OLS) regression models the relationship
or linear dependence between a dependent (predictand) variable and and one or more
independent (predictor) variables. The procedure estimates the linear regression
coefficients by minimising the sum of squared residuals.

The estimated regression model is of the form:

$$
y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_m x_m + \dots + \beta_M x_M
$$

where \( y \) is the predictand variable, \( x_m \ (m = 1 \dots M) \) are the predictor variables (`x`)
with `nd` observations, \( \beta_0 \) is the y-intercept (`b0`), \( \beta_m \ (m = 1 \dots M) \) (`b`)
are the regression coefficients, and \( M \) (`nv`) is the number of predictors (excluding the intercept).

The subroutine constructs a full matrix internally by prepending a column of ones to account for
the intercept. The regression coefficients are estimated as:

$$
\hat{\beta} = (X^\top X)^{-1} X^\top y
$$

where \( X \) is the extended design matrix including the intercept term.

The coefficient of determination \( R^2 \) (`r2`) which represents the proportion of the total
variance of \( y \)  (`y`) explained by the predictors. The predicted values (`y_hat`),
standard errors (`se`) of the coefficients, and the covariance matrix of the predictors (`cov_b`)
can optionally be returned by the procedure, too.

**Note:** This procedure uses `eigh` from the `stdlib_linalg` module.

**Note:** The intercept and predictor coefficients are computed separately and returned explicitly.

@note
The procesure has no pure equivalent, because it relies on an
impure procedure for eigen-decomposition (`eigh`).
@endnote

### Syntax
`call` [[fsml(module):fsml_ols(interface)]]`(x, y, nd, nv, b0, b, r2 [, y_hat, se, covb])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`y`: A rank-1 array of type `real` with dimension `nd`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

Invalid argument values will result in the return of a sentinel value.

### Returns
`b0`: A scalar of type `real`.

`b`: A rank-1 array of type `real` with dimension `nv`.

`r2`: A scalar of type `real`.

`y_hat` An optional return and rank-1 array of type `real` with dimension `nd`.

`se` An optional return and rank-1 array of type `real` with dimension `nv`.

`cov_b` An optional return and rank-2 array of type `real` with dimensions `nv`, `nv`.


<br>
# Ridge Regression

## `fsml_ridge`

### Description
The multiple linear Ridge regression models the relationship
or linear dependence between a dependent (predictand) variable and one or more
independent (predictor) variables, incorporating a penalty term on the size of the
regression coefficients to reduce multicollinearity and overfitting.

The procedure estimates the linear regression coefficients by minimising the sum of squared residuals plus a penalty proportional to the square of the magnitude of coefficients:

$$
\hat{\beta} = (X^\top X + \lambda I)^{-1} X^\top y
$$

where \( \lambda \) (`lambda`) is the ridge penalty parameter, and \( I \) is the
identity matrix with the first diagonal element corresponding to the intercept set to zero
(no penalty on intercept).

The estimated regression model is of the form:

$$
y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_m x_m + \dots + \beta_M x_M
$$

where \( y \) is the predictand variable, \( x_m \ (m = 1 \dots M) \) are the predictor variables (`x`)
with `nd` observations, \( \beta_0 \) is the y-intercept (`b0`), \( \beta_m \ (m = 1 \dots M) \) (`b`)
are the ridge regression coefficients, and \( M \) (`nv`) is the number of predictors
(excluding the intercept).

The subroutine constructs a full matrix internally by prepending a column of ones to account for
the intercept. The coefficient of determination \( R^2 \) (`r2`), predicted values (`y_hat`),
ridge-adjusted standard errors (`se`) of the coefficients, and the ridge-adjusted covariance
matrix of the predictors (`cov_b`) can optionally be returned. The covariance matrix and standard
errors are adjusted for the ridge penalty as:

$$
\mathrm{cov}(\hat{\beta}) = \sigma^2 (X^\top X + \lambda I)^{-1} X^\top X (X^\top X + \lambda I)^{-1}
$$

where \( \sigma^2 \) is the residual variance estimate.

**Note:** This procedure uses `eigh` from the `stdlib_linalg` module.

@note
The procesure has no pure equivalent, because it relies on an
impure procedure for eigen-decomposition (`eigh`).
@endnote

### Syntax
`call` [[fsml(module):fsml_ridge(interface)]]`(x, y, nd, nv, lambda, b0, b, r2 [, y_hat, se, covb])`

### Parameters
`x`: A rank-2 array of type `real` with dimensions `nd`, `nv`.

`y`: A rank-1 array of type `real` with dimension `nd`.

`nd`: A scalar of type `integer`.

`nv`: A scalar of type `integer`.

`lambda`: A scalar of type `real`. Must be non-zero positive.

Invalid argument values will result in the return of a sentinel value.

### Returns
`b0`: A scalar of type `real`.

`b`: A rank-1 array of type `real` with dimension `nv`.

`r2`: A scalar of type `real`.

`y_hat`: An optional return and rank-1 array of type `real` with dimension `nd`.

`se`: An optional return and rank-1 array of type `real` with dimension `nv`.

`cov_b`: An optional return and rank-2 array of type `real` with dimensions `nv`, `nv`.


<br>
# Mahalanobis Distance

## `fsml_mahalanobis`

### Description
Computes the Mahalanobis distance between two input feature vectors `x` and `y`.
If a covariance matrix `cov` is provided, it is used directly in the calculation.
Otherwise, the procedure estimates the covariance matrix from the two-sample dataset
formed by `x` and `y`. A Cholesky-based solver is used to perform the distance
calculation.

The Mahalanobis distance is defined as:

$$
D_M(x, y) = \sqrt{ (x - y)^\top \Sigma^{-1} (x - y) }
$$

where \( \Sigma \) is the covariance matrix. The inverse is applied via
the Cholesky decomposition for numerical stability.

**Note:** If passed, the covariance matrix (`cov`) must be positive definite
for the factorisation to succeed.

**Note:** This procedure uses `chol` from the `stdlib_linalg` module.

### Syntax
`dist =` [[fsml(module):fsml_mahalanobis(interface)]]`(x, y [, cov])`

### Parameters
`x`: A rank-1 array of type `real` with dimension `m`.

`y`: A rank-1 array of type `real` with dimension `m`.

`cov`: An optional rank-2 array of type `real` with dimensions `m`, `m`.

### Returns
`dist`: A scalar of type `real`.


<br>
# Examples

```fortran
{!example/example_lin.f90!}
```


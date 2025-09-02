---
title: STS: Basic Statistics
---

# Overview

This is the API documentation for all population and sample statistics procedures.

[TOC]


<br>
# Mean

## `fsml_mean`

### Description
The procedure computes the arithmetic mean.
$$ \bar{x} = \frac{1}{n} \cdot \sum_{i=1}^{n} x_i $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \bar{x} \) is the arithmetic mean of `x`.

### Syntax
`result =` [[fsml(module):fsml_mean(interface)]]`(x)`

### Parameters
`x`: A rank-1 array of type `real`. Must have a length of at least 2.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Median

## `fsml_median`

### Description
The procedure computes median of vector `x` and handles tied ranks.

### Syntax
`result =` [[fsml(module):fsml_median(interface)]]`(x)`

### Parameters
`x`: A rank-1 array of type `real`. Must have a length of at least 2.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.

<br>
# Variance

## `fsml_var`

### Description
The procedure computes the population or sample variance (depending on passed arguments).
$$ \operatorname{var}(x) = \frac{1}{n - \nu} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`,
\( \nu \) (`ddf`) is a degrees of freedom adjustment
(`ddf = 0.0` for population variance, `ddf = 1.0` for sample variance), and
\( \bar{x} \) is the arithmetic mean of `x`.

### Syntax
`result =` [[fsml(module):fsml_var(interface)]]`(x [,ddf])`

### Parameters
`x`: A rank-1 array of type `real`.

`ddf`: An optional argument and scalar of type `real`. If passed, it must be either *0.0* or *1.0*. It will default to *0.0* if not passed.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Standard Deviation

## `fsml_std`

### Description
The procedure computes the population or sample standard deviation (depending on passed arguments).
$$ \sigma = \sqrt{\operatorname{var}(x)} = \sqrt{ \frac{1}{n - \nu} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 } $$
where \( \operatorname{var}(x) \) is the variance of vector `x`,
\( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \nu \) (`ddf`) is a degrees of freedom adjustment
(`ddf = 0.0` for population variance, `ddf = 1.0` for sample variance).

### Syntax
`result =` [[fsml(module):fsml_std(interface)]]`(x [,ddf])`

### Parameters
`x`: A rank-1 array of type `real`.

`ddf`: An optional argument and scalar of type `real`. If passed, it must be either *0.0* or *1.0*. It will default to *0.0* if not passed.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Covariance

## `fsml_cov`

### Description
The procedure computes the population or sample covariance (depending on passed arguments).
$$ \operatorname{cov}(x, y) = \frac{1}{n - \nu} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
where \( n \) is the size of (or number of observations in) vectors `x` and `y`,
\( x_i \) and \( y_i \) are individual elements in `x` and `y`,
\( \nu \) (`ddf`) is a degrees of freedom adjustment
(`ddf = 0.0` for population variance, `ddf = 1.0` for sample variance), and
\( \bar{x} \) and \( \bar{y} \) are the arithmetic means of `x` and `y`.

Vectors `x` and `y` must be the same size.

### Syntax
`result =` [[fsml(module):fsml_cov(interface)]]`(x, y [,ddf])`

### Parameters
`x`: A rank-1 array of type `real`. It must be the same size as `y`.

`y`: A rank-1 array of type `real`. It must be the same size as `x`.

`ddf`: An optional argument and scalar of type `real`. If passed, it must be either *0.0* or *1.0*. It will default to *0.0* if not passed.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x` and `y`.


<br>
# Linear Trend

## `fsml_trend`

### Description
The procedure computes the regression coefficient/trend.
$$ m = \frac{\operatorname{cov}(x, y)}{\operatorname{var}(x)} $$
where \( m \) is the slope of the regression line (linear trend),
\( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
\( \operatorname{var}(x) \) is the variance of `x`.

Vectors `x` and `y` must be the same size.

### Syntax
`result =` [[fsml(module):fsml_trend(interface)]]`(x, y)`

### Parameters
`x`: A rank-1 array of type `real`. It must be the same size as `y`.

`y`: A rank-1 array of type `real`. It must be the same size as `x`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x` and `y`.


<br>
# Pearson Correlation Coefficient

## `fsml_pcc`

### Description
The procedure computes the Pearson correlation coefficient (PCC).
$$ \rho_{x,y} = \frac{\operatorname{cov}(x, y)}{\sigma_x \cdot \sigma_y} $$
where \( \rho_{x,y} \) is the Pearson correlation coefficient for vectors `x` and `y`,
\( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
\( \sigma_{x} \) and \( \sigma_{y} \) are the standard deviations of `x` and `y`.

Vectors `x` and `y` must be the same size.

### Syntax
`result =` [[fsml(module):fsml_pcc(interface)]]`(x, y)`

### Parameters
`x`: A rank-1 array of type `real`. It must be the same size as `y`.

`y`: A rank-1 array of type `real`. It must be the same size as `x`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x` and `y`.


<br>
# Spearman Rank Correlation Coefficient

## `fsml_scc`

### Description
The procedure computes the Spearman rank correlation coefficient (SCC).
The procedure gets the ranks of cectors `x` and `y`, then
calculates the Pearson correlation coefficient on these ranks.

Vectors `x` and `y` must be the same size.

### Syntax
`result =` [[fsml(module):fsml_scc(interface)]]`(x, y)`

### Parameters
`x`: A rank-1 array of type `real`. It must be the same size as `y`.

`y`: A rank-1 array of type `real`. It must be the same size as `x`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x` and `y`.


<br>
# Examples

```fortran
{!example/modules/example_sts.f90!}
```


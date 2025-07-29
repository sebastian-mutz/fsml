---
title: STS: Sample Statistics
---

# Overview

This is the API documentation for all sample statistics procedures.

[TOC]

<br>
# Mean

## `fsml_mean`
Computes arithmetic mean.
$$ \bar{x} = \frac{1}{n} \cdot \sum_{i=1}^{n} x_i $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \bar{x} \) is the arithmetic mean of `x`.

### Description

### Syntax
`result =` [[fsml(module):fsml_mean(interface)]]`(x)`

### Parameters
`x`: A rank-1 array of type `real`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Variance

## `fsml_var`

### Description
Computes variance.
$$ \operatorname{var}(x) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \bar{x} \) is the arithmetic mean of `x`.

### Syntax
`result =` [[fsml(module):fsml_var(interface)]]`(x)`

### Parameters
`x`: A rank-1 array of type `real`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Standard Deviation

## `fsml_std`

### Description
Computes standard deviation.
$$ \sigma = \sqrt{\operatorname{var}(x)} $$
where \( \operatorname{var}(x) \) is the variance of vector `x`.

### Syntax
`result =` [[fsml(module):fsml_std(interface)]]`(x)`

### Parameters
`x`: A rank-1 array of type `real`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


<br>
# Covariance

## `fsml_cov`

### Description
Computes covariance.
$$ \operatorname{cov}(x, y) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
where \( n \) is the size of (or number of observations in) vectors `x` and `y`,
\( x_i \) and \( y_i \) are individual elements in `x` and `y`, and
\( \bar{x} \) and \( \bar{y} \) are the arithmetic means of `x` and `y`.

Vectors `x` and `y` must be the same size.

### Syntax
`result =` [[fsml(module):fsml_cov(interface)]]`(x, y)`

### Parameters
`x`: A rank-1 array of type `real`. It must be the same size as `y`.

`y`: A rank-1 array of type `real`. It must be the same size as `x`.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x` and `y`.


<br>
# Linear Trend

## `fsml_trend`

### Description
Computes regression coefficient/trend.
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
Computes Pearson correlation coefficient (PCC).
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


# Examples


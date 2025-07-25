---
title: STS: Sample Statistics
---

# <span style="color:#734f96">Overview</span>

This is the API documentation for all sample statistics procedures.

[TOC]

<br>
# <span style="color:#734f96">Mean</span>

## `fsml_mean`
Computes arithmetic mean.
$$ \bar{x} = \frac{1}{n} \cdot \sum_{i=1}^{n} x_i $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \bar{x} \) is the arithmetic mean of `x`.

### Description

### Syntax
`result =` [[fsml(module):fsml_mean(interface)]] `(x)`

### Parameters
`x`: A rank-1 array of type `real`.

### Returns
The result is a scalar and the same type as `x`.


<br>
# <span style="color:#734f96">Variance</span>

## `fsml_var`

### Description
Computes variance.
$$ \operatorname{var}(x) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
where \( n \) is the size of (or number of observations in) vector `x`,
\( x_i \) are individual elements in `x`, and
\( \bar{x} \) is the arithmetic mean of `x`.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Standard Deviation</span>

## `fsml_std`

### Description
Computes standard deviation.
$$ \sigma = \sqrt{\operatorname{var}(x)} $$
where \( \operatorname{var}(x) \) is the variance of vector `x`.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Covariance</span>

## `fsml_cov`

### Description
Computes covariance.
$$ \operatorname{cov}(x, y) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
where \( n \) is the size of (or number of observations in) vectors `x` and `y`,
\( x_i \) and \( y_i \) are individual elements in `x` and `y`, and
\( \bar{x} \) and \( \bar{y} \) are the arithmetic means of `x` and `y`.

Vectors `x` and `y` must be the same size.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Linear Trend</span>

## `fsml_trend`

### Description
Computes regression coefficient/trend.
$$ m = \frac{\operatorname{cov}(x, y)}{\operatorname{var}(x)} $$
where \( m \) is the slope of the regression line (linear trend),
\( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
\( \operatorname{var}(x) \) is the variance of `x`.

Vectors `x` and `y` must be the same size.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Pearson Correlation Coefficient</span>

## `fsml_pcc`

### Description
Computes Pearson correlation coefficient (PCC).
$$ \rho_{x,y} = \frac{\operatorname{cov}(x, y)}{\sigma_x \cdot \sigma_y} $$
where \( \rho_{x,y} \) is the Pearson correlation coefficient for vectors `x` and `y`,
\( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
\( \sigma_{x} \) and \( \sigma_{y} \) are the standard deviations of `x` and `y`.

Vectors `x` and `y` must be the same size.

### Syntax

### Parameters

### Returns



# <span style="color:#734f96">Examples</span>

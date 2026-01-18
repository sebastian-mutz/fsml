---
title: DST: Statistical Distributions
---

# Overview

This is the API documentation for all statistical distribution procedures.
There is a probability density function (PDF), cumulative distribution function (CDF),
and quantile function or percent point function (PPF) for each distribution.

[TOC]

<br>
# Normal Distribution


<br>
## `fsml_norm_pdf`

### Description
Probability density function for the normal distribution.
$$ f(x) = \frac{1}{\sigma \cdot \sqrt{2 \cdot \pi}} e^{ -\frac{1}{2} \cdot \left( \frac{x - \mu}{\sigma} \right)^2 } $$

The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_norm_pdf(interface)]]`(x [, mu, sigma])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.



<br>
## `fsml_norm_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for normal distribution.

The location parameter \( \mu \) (`mu`), scale parameter \( \sigma \) (`sigma`), and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_norm_cdf(interface)]]`(x [, mu, sigma, tail])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_norm_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for normal distribution.

It computes the position (`x`) based on the probability (`p`).
The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

**Note:** The procedure uses bisection method.
Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
will return large negative or positive numbers (highly dependent on the tolerance threshold).

### Syntax
`result =` [[fsml(module):fsml_norm_ppf(interface)]]`(p [, mu, sigma])`

### Parameters
`p`: A scalar of type `real`. It must be between *0.0* and *1.0*.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Student's t Distribution


<br>
## `fsml_t_pdf`

### Description
Probability density function for student t distribution.
Uses intrinsic gamma function (Fortran 2008 and later).
$$ f(x) = \frac{\Gamma\left(\frac{\nu + 1}{2}\right)}{\sqrt{\nu \cdot \pi}\, \cdot \Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu + 1}{2}} $$
where  \(v\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.

The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_t_pdf(interface)]]`(x, df [, mu, sigma])`

### Parameters
`x`: A scalar of type `real`.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_t_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for student t distribution.

The location parameter \( \mu \) (`mu`), scale parameter \( \sigma \) (`sigma`), and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_t_cdf(interface)]]`(x, df [, mu, sigma, tail])`

### Parameters
`x`: A scalar of type `real`.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_t_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for t distribution.

It computes the position (`x`) based on the probability (`p`) and degrees of freedom (`df`).
The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

**Note:** The procedure uses bisection method.
Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
will return large negative or positive numbers (highly dependent on the tolerance threshold).

### Syntax
`result =` [[fsml(module):fsml_t_ppf(interface)]]`(p, df [, mu, sigma])`

### Parameters
`p`: A scalar of type `real`. It must be between *0.0* and *1.0*.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Gamma Distribution


<br>
## `fsml_gamma_pdf`

### Description
Probability density function for gamma distribution.
Uses intrinsic exp function.
$$ f(x) = \frac{\lambda^\alpha}{\Gamma(\alpha)} \cdot x^{\alpha - 1} \cdot e^{-\lambda \cdot x}, \quad x > 0, \ \alpha > 0, \ \lambda > 0 $$

The equation can also be expressed with the scale parameter \( \beta \), which is the inverse of the rate parameter \( \lambda \), so that \( \beta =  \frac{1}{\lambda} \).

$$ f(x) = \frac{1}{\Gamma(\alpha) \, \beta^\alpha} \cdot x^{\alpha - 1} \cdot e^{-x/\beta}, \quad x > 0,\ \alpha > 0,\ \beta > 0 $$

The scale parameters \( \alpha \) (`alpha`) and \( \beta \) (`beta`) and the location parameter (`loc`) are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.


### Syntax
`result =` [[fsml(module):fsml_gamma_pdf(interface)]]`(x [, alpha, beta, loc])`

### Parameters
`x`: A scalar of type `real`.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_gamma_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for gamma distribution.

The scale parameters \( \alpha \) (`alpha`) and \( \beta \) (`beta`), the location parameter (`loc`), and tail option (`tail`) are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.

### Syntax
`result =` [[fsml(module):fsml_gamma_cdf(interface)]]`(x [, alpha, beta, loc, tail])`

### Parameters
`x`: A scalar of type `real`.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_gamma_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for gamma distribution.

It computes the position (`x`) based on the probability (`p`).
The scale parameters \( \alpha \) (`alpha`) and \( \beta \) (`beta`) and the location parameter (`loc`)  are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.

**Note:** The procedure uses bisection method.

### Syntax
`result =` [[fsml(module):fsml_gamma_ppf(interface)]]`(p [, alpha, beta, loc])`

### Parameters
`p`: A scalar of type `real`. It must be between *0.0* and *1.0*.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Exponential Distribution


<br>
## `fsml_exp_pdf`

### Description
Probability density function for exponential distribution.
Uses intrinsic exp function.
$$ f(x) = \lambda \cdot e^{-\lambda \cdot x}, \quad x \geq 0, \ \lambda > 0 $$

The rate parameter \( \lambda \) (`lambda`) and location parameter (`loc`) are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.

### Syntax
`result =` [[fsml(module):fsml_exp_pdf(interface)]]`(x [, lambda, loc])`

### Parameters
`x`: A scalar of type `real`.

`lambda`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_exp_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for exponential distribution.

The rate parameter \( \lambda \) (`lambda`), location parameter (`loc`), and tail option (`tail`) are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.

### Syntax
`result =` [[fsml(module):fsml_exp_cdf(interface)]]`(x [, lambda, loc, tail])`

### Parameters
`x`: A scalar of type `real`.

`lambda`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_exp_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for exponential distribution.
Procedure uses bisection method. `p` should be between 0.0 and 1.0.

It computes the position (`x`) based on the probability (`p`).
The rate parameter \( \lambda \) (`lambda`) and location parameter (`loc`) are optional arguments. The location parameter will shift the distribution in the manner that \( \mu \) does for the normal distribution.

### Syntax
`result =` [[fsml(module):fsml_exp_ppf(interface)]]`(p [, lambda, loc])`

### Parameters
`p`: A scalar of type `real`.

`lambda`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`:  A scalar of type `real`.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Chi-Squared Distribution


<br>
## `fsml_chi2_pdf`

### Description
Probability density function for the chi-squared distribution.
Uses intrinsic exp and gamma function.
$$ f(x) = \frac{x^{\frac{k}{2} - 1} \cdot e^{ - \frac{x}{2} }}{2^{\frac{k}{2}} \cdot \Gamma\left(\frac{k}{2}\right)}, \quad x \geq 0, \ k > 0 $$
where \(k\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.

The procedure calculates the probability based on the provided parameters `x` and `df`.
The location (`loc`) and scale parameter (`scale`) are optional arguments. The parameters `loc` and `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_chi2_pdf(interface)]]`(x, df [, loc, scale])`

### Parameters
`x`: A scalar of type `real`.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_chi2_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for the chi-squared distribution.

The procedure calculates the probability based on the provided parameters `x` and `df`.
The location (`loc`), scale parameter (`scale`), and tail option (`tail`) arguments are optional arguments. The parameters `loc` and `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_chi2_cdf(interface)]]`(x, df [, loc, scale, tail])`

### Parameters
`x`: A scalar of type `real`.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_chi2_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for the chi-squared distribution.
Uses the bisection method for numerical inversion of the CDF.

It computes the position (`x`) based on the probability (`p`) and degrees of freedom (`df`).
The location (`loc`) and scale parameter (`scale`) are optional arguments. The parameters `loc` and  `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_chi2_ppf(interface)]]`(p, df [, loc, scale])`

### Parameters
`p`: A scalar of type `real`.

`df`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# F Distribution

<br>
## `fsml_f_pdf`

### Description
Probability density function for the F distribution.
$$ f(x) = \frac{1}{\mathrm{B}\left(\frac{d_1}{2}, \frac{d_2}{2}\right)} \cdot \left( \frac{d_1}{d_2} \right)^{ \frac{d_1}{2} } \cdot \frac{x^{ \frac{d_1}{2} - 1}}{\left(1 + \frac{d_1}{d_2} x \right)^{ \frac{d_1 + d_2}{2} }}, \quad x > 0 $$
where \(d_1\) = numerator degrees of freedom, \(d_2\) = denominator degrees of freedom and \( B \) is the complete beta function.
(Uses intrinsic gamma function for beta.)

The F distribution is the distribution of \( X = \frac{U_1/d_1}{U_2/d_2} \), where \( U_1 \) and \( U_2 \) are are random variables with chi-square distributions with \( d_1 \) (`d1`) and \( d_2 \) (`d1`) degrees of freedom, respectively.

The procedure calculates the probability based on the provided parameters `x`, `d1`, and `d2`.
The location (`loc`) and scale parameter (`scale`) are optional arguments. The parameters `loc` and `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_f_pdf(interface)]]`(x, d1, d2 [, loc, scale])`

### Parameters
`x`: A scalar of type `real`.

`d1`: A scalar of type `real`. It must be *1.0* or higher.

`d2`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.

<br>
## `fsml_f_cdf`

### Description
Cumulative density function \(F(x) = \mathbb{P}(X \leq x)\) for the F distribution.

The procedure calculates the probability based on the provided parameters `x`, `d1`, and `d2`.
The location (`loc`), scale parameter (`scale`), and tail option (`tail`) arguments are optional arguments. The parameters `loc` and `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_f_cdf(interface)]]`(x, d1, d2 [, loc, scale, tail])`

### Parameters

`x`: A scalar of type `real`.

`d1`: A scalar of type `real`. It must be *1.0* or higher.

`d2`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_f_ppf`

### Description
Percent point function / quantile function \( Q(p) = F^{-1}(p) \) for the F distribution.
Uses the bisection method to numerically invert the CDF.

It computes the position (`x`) based on the probability (`p`) and degrees of freedom (`d1` and `d2`).
The location (`loc`) and scale parameter (`scale`) are optional arguments. The parameters `loc` and  `scale` shift and scale the distribution, respectively, by a specified factor.

### Syntax
`result =` [[fsml(module):fsml_f_ppf(interface)]]`(p, d1, d2 [, loc, scale])`

### Parameters
`p`: A scalar of type `real`.

`d1`: A scalar of type `real`. It must be *1.0* or higher.

`d2`: A scalar of type `real`. It must be *1.0* or higher.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.

<br>
# Generalised Pareto Distribution


<br>
## `fsml_gpd_pdf`

### Description
Probability density function for generalised pareto distribution.
$$ f(x) = \frac{1}{\sigma} \cdot \left ( 1 + \frac{\xi \cdot (x - \mu)}{\sigma} \right)^{-\frac{1}{\xi} - 1}, \quad x \geq \mu, \ \sigma > 0, \ \xi \in \mathbb{R} $$
where \(\xi\)(`xi`) is a shape parameter, \(\mu\) (`mu`) is the location, and \(\sigma\) (`sigma`) is the scale parameter.

The procedure calculates the probability based on the provided parameters `x` and `xi`.
The location (`mu`) and scale parameter (`sigma`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_gpd_pdf(interface)]]`(x, xi [, mu, sigma])`

### Parameters
`x`: A scalar of type `real`.

`xi`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_gpd_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for generalised pareto distribution.

The procedure calculates the probability based on the provided parameters `x` and `xi`.
The location (`mu`) and scale parameter (`sigma`), and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_gpd_cdf(interface)]]`(x, xi [, mu, sigma, tail])`

### Parameters
`x`: A scalar of type `real`.

`xi`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_gpd_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for generalised pareto distribution.
Procedure uses bisection method. `p` must be between 0.0 and 1.0.

It computes the position (`x`) based on the probability (`p`) and `xi`.
The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_gpd_ppf(interface)]]`(p, xi [, mu, sigma])`

### Parameters
`p`: A scalar of type `real`.

`xi`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Logistic Distribution


<br>
## `fsml_logistic_pdf`

### Description
Probability density function for logistic distribution.
$$ f(x) = \frac{e^{-(x-\mu)/s}}{s\left(1 + e^{-(x-\mu)/s}\right)^2} $$
where \(\mu\) (`mu`) is the location and mean, and \(s\) (`scale`) is the scale parameter.

The procedure calculates the probability based on the provided parameter `x`.
The location (`mu`) and scale parameter (`scale`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_logistic_pdf(interface)]]`(x, [, mu, scale])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_logistic_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for logistic distribution.

The procedure calculates the probability based on the provided parameter `x`.
The location (`mu`) and scale parameter (`scale`), and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_logistic_cdf(interface)]]`(x, [, mu, scale, tail])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_logistic_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for logistic distribution.

It computes the position (`x`) based on the probability (`p`). `p` must be between 0.0 and 1.0.
The location parameter (`mu`) and scale parameter (`scale`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_logistic_ppf(interface)]]`(p, [, mu, scale])`

### Parameters
`p`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`scale`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Log-Logistic Distribution


<br>
## `fsml_llogistic_pdf`

### Description
Probability density function for log-logistic distribution.
$$ f(x) = \frac{\frac{\beta}{\alpha}\left(\frac{x}{\alpha}\right)^{\beta-1}}{\left(1 + \left(\frac{x}{\alpha}\right)^{\beta}\right)^2} $$
where \(\alpha\) (`alpha`) is the scale parameter, and \(\beta\) (`beta`) is the shape parameter.


The procedure calculates the probability based on the provided parameter `x`.
The scale (`alpha`), shape (`beta`) and location (`loc`) parameters are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_llogistic_pdf(interface)]]`(x, [, alpha, beta, loc])`

### Parameters
`x`: A scalar of type `real`.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_llogistic_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for logistic distribution.

The procedure calculates the probability based on the provided parameter `x`.
The scale (`alpha`), shape (`beta`) and location (`loc`) parameters and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_llogistic_cdf(interface)]]`(x, [, alpha, beta, loc, tail])`

### Parameters
`x`: A scalar of type `real`.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`tail`: An optional argument and `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `x`.


<br>
## `fsml_llogistic_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for logistic distribution.

It computes the position (`x`) based on the probability (`p`). `p` must be between 0.0 and 1.0.
The scale (`alpha`), shape (`beta`) and location (`loc`) parameters are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_llogistic_ppf(interface)]]`(p, [, alpha, beta, loc])`

### Parameters
`p`: A scalar of type `real`.

`alpha`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`beta`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`loc`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

Invalid argument values will result in the return of a sentinel value (NaN).

### Returns
The result is a scalar of the same type as `p`.


<br>
# Examples

```fortran
{!example/example_dst.f90!}
```

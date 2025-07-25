---
title: DST: Statistical Distributions
---

# <span style="color:#734f96">Overview</span>

This is the API documentation for all statistical distribution procedures.
There is a probability density function (PDF), cumulative distribution function (CDF),
and quantile function or percent point function (PPF) for each distribution.

[TOC]

<br>
# <span style="color:#734f96">Normal Distribution</span>


## `fsml_norm_pdf`

### Description
Probability density function for the normal distribution.
$$ f(x) = \frac{1}{\sigma \cdot \sqrt{2 \cdot \pi}} e^{ -\frac{1}{2} \cdot \left( \frac{x - \mu}{\sigma} \right)^2 } $$

The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_norm_pdf(interface)]] `(x, [,mu, sigma])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


## `fsml_norm_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for normal distribution.

The location parameter \( \mu \) (`mu`), scale parameter \( \sigma \) (`sigma`), and tail option (`tail`) are optional arguments.

### Syntax
`result =` [[fsml(module):fsml_norm_cdf(interface)]] `(x, [,mu, sigma, tail])`

### Parameters
`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

`tail`: An optional argument and positive `character` string. If passed, it must be one of the following: *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to *"left"*.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `x`.


## `fsml_norm_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for normal distribution.

It computes the position (`x`) based on the probability (`p`).
The location parameter \( \mu \) (`mu`) and scale parameter \( \sigma \) (`sigma`) are optional arguments.

**Note:** The procedure uses bisection method.
Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
will return large negative or positive numbers (highly dependent on the tolerance threshold).

### Syntax
`result =` [[fsml(module):fsml_norm_ppf(interface)]] `(p, [,mu, sigma])`

### Parameters

`p`: A scalar of type `real`. It must be between *0.0* and *1.0*.

`mu`: An optional argument and scalar of type `real`. It will default to *0.0* if not passed.

`sigma`: An optional argument and positive scalar of type `real`. If passed, it must be non-zero positive. It will default to *1.0* if not passed.

Invalid argument values will result in the return of a sentinel value.

### Returns
The result is a scalar and the same type as `p`.


<br>
# <span style="color:#734f96">Student's t Distribution</span>

## `fsml_t_pdf`

### Description
Probability density function for student t distribution.
Uses intrinsic gamma function (Fortran 2008 and later).
$$ f(x) = \frac{\Gamma\left(\frac{\nu + 1}{2}\right)}{\sqrt{\nu \cdot \pi}\, \cdot \Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu + 1}{2}} $$
where  \(v\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.

The value for degrees of freedom (`df`) must be 1.0 or higher.
The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.

### Syntax

### Parameters

### Returns

## `fsml_t_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for student t distribution.

The value for degrees of freedom (`df`) must be 1.0 or higher.
The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.
The tail option (`tail`) is an optional argument. If passed, it must be one of the following:
*"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to "left".

### Syntax

### Parameters

### Returns

## `fsml_t_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for t distribution.

Procedure uses bisection method.
Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
will return large negative or positive numbers (highly dependent on the tolerance threshold).

The value for degrees of freedom (`df`) must be 1.0 or higher.
The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Gamma Distribution</span>

## `fsml_gamma_pdf`

### Description
Probability density function for gamma distribution.
Uses intrinsic exp function.
$$ f(x) = \frac{\lambda^\alpha}{\Gamma(\alpha)} \cdot x^{\alpha - 1} \cdot e^{-\lambda \cdot x}, \quad x > 0, \ \alpha > 0, \ \lambda > 0 $$

### Syntax

### Parameters

### Returns

## `fsml_gamma_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for gamma distribution.

### Syntax

### Parameters

### Returns

## `fsml_gamma_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for gamma distribution.
Procedure uses bisection method. `p` should be between 0.0 and 1.0.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Exponential Distribution</span>

## `fsml_exp_pdf`

### Description
Probability density function for exponential distribution.
Uses intrinsic exp function.
$$ f(x) = \lambda \cdot e^{-\lambda \cdot x}, \quad x \geq 0, \ \lambda > 0 $$

### Syntax

### Parameters

### Returns

## `fsml_exp_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for exponential distribution.

### Syntax

### Parameters

### Returns

## `fsml_exp_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for exponential distribution.
Procedure uses bisection method. `p` should be between 0.0 and 1.0.

### Syntax

### Parameters

### Returns

<br>
# <span style="color:#734f96">Chi-Squared Distribution</span>

## `fsml_chi2_pdf`

### Description
Probability density function for the chi-squared distribution.
Uses intrinsic exp and gamma function.
$$ f(x) = \frac{x^{\frac{k}{2} - 1} \cdot e^{ - \frac{x}{2} }}{2^{\frac{k}{2}} \cdot \Gamma\left(\frac{k}{2}\right)}, \quad x \geq 0, \ k > 0 $$
where \(k\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.

### Syntax

### Parameters

### Returns

## `fsml_chi2_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for the chi-squared distribution.

### Syntax

### Parameters

### Returns

## `fsml_chi2_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for the chi-squared distribution.
Uses the bisection method for numerical inversion of the CDF.

### Syntax

### Parameters

### Returns

<br>
# <span style="color:#734f96">F Distribution</span>

## `fsml_f_pdf`

### Description
Probability density function for the F distribution.
$$ f(x) = \frac{1}{\mathrm{B}\left(\frac{d_1}{2}, \frac{d_2}{2}\right)} \cdot \left( \frac{d_1}{d_2} \right)^{ \frac{d_1}{2} } \cdot \frac{x^{ \frac{d_1}{2} - 1}}{\left(1 + \frac{d_1}{d_2} x \right)^{ \frac{d_1 + d_2}{2} }}, \quad x > 0 $$
where \(d_1\) = numerator degrees of freedom, \(d_2\) = denominator degrees of freedom and \( B \) is the complete beta function.
(Uses intrinsic gamma function for beta.)

The F distribution is the distribution of \( X = \frac{U_1/d_1}{U_2/d_2} \), where \( U_1 \) and \( U_2 \) are are random variables with chi-square distributions with \( d_1 \) and \( d_2 \) degrees of freedom, respectively.

### Syntax

### Parameters

### Returns

## `fsml_f_cdf`

### Description
Cumulative density function \(F(x) = \mathbb{P}(X \leq x)\) for the F distribution.

### Syntax

### Parameters

### Returns

## `fsml_f_ppf`

### Description
Percent point function / quantile function \( Q(p) = F^{-1}(p) \) for the F distribution.
Uses the bisection method to numerically invert the CDF.

### Syntax

### Parameters

### Returns


<br>
# <span style="color:#734f96">Generalised Pareto Distribution</span>

## `fsml_gpd_pdf`

### Description
Probability density function for generalised pareto distribution.
$$ f(x) = \frac{1}{\sigma} \cdot \left ( 1 + \frac{\xi \cdot (x - \mu)}{\sigma} \right)^{-\frac{1}{\xi} - 1}, \quad x \geq \mu, \ \sigma > 0, \ \xi \in \mathbb{R} $$
where \(\xi\) is a shape parameter (xi), \(\sigma\) is the scale parameter (sigma), \(\mu\) (mu) is the location (not mean).

### Syntax

### Parameters

### Returns

## `fsml_gpd_cdf`

### Description
Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for generalised pareto distribution.

### Syntax

### Parameters

### Returns

## `fsml_gpd_ppf`

### Description
Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for generalised pareto distribution.
Procedure uses bisection method. `p` must be between 0.0 and 1.0.

### Syntax

### Parameters

### Returns

# <span style="color:#734f96">Examples</span>

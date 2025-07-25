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

<br>
## PDF: `fsml_norm_pdf`
<br>

### Description

Probability density function for the normal distribution.
$$ f(x) = \frac{1}{\sigma \cdot \sqrt{2 \cdot \pi}} e^{ -\frac{1}{2} \cdot \left( \frac{x - \mu}{\sigma} \right)^2 } $$

The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.

### Syntax

`result =` [[fsml(module):fsml_norm_pdf(interface)]] `(x, [,mu, sigma])`


### Input

`x`: A scalar of type `real`.

`mu`: An optional argument and scalar of type `real`.

`sigma`: An optional argument and positive scalar of type `real`.

### Output

The result is a scalar and the same type as `x`.


<br>
## CDF: `fsml_norm_cdf`
<br>

### Description

Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for normal distribution.

The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.
The tail option (`tail`) is an optional argument. If passed, it must be one of the following:
*"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to "left".


<br>
## PPF: `fsml_norm_ppf`
<br>

### Description

Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for normal distribution.

The probability (`p`)must be between 0.0 and 1.0.
The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
It will default to 1.0 if not passed.

The procedure uses bisection method.
Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
will return large negative or positive numbers (highly dependent on the tolerance threshold).

<br>
# <span style="color:#734f96">Student's t Distribution</span>

<br>
# <span style="color:#734f96">Examples</span>

---
title: TST: Statistical Tests
---

# Overview

This is the API documentation for statistical hypothesis tests.

[TOC]

<br>
# Student t-test (1 sample)

## `fsml_ttest_1sample`

### Description
The 1-sample t-test determines if the sample mean has the value specified in the null hypothesis.

The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
\( H_{0} \): \( \bar{x}  =  \mu_0 \), and \( H_{1} \): \( \bar{x} \neq \mu_0 \)

The test statstic \( t \) is calculated as follows:
$$ t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$
where \( \bar{x} \) is the sample mean,
\( s \) is the sample standard deviation,
\( n \) is the sample size, and
\( \mu_0 \) is the population mean.

The degrees of freedom \( \nu \) is calculated as follows:
$$ \nu = n -1 $$

### Syntax

### Parameters

### Returns

<br>
# Paired sample t-test

## `fsml_ttest_paired`
The paired sample t-test (or  dependent sample t-test) determines if
the mean difference between two sample sets are zero.
It is mathematically equivalent to the 1-sample t-test conducted
on the difference vector \( d \) with \( \mu_0 = 0 \).

The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
\( H_{0} \): \( \bar{d}  =  0 \), and \( H_{1} \): \( \bar{d} \neq 0 \)

The test statstic \( t \) is calculated as follows:
$$ t = \frac{\bar{d} - 0}{s_d / \sqrt{n}}$$
where \( \bar{d} \) is the mean of the differences between the sample sets,
\( s_d \) is the standard deviation of the differences, and
\( n \) is the number of paired samples.

The degrees of freedom \( \nu \) is calculated as follows:
$$ \nu = n -1 $$

### Description

### Syntax

### Parameters

### Returns

<br>
# Pooled and Welch's t-test (2 sample)

## `fsml_ttest_2sample`

### Description
The 2-sample t-test determines if two population means \( \mu_1 \) and \( \mu_2\) are the same.
The procedure can handle 2-sample t-tests for equal variances and Welch's t-tests for unequal variances.

The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
\( H_{0} \): \( \mu_1 \ = \mu_2 \), and \( H_{1} \): \( \mu_1 \ \neq \mu_2\)

The procedure defaults to Welch's t-test for unequal variances if `eq_var` is not specified.
In this case, the test statstic \( t \) is calculated as follows:
$$t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{\frac{s^2_1}{n_1} + \frac{s^2_2}{n_2} }} $$
where \( \bar{x}_1 \) and \( \bar{x}_2 \) are the sample means
\( s^2_1 \) and \( s^2_2 \) are the sample standard deviations, and
\( n_1 \) and \( n_1 \) are the sample sizes.
The degrees of freedom \( \nu \) is approximated with the Welch–Satterthwaite equation:
$$ \nu = \frac{\left( \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2} \right)^2} {\frac{\left( \frac{s_1^2}{n_1} \right)^2}{n_1 - 1} + \frac{\left( \frac{s_2^2}{n_2} \right)^2}{n_2 - 1}} $$

If variances are assumed to be equal (`eq_var = .true.`),
the procedure conducts a 2 sample t-test for equal variances, using the pooled standard
deviation \( s_p \) to calculate the t-statistic:
$$ t = \frac{\bar{x}_1 - \bar{x}_2}{s_p \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}} $$

In case of assumed equal variances, the degrees of freedom is calculated as follows:
$$ \nu = n_1 + n_2 - 2 $$

### Syntax

### Parameters

### Returns

<br>
# Analysis of variance (one way)

## `fsml_anova_1way`

### Description
The one-way ANOVA (Analysis of Variance) tests whether three or more population means
\( \mu_1, \mu_2, \dots, \mu_k \) are equal.

The null hypothesis \( H_0 \) and alternative hypothesis \( H_1 \) are defined as:
\( H_0 \): \( \mu_1 = \mu_2 = \cdots = \mu_k \), and
\( H_1 \): At least one \( \mu_j \) differs from the others.

The data is passed to the procedure as a rank-2 array `x`, where each column is a group of observations.
The procedure partitions then the total variability in the data (  \( SS_{total} \) ) into
variability between groups ( \( SS_{between} \); variability explained by groups ), and
variability within groups ( \( SS_{within} \); unexplained or residual variability ), so that
$$ SS_{total} =  SS_{between} + SS_{within} $$

The F-statistic (`f`) is the ratio of the mean sum of squares between groups
to the mean sum of squares within groups:
$$ F = \frac{SS_{between} / (k - 1)}{SS_{within} / (n - k)} $$
where \( k \) is the number of groups, \( n \) is the total number of observations,
\( SS_{between} \) is the sum of squares between groups, and
\( SS_{within} \) is the sum of squares within groups.

The degrees of freedom are \( \nu_1 = k - 1 \) between groups (`df_b`)
and \( \nu_2 = n - k \) within groups (`df_w`).

The resulting p-value (`p`) is computed from the F-distribution:
$$ p = P(F_{\nu_1, \nu_2} > F_{observed}) = 1 - \text{CDF}(F_{observed}) $$

It is computed with the elemental procedure `f_dst_f_cdf_core`.

The ANOVA makes the assumptions that
a) the groups are independent,
b) the observations within each group are normally distributed, and
c) The variances within groups are equal.

### Syntax

### Parameters

### Returns

<br>
# Wilcoxon signed-rank (1 sample)

## `fsml_signedrank_1sample`

### Description
The 1-sample Wilcoxon signed rank test is a non-parametric test that
determines if data comes from a symmetric population with centre \( mu_0 \).
It can be regarded as a non-parametric version of the 1-sample t-test.

If the data consists of independent and similarly distributed samples
from distribution \( D \), the null hypothesis \( H_0 \) can be expressed as:

\( D \) is symmetric around \( \mu = \mu_0 \).

The default alternative hypothesis \( H_1 \) is two-sided and also be
set explicitly (`h1 = "two"`). It can be expressed as:

\( D \) is symmetric around \( \mu \neq \mu_0 \)

If the alternative hypothesis is set to "greater than" (`h1 = "gt"`), it is:

\( D \) is symmetric around \( \mu > \mu_0 \)

If the alternative hypothesis is set to "less than" (`h1 = "lt"`), it is:

\( D \) is symmetric around \( \mu < \mu_0 \)

The test statistic \( W \) is the smaller of the sum of positive and
negative signed ranks:
$$ W = \min \left( \sum_{d_i > 0} R_i, \sum_{d_i < 0} R_i \right) $$

The procedure takes into consideration tied ranks.

### Syntax

### Parameters

### Returns

<br>
# Wilcoxon signed-rank (paired)

## `fsml_signedrank_paired`

### Description
The Wilcoxon signed rank test is a non-parametric test that determines
if two related paired samples come from the same distribution.
It can be regarded as a non-parametric version of the paired t-test.

The Wilcoxon signed rank test is mathematically equivalent to the
1-sample Wilcoxon signed rank test conducted on the difference vector
\( d = x_1 - x_2 \) with \( mu_0 \) set to zero. Consequently, the
the null hypothesis \( H_0 \)  can be expressed as:

Samples \( x_1 - x_2 \) are symmetric around \( \mu = 0 \).

The default alternative hypothesis \( H_1 \) is two-sided and also be
set explicitly (`h1 = "two"`). It can be expressed as:

Samples \( x_1 - x_2 \) are symmetric around \( \mu \neq 0 \)

If the alternative hypothesis is set to "greater than" (`h1 = "gt"`), it is:

Samples \( x_1 - x_2 \) are symmetric around \( \mu > 0 \)

If the alternative hypothesis is set to "less than" (`h1 = "lt"`), it is:

Samples \( x_1 - x_2 \) are symmetric around \( \mu < 0 \)

The procedure takes into consideration tied ranks.

### Syntax

### Parameters

### Returns

<br>
# Mann–Whitney U rank-sum (2 sample)

## `fsml_ranksum`

### Description
The ranks sum test (Wilcoxon rank-sum test or Mann–Whitney U test) is a
non-parametric test to determine if two independent samples \( x_1 \) and
\( x_2 \) are have the same distribution. It can be regarded as the non-parametric
equivalent of the 2-sample t-test.

The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
\( H_0 \): the distributions of \( x_1 \) and \( x_2 \) are equal.
\( H_1 \): the distributions of \( x_1 \) and \( x_2 \) are not equal.

The Mann–Whitney U statistic is calculated for each sample as follows:
$$ U_i = R_i - \frac{n_i \cdot (n_i + 1)}{2} $$
where \( R_i \) is the sum of ranks of sample set \( i \)
and \( n_i \) is the sample size of sample set \( i \).
The final U statistic is:
$$ U = \min(U_1, U_2) $$

The procedure takes into consideration tied ranks.

### Syntax

### Parameters

### Returns

<br>
# Kruskall Wallis H test

## `fsml_kruskalwallis`

### Description
The Kruskal-Wallis H-test is used to determine whether samples originate from the same
distribution without assuming normality. It is therefore considered a nonparametric
alternative to the one-way ANOVA (Analysis of Variance).

The null hypothesis \( H_0 \) and alternative hypothesis \( H_1 \) are defined as:
\( H_0 \): The populations have the same distribution (medians are equal), and
\( H_1 \): At least one population differs from the others.

The data is passed to the procedure as a rank-2 array `x`, where each column is a group of observations.
All values are ranked across the entire dataset, with tied values assigned the average rank.

The Kruskal-Wallis H-statistic (`h`) is computed as:
$$
H = \frac{12}{n(n+1)} \sum_{j=1}^{k} \frac{R_j^2}{n_j} - 3(n+1)
$$
where:
- \( n \) is the total number of observations,
- \( k \) is the number of groups,
- \( n_j \) is the number of observations in group \( j \), and
- \( R_j \) is the sum of ranks in group \( j \).

The degrees of freedom are:
$$ \nu = k - 1 $$
and returned as `df`.

The p-value (`p`) is computed from the chi-squared distribution:
$$ p = P(\chi^2_{\nu} > H_{observed}) = 1 - \text{CDF}(H_{observed}) $$

It is computed using the elemental procedure `f_dst_chi2_cdf_core`.

The Kruskal-Wallis test assumes that:
a) all groups are independent,
b) the response variable is ordinal or continuous,
c) the group distributions have the same shape,
and d) observations are independent both within and between groups.

### Syntax

### Parameters

### Returns







# Examples

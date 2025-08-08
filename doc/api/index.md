---
title: API Reference
---

# API Reference

This is a guide to public-facing API (Application Programming Interface) of FSML!

# Structure

The FSML procedures are categorised into several thematic modules:

- `STS`: Basic (population and sample) statistics (e.g., mean, variance)
- `TST`: Statistical hypothesis testing (e.g., ANOVA)
- `LIN`: Statistical procedures relying heavily on linear algebra (e.g., PCA, OLS regression)
- `NLP`: Nonlinear and algorithmic procedures (e.g., k-means clustering)
- `DST`: Statistical distributions (e.g., Student's t distribution) probability density function (PDF), cumulative distribution function (CDF), and percent point function (PPF)

While the public interfaces do not include these as prefixes, the handbook makes use of these categories to give it more structure.
If you are interested in exploring the code, you will recognise these in module and procedure names.

# Coverage

The following procedures are currently covered and have a public-facing interface. The links will take you directly to the documentation for the API.

| Population and Sample Statistics (STS)             |
| -------------------------------------------------- |
| [Mean](./sts.html#fsml_mean)                       |
| [Median](./sts.html#fsml_median)                   |
| [Variance](./sts.html#fsml_var)                    |
| [Standard deviation](./sts.html#fsml_std)          |
| [Covariance](./sts.html#fsml_cov)                  |
| [Linear trend](./sts.html#fsml_trend)              |
| [Correlation (Pearson)](./sts.html#fsml_pcc)       |
| [Correlation (Spearman rank)](./sts.html#fsml_scc) |

<br>

| Statisical Hypothesis Tests (TST)                                     |
| --------------------------------------------------------------------- |
| [Student t-test (1 sample)](./tst.html#fsml_ttest_1sample)            |
| [Paired sample t-test](./tst.html#fsml_ttest_paired)                  |
| [Pooled t-test (2 sample)](./tst.html#fsml_ttest_2sample)             |
| [Welch's t-test (2 sample)](./tst.html#fsml_ttest_2sample)            |
| [Analysis of variance (one way)](./tst.html#fsml_anova_1way)          |
| [Wilcoxon signed-rank (1 sample)](./tst.html#fsml_signedrank_1sample) |
| [Wilcoxon signed-rank (paired)](./tst.html#fsml_signedrank_paired)    |
| [Mann–Whitney U rank-sum (2 sample)](./tst.html#fsml_ranksum)         |
| [Kruskall Wallis H test](./tst.html#fsml_kruskalwallis)               |

<br>

| Linear Procedures (LIN)                                             |
| ------------------------------------------------------------------- |
| [Principal Component Analysis](./lin.html#fsml_pca)                 |
| [Empirical Orthogonal Functions](./lin.html#fsml_eof)               |
| [Linear Discriminant Analysis (2-Class)](./lin.html#fsml_lda_2class)|
| [Ordinary Least Squares Regression](./lin.html#fsml_ols)            |
| [Ridge Regression](./lin.html#fsml_ridge)                           |

<br>

| Nonlinear Procedures (NLP)                      |
| ----------------------------------------------- |
| Hierarchical Clustering (in progress)           |
| K-Means Clustering (in progress)                |

<br>

| Statistical Distribution Functions (DST)          |
| ------------------------------------------------- |
| [Normal PDF](./dst.html#fsml_norm_pdf)            |
| [Normal CDF](./dst.html#fsml_norm_cdf)            |
| [Normal PPF](./dst.html#fsml_norm_ppf)            |
| [Student's t PDF](./dst.html#fsml_t_pdf)          |
| [Student's t CDF](./dst.html#fsml_t_cdf)          |
| [Student's t PPF](./dst.html#fsml_t_ppf)          |
| [Gamma PDF](./dst.html#fsml_gamma_pdf)            |
| [Gamma CDF](./dst.html#fsml_gamma_cdf)            |
| [Gamma PPF](./dst.html#fsml_gamma_ppf)            |
| [Exponential PDF](./dst.html#fsml_exp_pdf)        |
| [Exponential CDF](./dst.html#fsml_exp_cdf)        |
| [Exponential PPF](./dst.html#fsml_exp_ppf)        |
| [Chi-squared PDF](./dst.html#fsml_chi2_pdf)       |
| [Chi-squared CDF](./dst.html#fsml_chi2_cdf)       |
| [Chi-squared PPF](./dst.html#fsml_chi2_ppf)       |
| [F PDF](./dst.html#fsml_f_pdf)                    |
| [F CDF](./dst.html#fsml_f_cdf)                    |
| [F PPF](./dst.html#fsml_f_ppf)                    |
| [Generalised Pareto PDF](./dst.html#fsml_gpd_pdf) |                                                                          |
| [Generalised Pareto CDF](./dst.html#fsml_gpd_cdf) |
| [Generalised Pareto PPF](./dst.html#fsml_gpd_ppf) |

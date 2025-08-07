---
title: API Reference
---

# API Reference

This is a guide to public-facing API (Application Programming Interface) of FSML!

# Structure

The FSML procedures are categorised into several thematic modules:

- `STS`: Basic (sample) statistics (e.g., mean, variance)
- `DST`: Statistical distributions (e.g., Student's t distribution) probability density function (PDF), cumulative distribution function (CDF), and percent point function (PPF)
- `TST`: Statistical hypothesis testing (e.g., ANOVA)
- `LIN`: Statistical procedures relying heavily on linear algebra (e.g., PCA, OLS regression)
- `NLP`: Nonlinear and algorithmic procedures (e.g., k-means clustering)

While the public interfaces do not include these as prefixes, the handbook makes use of these categories to give it more structure.
If you are interested in exploring the code, you will recognise these in module and procedure names.

# Coverage

The following procedures are currently covered and have a public-facing interface. The links will take you directly to the documentation for the API.

| Sample Statistics (STS)                      | Distributions (DST)                               | Tests (TST)                                                              |
| -------------------------------------------- | ------------------------------------------------- | ------------------------------------------------------------------------ |
| [Mean](./sts.html#fsml_mean)                 | [Normal PDF](./dst.html#fsml_norm_pdf)            | [Student t-test (1 sample)](./tst.html#fsml_ttest_1sample)               |
| [Variance](./sts.html#fsml_var)              | [Normal CDF](./dst.html#fsml_norm_cdf)            | [Paired sample t-test](./tst.html#fsml_ttest_paired)                     |
| [Standard deviation](./sts.html#fsml_std)    | [Normal PPF](./dst.html#fsml_norm_ppf)            | [Pooled t-test (2 sample)](./tst.html#fsml_ttest_2sample)                |
| [Covariance](./sts.html#fsml_cov)            | [Student's t PDF](./dst.html#fsml_t_pdf)          | [Welch's t-test (2 sample)](./tst.html#fsml_ttest_2sample)               |
| [Linear trend](./sts.html#fsml_trend)        | [Student's t CDF](./dst.html#fsml_t_cdf)          | [Analysis of variance (one way)](./tst.html#fsml_anova_1way)             |
| [Correlation (Pearson)](./sts.html#fsml_pcc) | [Student's t PPF](./dst.html#fsml_t_ppf)          | [Wilcoxon signed-rank (1 sample)](./tst.html#fsml_signedrank_1sample)    |
|                                              | [Gamma PDF](./dst.html#fsml_gamma_pdf)            | [Wilcoxon signed-rank (paired)](./tst.html#fsml_signedrank_paired)       |
|                                              | [Gamma CDF](./dst.html#fsml_gamma_cdf)            | [Mann–Whitney U rank-sum (2 sample)](./tst.html#fsml_ranksum)            |
|                                              | [Gamma PPF](./dst.html#fsml_gamma_ppf)            | [Kruskall Wallis H test](./tst.html#fsml_kruskalwallis)                  |
|                                              | [Exponential PDF](./dst.html#fsml_exp_pdf)        |                                                                          |
|                                              | [Exponential CDF](./dst.html#fsml_exp_cdf)        |                                                                          |
|                                              | [Exponential PPF](./dst.html#fsml_exp_ppf)        |                                                                          |
|                                              | [Chi-squared PDF](./dst.html#fsml_chi2_pdf)       |                                                                          |
|                                              | [Chi-squared CDF](./dst.html#fsml_chi2_cdf)       |                                                                          |
|                                              | [Chi-squared PPF](./dst.html#fsml_chi2_ppf)       |                                                                          |
|                                              | [F PDF](./dst.html#fsml_f_pdf)                    |                                                                          |
|                                              | [F CDF](./dst.html#fsml_f_cdf)                    |                                                                          |
|                                              | [F PPF](./dst.html#fsml_f_ppf)                    |                                                                          |
|                                              | [Generalised Pareto PDF](./dst.html#fsml_gpd_pdf) |                                                                          |
|                                              | [Generalised Pareto CDF](./dst.html#fsml_gpd_cdf) |                                                                          |
|                                              | [Generalised Pareto PPF](./dst.html#fsml_gpd_ppf) |                                                                          |

<br>

| Linear Procedures (LIN)                                             | Nonlinear Procedures (NLP)                      |
| ------------------------------------------------------------------- | ----------------------------------------------- |
| [Principal Component Analysis](./lin.html#fsml_pca)                 | Hierarchical Clustering (in progress)           |
| [Empirical Orthogonal Functions](./lin.html#fsml_eof)               | K-Means Clustering (in progress)                |
| [Linear Discriminant Analysis (2-Class)](./lin.html#fsml_lda_2class)|                                                 |
| [Multiple Linear Regression (OLS)](./lin.html#fsml_ols)                   |                                                 |

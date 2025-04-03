# <span style="color:#734f96">FSML - Fortran Statistics and Machine Learning Library</span>

[![GitHub](https://img.shields.io/github/license/sebastian-mutz/fsml)](https://github.com/sebastian-mutz/fsml/blob/main/LICENCE)
![0%](https://progress-bar.xyz/0?title=Alpha)

> [!IMPORTANT]
> FSML is in a pre-alpha state, and only suitable for developers at this point.
>

 <!--
@warning
FSML is in a pre-alpha state, and only suitable for developers at this point.
@endwarning
 -->

## <span style="color:#734f96">Description</span>

FSML is a scientific toolkit consisting of common statistical and machine learning procedures, including basic descriptive statistics (mean, variance, correlation, variance, standard deviation), common univariate and multivariate statistical tests (t-test, ANOVA, Kruskal Wallis, KS), classic statistical analysis (principal component analysis, discriminant analysis, multiple regression), and classic parametric and non-parametric machine learning procedures (k-means clustering, hierarchical clustering, random forests, LASSO- and ridge regression in bootstrap and cross-validation settings), as well as Bayesian classifiers.

## <span style="color:#734f96">Development</span>

FSML is an effort to rewrite, re-structure, clean-up, and enhance old Fortran code I've written for my research in the past 15 years, and to bundle and publish it as a well organised and well documented library.

The published research below uses some of the to-be-reworked code and demonstrates some applications of the above-mentioned methods:

- [Mutz and Ehlers (2019)](https://doi.org/10.5194/esurf-7-663-2019) (k-means and hierarchical clustering, and discriminant analysis)
- [Mutz et al. (2015)](https://doi.org/10.1007/s00382-015-2663-5) (multiple regression in cross validation and bootstrap setting, principal component analysis, and Bayesian classifier)

### <span style="color:#734f96">Alpha</span>

I will consider the library to be in "alpha" once FSML covers the functionality needed to reproduce ~80% of all the Fortran-based data analysis I've conducted (and published) in the past ~15 years.

### <span style="color:#734f96">Beta</span>

This stage is reached once:

- FSML is able to reproduce all of the above mentioned type of analyses without issues.
- FSML fully works with GFortran, LFortran, and Flang compilers.
- FSML has proper documentation.

### <span style="color:#734f96">Progress</span>

#### Basic Statistics

Basic Statistics (Descriptive measures for understanding data).

| Basic Statistics (STS) | Covered |
| ---------------------- | ------- |
| Mean                   | ✓       |
| Variance               | ✓       |
| Standard Deviation     | ✓       |
| Covariance             | ✓       |
| Regression             | ✓       |
| Correlation            | ✓       |

#### Hypothesis Testing

Hypothesis Testing (Statistical tests for inference and comparing groups).

| Hypothesis Testing (TST) | Covered |
| ------------------------ | ------- |
| Student T-Test           | -       |
| Kolmogorov Smirnov Test  | -       |
| Kruskall Wallis Test     | -       |
| Analysis of Variance     | -       |

#### Linear Parametric Models (LPM)

 Models that assume a linear relationship between the features and target, and estimate parameters (coefficients). Methods in brackets are optional, new implementations (rather than reworked old code).

| Linear Parametric Models (LPM)| Covered |
| ----------------------------- | ------- |
| OLS Regression                | -       |
| (LASSO Regression)            | -       |
| (Ridge Regression)            | -       |
| Pincipal Component Analysis   | -       |
| Discriminant Analysis (LDA)   | -       |
| Bayesian Classification       | -       |

#### Non-Linear Models (NPM)

Models for clustering and/or capture non-linear relationships, either explicitly or through flexible structures (such as decision trees). Methods in brackets are optional, new implementations (rather than reworked old code).

| Non-Linear Models (NLM)     | Covered |
| --------------------------- | ------- |
| (Random Forests Regression) | -       |
| Hierarchical Clustering     | -       |
| K-Means Clustering          | -       |
| (Simple Neural Networks)    | -       |


## <span style="color:#734f96">Installation</span>

FSML can be installed/compiled with the [fortran package manager (fpm)](https://github.com/fortran-lang/fpm).



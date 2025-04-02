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

- [Mutz and Ehlers (2019) Detection and Explanation of Spatiotemporal Patterns in Late Cenozoic Palaeoclimate Change Relevant to Earth Surface Processes](https://doi.org/10.5194/esurf-7-663-2019) (k-means and hierarchical clustering, and discriminant analysis)
- [Mutz et al. (2015)  Modelling of future mass balance changes of Norwegian glaciers by application of a dynamical-statistical model](https://doi.org/10.1007/s00382-015-2663-5) (multiple regression in cross validation and bootstrap setting, principal component analysis, and Bayesian classifier)

### <span style="color:#734f96">Alpha</span>

I will consider the library to be in "alpha" once FSML covers the functionality needed to reproduce ~80% of all the analyses I've conducted (and published as first author) in the past ~15 years.

### <span style="color:#734f96">Beta</span>

This stage is reached once:

- FSML is able to reproduce all of the above mentioned analyses without issues.
- FSML fully works with GFortran, LFortran, and Flang compilers.
- FSML has proper documentation.

### <span style="color:#734f96">Progress</span>

| Statistics  | Covered |
| ----------- | ------- |
| Mean        | ✓       |
| Variance    | ✓       |
| Std Dev.    | -       |
| Covariance  | ✓       |
| Regression  | ✓       |
| Correlation | ✓       |

| Tests   | Covered |
| ------- | ------- |
| T-Test  | -       |
| KS-Test | -       |
| KW-Test | -       |
| ANOVA   | -       |

| S. Analyses | Covered |
| ----------- | ------- |
| PCA         | -       |
| 2C-LDA      | -       |
| 2C-MDA      | -       |
| Bayes Class.| -       |
| Multireg.   | -       |
| LASSO       | -       |
| Ridge       | -       |

| Machine Learning  | Covered |
| ----------------- | ------- |
| RF                | -       |
| K-Means Cl.       | -       |
| Hierarchical Cl.  | -       |
| CV setting        | -       |
| Bootstrap setting | -       |

## <span style="color:#734f96">Installation</span>

FSML can be installed/compiled with the [fortran package manager (fpm)](https://github.com/fortran-lang/fpm).



---
title: 'FSML – A Modern Fortran Statistics and Machine Learning Library'
tags:
  - statistics
  - machine learning
authors:
  - name: Sebastian Gerhard Mutz
    orcid: 0000-0001-8180-6150
    affiliation: 1
affiliations:
 - index: 1
   name: School of Geographical and Earth Sciences, University of Glasgow
date: 05 September 2025
bibliography: refs.bib
---

# Summary

`FSML` is a modern Fortran statistics and machine learning library suitable for contemporary research problems and teaching. It includes procedures for basic statistics, hypothesis tests, linear and non-linear methods, and statistical distribution functions.

# Statement of Need

The advances in computing technology over the past two decades have expanded the practical scope of statistics and allowed the widespread use of machine learning (ML). This also transformed research practices and enhanced predictive modelling across many disciplines, including Earth sciences [@boateng:2023], operational weather forecasting [@lang:2024], and more.

Fortran is a well-established general purpose programming language that is commonly adopted in science due to its stability, reliability, performance, and array functionality. It is widely used for parallelised high-performance computing and numerical modelling [e.g., @giorgetta:2018]. The same strengths make it suitable for computationally demanding ML procedures and data-driven predictions. However, despite Fortran’s long history in data-driven prediction and ML [e.g., @breiman:2001; @tomasetti:2009; @gutmann:2022], it has not been as widely adopted in these fields as other languages and lacks well documented, accessible toolkits for statistics and classic ML.

Although projects like Neural-Fortran [@curcic:2019], ATHENA [@taylor:2024], and [FStats](https://github.com/jchristopherson/fstats) cover some important procedures for deep-learning and classic statistics, the Fortran statistics and ML ecosystem remains relatively small. This potentially deters from the use of Fortran, which is already perceived as less accessible than other popular languages due to 1) the lack of familiarity with modern Fortran features, which is exacerbated by stagnating adoption of Fortran at universities, and 2) shortcomings that are currently being addressed by Fortran-lang community projects [@kedward:2022].


`FSML` (Fortran Statistics and Machine Learning) purposefully integrates these projects: It uses [stdlib](https://github.com/fortran-lang/stdlib) for linear algebra, leverages [fpm](https://github.com/fortran-lang/fpm) for easier building and distribution, and is developed to support compilation with the interactive [LFortran](https://github.com/lfortran/lfortran) compiler in addition to GFortran. As such, it builds on recent community efforts and addresses two needs:

1. It adds to the modern Fortran statistics and ML software ecosystem; a richer ecosystem makes Fortran a more attractive choice as a well-established, robust, high-performance language.

2. The use of fpm, the support of free open-source compilers, the extensive documentation, and its permissive license (MIT) facilitate its early adoption and integration into various statistics and ML projects by students, early career researchers, and teachers. It can thus promote the adoption of Fortran.


# Software Description

## Scope

FSML consists of a set of accessible and well-documented statistics and ML procedures, suitable for many contemporary research problems and teaching. These procedures are subdivided into five categories:

- DST: Statistical distribution functions (e.g., the probability density, cumulative distribution, and quantile functions of the Student’s t and generalised Pareto distributions).
- STS: Basic statistics for describing and understanding data (e.g., mean, variance, correlation).
- TST: Parametric and non-parametric hypothesis tests (e.g., analysis of variance, Mann–Whitney U).
- LIN: Statistical procedures relying heavily on linear algebra (e.g., principal component analysis, ridge regression, linear discriminant analysis).
- NLP: Non-linear and algorithmic procedures (e.g., k-means clustering).

![FSML has five thematic modules: Basic statistics (STS), hypothesis tests (TST), linear procedures (LIN), non-linear procedures (NLP), and statistical distribution functions (DST). \label{fig:fig1}](figs/modules.png){ width=50% }

FSML has minimal requirements. It uses Fortran 2008 features, Fortran-lang stdlib for linear algebra, and fpm for easy building and distribution.

**Note:** At the time of publication, LFortran does not reliably compile stdlib. Therefore, early users of FSML are advised to use GFortran.

## Documentation

The FSML handbook is hosted on [fsml.mutz.science](https://fsml.mutz.science/) and can be re-generated from its source files. It includes a detailed, example-rich documentation of the covered procedures, as well as installation instructions and information for contributors.

# Examples

The examples below demonstrate the use of FSML interfaces, using double precision (dp):

* statistical distribution functions:
```fortran
       ! exponential distribution probability density function
       ! with x=0.8 and lambda=0.5
       fx = fsml_exp_pdf(0.8_dp, lambda=0.5_dp)
       ! generalised Pareto cumulative distribution function
       ! with modified shape (xi) and location (mu) parameters
       fx = fsml_gpd_cdf(1.9_dp, xi=1.2_dp, mu=0.6_dp)
```

* sample statistics and dependency measures:
```fortran
       ! mean of vector x
       mean = fsml_mean(x)
       ! sample standard deviation of vector x
       std = fsml_std(x, ddf=1.0_dp)
       ! Pearson correlation coefficient for vectors x1 and x2
       pcc = fsml_pcc(x1, x2)
```

* hypothesis tests:
```fortran
       ! two-sample t-test for unequal variances (Welch's t-test);
       ! returns test statistic (t), degrees of freedom (df), and p
       call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.false.)
       ! one-way ANOVA on a rank-2 array (x2d);
       ! returns f-statistic (f), degrees of freedom (df1, df2) and p
       call fsml_anova_1way(x2d, f, df1, df2, p)
```

* multiple linear ridge regression:
```fortran
       ! ridge regression for 100 data points, 5 variables, and lambda=0.2;
       ! returns y intercept (b0), regression coefficients (b), and R^2 (rsq)
       call fsml_ridge(x, y, 100, 5, 0.2_dp, b0, b, rsq)
```

`FSML`'s repository and handbook includes examples for every public interface.

# Past and Ongoing FSML Projects

The FSML procedures for clustering and linear discriminant analysis were reworked from the code used for climate pattern detection and explanation [@mutz:2018; @mutz:2019]. FSML's empirical orthogonal functions and analysis of variance were used in @mutz:2025. FSML's distribution functions are currently used for modelling climate extremes.

# Acknowledgements

I gratefully acknowledge the Fortran-lang community efforts that this project integrates (fpm, stdlib, and LFortran), as well as the always helpful discussions the with the same community on [Fortran Discourse](https://fortran-lang.discourse.group/) and GitHub. I thank the editor (Jack Atkinson) and reviewers (Ivan Pribec and Mikolaj A. Kowalski) for their time and effort. I also extend my gratitude to Herbert Peck.

# References

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

`FSML` is a set of modern Fortran tools for statistics and machine learning that are commonly used for contemporary research problems. The library includes procedures for descriptive statistics, hypothesis tests, linear and non-linear methods, and statistical distribution functions. Its intuitive modular structure and familiar API also make it suitable for teaching and faciliate easy adoption of Fortran in mathematical and empirical disciplines.

# Statement of Need

The advances in computing technology over the past two decades have expanded the practical scope for the application of statistical methods and made the widespread use of computationally intensive machine learning (ML) techniques feasible. This also transformed research practices and enhanced predictive modelling across many disciplines, including physics [ref], geosciences [@boateng2023], operational weather forecasting [@lang2024], and more.

Fortran is a well-established general purpose programming language that is commonly adopted in science due to its stability, reliability, performance, and array functionality. It is widely used for parallelised high-performance computing and numerical modelling [e.g., @giorgetta2018]. The same strenghts make it equally suitable for computationally demanding ML procedures and data-driven predictions. Furthermore, it is significantly more energy-efficient than other high-level programming languages [@pereira2021], which is another factor to consider as the widespread adoption of computationally demanding ML techniques also increases electricity consumption [@jia2024], adds more stress on Earth’s climate and environments, and creates new challenges as a consequence [e.g., @freitag2021,  @dodge2022]. Despite Fortran’s long history in modern, data-driven prediction and ML [e.g., @breiman2001, @tomasetti2009, @gutmann2022], it has not been as widely adopted in the field as other languages and lacks a well documented, accessible toolkit for statistics and classic ML. While projects like Neural-Fortran, ATHENA, and FStats do cover important procedures for deep-learning and classic statistics, the Fortran ecosystem in this discipline remains relatively small. This potentially deters from the use of Fortran, which is – despite its strenghts - already perceived as less accessible than other popular languages. This can be attributed to 1) shortcomings that are currently being addressed by Fortran-lang community projects like stdlib, fpm, LFortran, the https://fortran-lang.org website and discourse, and 2) the lack of familiarity with modern Fortran features [@kedward2022], which is exacerbated by stagnating adoption of Fortran at universities.

`FSML` is an accessible statistics and ML toolkit suitable for many contemporary research problems. It purposefully integrates Fortran-lang community efforts. Specifically, FSML only relies on stdlib’s interfaces for linear algebra, uses fpm for easier building and distribution, and is developed to support compilation with the community-maintained LFortran in addition to GFortran. As such, it builds on recent community efforts and addresses two gaps:

1. It complements projects like FStats and Neural-Fortran and adds to the Fortran statistics and ML software ecosytem. A richer ecosystem makes Fortran a more attractive choice as a robust, high-performance, energy-efficient option.

2. The use of fpm, the support of free open-source compilers, the extensive user-facing documentation (fsml.mutz.science), and its permissive license (MIT) faciliate its early adoption and integration into various statistics and ML projects by (university) students, early career researchers, and teachers. It can thus help counter the stagnating adoption of Fortran and prevent harmful language monoculture.


# Software Description

## Scope

FSML is a scientific toolkit consisting statistical and machine learning procedure that are frequently needed for common research problems (see examples). These procedures are categorised into five thematic modules:

- STS: Basic (sample) statistics for describing and understanding data (e.g., mean, variance, correlation)
- TST: Parametric and non-parametric hypothesis tests (e.g., Mann–Whitney U , analysis of variance)
- LIN: Statistical procedures relying heavily on linear algebra (e.g., principal component analysis, ridge regression, linear discriminant analysis)
- NLP: Non-linear and algorithmic procedures (e.g., k-means clustering)
- DST: Statistical distributions (e.g., Student's t, generalised Pareto); probability density function, cumulative distribution function, and percent point function

![FSML modules. \label{fig:fig1}](figs/modules.png)

## Structure and Paradigms

1. Minimal and modern: FSML has minimal requirements. It uses Fortran 2008+ intrinsics, Fortran-lang stdlib for linear algebra, and fpm for easy and quick building and distribution. This faciliates easier adoption and encourages transition to a modern, community-driven ecosystem.

2. Modular and intuitive: FSML’s structure is modular, but only to the extend that does not compromise readability. Called procedures are not abstract, but correspond to clearly defined, mostly mathematical procedures (e.g., cdf of a distribution, calling covariance function to construct a covariance matrix, etc.).


## Documentation

The FSML handbook that is hosted on [fsml.mutz.science](http://fsml.mutz.science/). It includes detailed, example-rich documentation of covered procedures, as well as installation instructions and information for contributors.

# Example


API dmonstration for statistical distribution functions:

```fortran

  ! exponential pdf (x=0.8, lambda=0.5)
  print*, fsml_exp_pdf(0.8_wp, lambda=0.5_wp)

  ! left-tailed p-value for t distribution with specified degrees of freedom
  print*, fsml_t_ppf(0.9_wp, df=20.0_wp, mu=0.2_wp, sigma=1.2_wp)

  ! genrealised pareto distribution cdf
  print*, fsml_gpd_cdf(1.9_wp, xi=1.2_wp, mu=0.6_wp, sigma=2.2_wp, tail="left")
```

API dmonstration for sample statistics and dependency measures:

```fortran

  ! mean of vector x1
  print*, "mean (x1): ", fsml_mean(x1)

  ! sample standard deviation of vector x1
  print*, "sample standard deviation (x1): ", fsml_std(x1, ddf=1.0_wp)

  ! Pearson correlation coefficient for x1 and x2
  print*, "Pearson R (x1, x2): ", fsml_pcc(x1, x2)

  ! Spearman rank correlation coefficient for x1 and x2
  print*, "Spearman R (x1, x2): ", fsml_scc(x1, x2)
```

API dmonstration for multiple linear regression with regularisation (ridge regression):

```fortran
  ! call with lambda set to 0.2
  call fsml_ridge(x3, y, nd, nv, 0.2_wp, b0, b, rsq, yhat, se, covb)
```

API dmonstration for hypothesis tests:

```fortran
  ! 2-sample t-test for unequal variances (Welch's t-test)
  call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.false., h1="two")

  ! 1-way ANOVA
  call fsml_anova_1way(x2d, t, df1, df2, p)
```


# Past and Ongoing Projects Using FSML

The FSML procedures for clustering and linear discriminant analysis were reworked from the code used for climate pattern detection and explanation [@mutz2018, @mutz2019]. FSML's empirical orthogonal functions and analysis of variance were used in [@mutz2025]. FSML's distribution function are currently used for modelling climate extremes in a parallel GPU computing environment.

# Future Development

The priority for future development is 1) an increase in scope (e.g., more regressors and distance measures), 2) the addition of helper procedures (reading, data transformation), and 3) the creation of more accessible tutorials.

# Acknowledgements

I gratefully acknowledge the Fortran-lang community efforts that this project benefits from (fpm, stdlib, and LFortran), as well as the always helpful discussions the with the same community on Fortran-lang  discourse and GitHub. I also extend my gratitude to Herbert Peck.


# References

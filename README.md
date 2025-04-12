# <span style="color:#734f96">FSML - Fortran Statistics and Machine Learning Library</span>

[![GitHub](https://img.shields.io/github/license/sebastian-mutz/fsml)](https://github.com/sebastian-mutz/fsml/blob/main/LICENCE)
![5%](https://progress-bar.xyz/5?title=Alpha)

> [!IMPORTANT]
> FSML is in a pre-alpha state, and only suitable for developers at this point.
>

 <!--
@warning
FSML is in a pre-alpha state, and only suitable for developers at this point.
@endwarning
 -->

## <span style="color:#734f96">Description</span>

![logo](assets/logo/FSML_small.png)

FSML is a scientific toolkit consisting of common statistical and machine learning procedures, including basic statistics (e.g., mean, variance, correlation), common statistical tests (e.g., t-test, Kolmogorov-Smirnov), linear parametric methods and models (e.g., principal component analysis, discriminant analysis, Bayesian classifier), and non-linear statistical and machine learning procedures (e.g., k-means clustering).

**Key features:**

 - Common statistics and machine learning techniques (as used in modern research).
 - Familiar/intuitive interface (similarities to popular Python or R libs).
 - Compromise between performance and readability (also suitable for demonstration, teaching, and tinkering).
 - Minimal requirements/dependencies (Fortran 2008 or later, and stdlib).


## <span style="color:#734f96">Example</span>

The example below loads data from a CSV file directly into a simple Fortran dataframe using *fsml_readcsv*. The file stores data for different variables in separate columns. *fsml_mean* and *fsml_var* calculate the mean and variance of a passed vector, respectively. *fsml_corr* computes the Pearson correlation coefficient from the vectors of column 1 and 2.

```fortran
program fortran_statistics
  use fsml
  use iso_fortran_env, dp => real64
  implicit none

  type(fsml_typ_df)  :: df
  character(len=128) :: infile

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  ! read data directly into dataframe
  call fsml_readcsv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

  ! mean of first variable (msl - mean sea level pressure)
  print*, "mean: ", fsml_mean(df%data(:,1))

  ! variance of second variable (t2m - 2m air temperature)
  print*, "variance: ", fsml_var(df%data(:,2))

  ! correlation of msl and t2m
  print*, "correlation coefficent: ", fsml_corr(df%data(:,1), df%data(:,2))

  ! normal pdf (x=0.8)
  print*, fsml_pdf_norm(0.8_dp)

  ! left-tailed p-value for normal distribution with specified mean and standard deviation
  print*, fsml_cdf_norm(2.0_dp, mu=0.3_dp, sigma=1.3_dp, tail="left")

  ! left-tailed p-value for t distribution with t=1.5 and 15 degrees of freedom
  print*, fsml_cdf_t(1.5_dp, df=15, tail="left")

end program fortran_statistics
```


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

> [!IMPORTANT]
> Uses double precision (real64) by default, but can be switched project-wide by changing working precision (wp) in the fsml_typ module.
>

#### Basic Statistics

Basic Statistics (descriptive measures for understanding data).

| Basic Statistics (STS) | Covered |
| ---------------------- | ------- |
| Mean                   | ✓       |
| Variance               | ✓       |
| Standard Deviation     | ✓       |
| Covariance             | ✓       |
| Regression             | ✓       |
| Correlation            | ✓       |

| Distribution Functions | Covered |
| ---------------------- | ------- |
| Normal PDF             | ✓       |
| Normal CDF             | ✓       |
| Normal PPF             | -       |
| Student T PDF          | ✓       |
| Student T CDF          | -       |
| Student T PPF          | -       |

#### Hypothesis Testing

Hypothesis Testing (statistical tests for inference and comparing groups).

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
| Multiple OLS Regression       | -       |
| (LASSO Regression)            | -       |
| (Ridge Regression)            | -       |
| Pincipal Component Analysis   | -       |
| Discriminant Analysis (LDA)   | -       |
| Bayesian Classification       | -       |

#### Non-Linear Models (NLM)

Models for clustering and/or capture non-linear relationships, either explicitly or through flexible structures (such as decision trees). Methods in brackets are optional, new implementations (rather than reworked old code).

| Non-Linear Models (NLM)     | Covered |
| --------------------------- | ------- |
| (Random Forests Regression) | -       |
| Hierarchical Clustering     | -       |
| K-Means Clustering          | -       |
| (Simple Neural Networks)    | -       |

#### Additional Functionality

| Additional Functionality    | Covered |
| --------------------------- | ------- |
| Read from CSV file          | ✓       |
| Read from netCDF file       | -       |
| Simple Fortran Dataframe    | ✓       |
| Cross-Validation            | -       |
| Bootstrapping               | -       |

## <span style="color:#734f96">Installation</span>

FSML can be installed/compiled with the [fortran package manager (fpm)](https://github.com/fortran-lang/fpm).



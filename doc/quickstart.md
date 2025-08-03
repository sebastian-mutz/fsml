---
title: Quick Start
---

# About this Guide

This guide lets you jump right into the action with minimal preparation.

@note If you already have the latest `GFortran` or `LFortran` compiler and `FPM`
, you can continue [here](./quickstart.html#get-the-code). @endnote

[TOC]


<br>
# Get the Tools

### Environment manager (recommended)

The use of [conda](https://docs.conda.io) is recommended for managing your
environment. While strictly speaking not essential, using `conda` is convenient,
clean, and a common way to manage coding environments and packages. You can get
the latest versions of all required packages from the `conda-forge`
channel. Once conda is installed, add the channel, then create a new
environment and activate it:

```
conda config --add channels conda-forge
conda create -n fsml_session
conda activate fsml_session
```

### Compilers

`FSML` supports the [GFortran](https://gcc.gnu.org/fortran/) and
[LFortran](https://lfortran.org/) open-source compilers, but it will most
likely also compile with Intel and other compilers (no compiler specific
extensions are used).

**GFortran** is a very mature and well-established compiler. Consequently, if
you're just looking for battle-tested compiler and smooth experience, GFortran
is your choice. You can install GFortran with conda as follows:

```
conda install -c conda-forge gfortran
```
**LFortran** is a new, interactive compiler that is still under heavy
development. If you're looking to use `FSML` in a jupyter notebook and you're
comfortable with the occasional tinkering (and bug reporting to help improve LFortran!),
LFortran is your choice. You can get LFortran with as follows:

```
conda install -c conda-forge lfortran
```
### FPM

The [Fortran Package Manager (FPM)](https://fpm.fortran-lang.org/) lets you
setup and manage your fortran projects easily. It is modelled after Rust's
Cargo, and many modern Fortran projects, like `FSML`, are offered as FPM
packages. The following will install FPM with conda:

```
conda install -c conda-forge fpm
```
### FORD (optional)

[FORD (FORtran Documenter)](https://forddocs.readthedocs.io/en/stable/) is a
tool that lets you document your Fortran code easily. `FSML`'s documentation is
generated with FORD. If you wish to change and document the code, add
documentation, or regenerate documentation, you can get FORD with conda:

```
conda install -c conda-forge ford
```
To generate `FSML`'s documentation, navigate to the folder containing `doc.md` and
simply run:

```
ford doc.md
```

<br>
# Get the Code

### Download and build with FPM

If you use `FPM`, you will not have to download the code manually. Instead, it
can be listed as a dependency and will be downloaded and compiled automatically
(see next section).

### Download manually

If you prefer to compile or use `FSML` any other way, you can directly [download the source code](https://github.com/sebastian-mutz/fsml/archive/refs/heads/master.zip) or visit the [GitHub repository](https://github.com/sebastian-mutz/fsml).


<br>
# Project Setup

This guide assues you already know how to set up a project with `FPM`. See the [FPM tutorial](https://fpm.fortran-lang.org/tutorial/index.html) for more details if you need an introduction to fpm - it's quick and easy to learn!

### Step 1: Create a new Fortran project

Create new project:

```
fpm new fsml_project
```

This will create a git repository with the `FPM` standard layout.

### Step 2: Add FSML to your project

In your fpm project toml file, add the following lines:

```toml
[dependencies]
fsml = { git = "https://github.com/sebastian-mutz/fsml.git" }
```

Include `FSML` in your fortran project just like you would include any other Fortran module:

```Fortran
use fsml
```

Now you have access to all the `FSML` interfaces. See the [API reference](./api/index.html) for a complete overview of interfaces.

### Step 3: Run your project

Run your project by executing the following from your `FPM` project folder:

```
fpm run
```

This will compile fsml alongside your own project and execure the programme! By default, it will use the `GFortran` compiler. If you want to use `LFortran` instead, run:

```
fpm run --compiler lfortran
```

### Watch your precision!

`FSML` uses double precision (real64) by default, but can be switched project-wide by changing the working precision (wp) in the `fsml_ini` module. The data type in your `FSML` application should match this. If you are happy using real64 (recommended), make sure you import it and declare your reals accordingly:

```Fortran
use iso_fortran_env, only: real64
real(real64) :: my_variable
```

<br>
# Examples


### Reading and Basic Statistics

The example below loads data from a CSV file directly into a simple Fortran dataframe using *fsml_read_csv*. The file stores data for different variables in separate columns. *fsml_mean* and *fsml_var* calculate the mean and variance of a passed vector, respectively. *fsml_pcc* computes the Pearson correlation coefficient from the vectors of column 1 and 2.

```fortran
program fsml_statistics
  use fsml
  use iso_fortran_env, dp => real64
  implicit none

  type(fsml_typ_df)  :: df     ! use fsml dataframe
  character(len=128) :: infile

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  call fsml_read_csv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

  ! mean of first variable (msl - mean sea level pressure)
  print*, "mean: ", fsml_mean(df%data(:,1))

  ! variance of second variable (t2m - 2m air temperature)
  print*, "variance: ", fsml_var(df%data(:,2))

  ! correlation of msl and t2m
  print*, "pearson correlation coefficent: ", fsml_pcc(df%data(:,1), df%data(:,2))

end program fsml_statistics
```

### Statistical Distribution Functions

The example below demonstrates the use of pdf (probability density function), cdf (cumulative distribution function), and ppf (percent point function/quantile function) procedures for different statistical distributions.

```fortran
program fsml_distributions
  use fsml
  use iso_fortran_env, dp => real64
  implicit none

  ! exponential pdf (x=0.8)
  print*, fsml_exp_pdf(0.8_dp)

  ! left-tailed p-value for normal distribution with specified mean and standard deviation
  print*, fsml_norm_cdf(2.0_dp, mu=0.3_dp, sigma=1.3_dp, tail="left")

  ! genrealised pareto distribution cdf
  print*, fsml_gpd_cdf(1.9_dp, xi=1.2_dp, mu=0.6_dp, sigma=2.2_dp, tail="left")

  ! chi square distribution ppf
  print*, fsml_chi2_ppf(0.2_dp, df=10, loc=2.0_dp, scale=1.2_dp)

end program fsml_distributions
```

### Statistical Tests

In the example below, we create 2 data vectors and demonstrate how to conduct a range of statistical tests on them. We then let the programme output the test results in a clean format. The commented values are the expected output.


```fortran
program fsml_tests
  use fsml
  use iso_fortran_env, dp => real64
  implicit none

  integer , parameter :: n1 = 10, n2 = 13
  real(dp), parameter :: mu = 1.25_dp
  real(dp), parameter :: x1(n1) = [1.10_dp, 1.02_dp, 1.12_dp, 2.01_dp, 1.92_dp, 1.01_dp, 1.10_dp, 1.26_dp, 1.51_dp, 1.01_dp]
  real(dp), parameter :: x2(n2) = [2.10_dp, 1.02_dp, 1.22_dp, 1.05_dp, 0.95_dp, 1.02_dp, 2.00_dp, 1.05_dp, 1.12_dp, 1.20_dp, 1.12_dp, 1.01_dp, 1.12_dp]
  real(dp)            :: t, p, df

  ! 1-sample t-test
  call fsml_ttest_1sample(x1, mu, t, df, p, h1="two")
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  ! test statistic (t):     0.4672
  ! degrees of freedom:     9.0000
  ! p-value:                0.6514

  ! 2-sample pooled t-test (assume equal variances)
  call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.true., h1="two")
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  ! test statistic (t):     0.4862
  ! degrees of freedom:    21.0000
  ! p-value:                0.6319

  ! 1-sample Wilcoxon signed rank test
  call fsml_signedrank_1sample(x1, mu, t, p, h1="two")
  write(*,'(A)') "> 1-sample Wilcoxon signed rank test"
  write(*,'(A,F10.4)')   "  test statistic (w): ", t
  write(*,'(A,F10.4)')   "  p-value:            ", p
  ! test statistic (w):    27.0000
  ! p-value:                0.9594

  ! 2-sample Wilcoxon Mann-Whitney U rank sum test
  call fsml_ranksum(x1, x2, t, p, h1="two")
  write(*,'(A)') "> Mann-Whitney U rank sum test"
  write(*,'(A,F10.4)')   "  test statistic (U): ", t
  write(*,'(A,F10.4)')   "  p-value:            ", p
  ! test statistic (U):    59.5000
  ! p-value:                0.7330

end program fsml_tests
```

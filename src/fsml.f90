module fsml

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Main module; provides interfaces.                                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! FSML interface module.

  ! load modules
  use :: fsml_ini
  use :: fsml_typ
  use :: fsml_dat
  use :: fsml_sts
  use :: fsml_dst
  use :: fsml_tst
  use :: fsml_utl

  ! basic options
  implicit none
  private

  ! public statistics procedures
  public :: fsml_mean, fsml_var, fsml_std, fsml_cov, fsml_trend, fsml_pcc
  ! public distribution procedures
  public :: fsml_norm_pdf, fsml_norm_cdf, fsml_norm_ppf
  public :: fsml_t_pdf, fsml_t_cdf, fsml_t_ppf
  public :: fsml_gamma_pdf, fsml_gamma_cdf, fsml_gamma_ppf
  public :: fsml_exp_pdf, fsml_exp_cdf, fsml_exp_ppf
  public :: fsml_chi2_pdf, fsml_chi2_cdf, fsml_chi2_ppf
  public :: fsml_f_pdf, fsml_f_cdf, fsml_f_ppf
  public :: fsml_gpd_pdf, fsml_gpd_cdf, fsml_gpd_ppf
  ! public statistical tests
  public :: fsml_ttest_1sample, fsml_ttest_paired, fsml_ttest_2sample
  public :: fsml_signedrank_1sample, fsml_signedrank_paired, fsml_ranksum
  ! public utility procedures
  public :: fsml_rank
  ! public data/io procedures
  public :: fsml_read_csv
  ! public derived types
  public :: fsml_typ_df

! ==== Interfaces
! Interfaces to offer simpler, consistent public procedure names


! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Basic Statistics

! arithmetic mean
interface fsml_mean
  !! Computes arithmetic mean.
  !! $$ \bar{x} = \frac{1}{n} \cdot \sum_{i=1}^{n} x_i $$
  !! where \( n \) is the size of (or number of observations in) vector `x`,
  !! \( x_i \) are individual elements in `x`, and
  !! \( \bar{x} \) is the arithmetic mean of `x`.
  module procedure f_sts_mean
end interface

! variance
interface fsml_var
  !! Computes variance.
  !! $$ \operatorname{var}(x) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
  !! where \( n \) is the size of (or number of observations in) vector `x`,
  !! \( x_i \) are individual elements in `x`, and
  !! \( \bar{x} \) is the arithmetic mean of `x`.
  module procedure f_sts_var
end interface

! standard deviation
interface fsml_std
  !! Computes standard deviation.
  !! $$ \sigma = \sqrt{\operatorname{var}(x)} $$
  !! where \( \operatorname{var}(x) \) is the variance of vector `x`.
  module procedure f_sts_std
end interface

! covariance
interface fsml_cov
  !! Computes covariance.
  !! $$ \operatorname{cov}(x, y) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
  !! where \( n \) is the size of (or number of observations in) vectors `x` and `y`,
  !! \( x_i \) and \( y_i \) are individual elements in `x` and `y`, and
  !! \( \bar{x} \) and \( \bar{y} \) are the arithmetic means of `x` and `y`.
  !!
  !! Vectors `x` and `y` must be the same size.
  module procedure f_sts_cov
end interface

! linear trend (regression coefficient)
interface fsml_trend
  !! Computes regression coefficient/trend.
  !! $$ m = \frac{\operatorname{cov}(x, y)}{\operatorname{var}(x)} $$
  !! where \( m \) is the slope of the regression line (linear trend),
  !! \( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
  !! \( \operatorname{var}(x) \) is the variance of `x`.
  !!
  !! Vectors `x` and `y` must be the same size.
  module procedure f_sts_trend
end interface

! Pearson correlation coefficient
interface fsml_pcc
  !! Computes Pearson correlation coefficient (PCC).
  !! $$ \rho_{x,y} = \frac{\operatorname{cov}(x, y)}{\sigma_x \cdot \sigma_y} $$
  !! where \( \rho_{x,y} \) is the Pearson correlation coefficient for vectors `x` and `y`,
  !! \( \operatorname{cov}(x, y) \) is the covariance of `x` and `y`, and
  !! \( \sigma_{x} \) and \( \sigma_{y} \) are the standard deviations of `x` and `y`.
  !!
  !! Vectors `x` and `y` must be the same size.
  module procedure f_sts_pcc
end interface


! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Statistical Distributions

! normal distribution pdf
interface fsml_norm_pdf
  !! Probability density function for normal distribution.
  !! $$ f(x) = \frac{1}{\sigma \cdot \sqrt{2 \cdot \pi}} e^{ -\frac{1}{2} \cdot \left( \frac{x - \mu}{\sigma} \right)^2 } $$
  !!
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  module procedure f_dst_norm_pdf
end interface

! normal distribution cdf
interface fsml_norm_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for normal distribution.
  !!
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  !! The tail option (`tail`) is an optional argument. If passed, it must be one of the following:
  !! *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to "left".
  module procedure f_dst_norm_cdf
end interface

! normal distribution ppf
interface fsml_norm_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for normal distribution.
  !!
  !! The probability (`p`)must be between 0.0 and 1.0.
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  !!
  !! The procedure uses bisection method.
  !! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
  !! will return large negative or positive numbers (highly dependent on the tolerance threshold).
  module procedure f_dst_norm_ppf
end interface

! t distribution pdf
interface fsml_t_pdf
  !! Probability density function for student t distribution.
  !! Uses intrinsic gamma function (Fortran 2008 and later).
  !! $$ f(x) = \frac{\Gamma\left(\frac{\nu + 1}{2}\right)}{\sqrt{\nu \cdot \pi}\, \cdot \Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu + 1}{2}} $$
  !! where  \(v\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.
  !!
  !! The value for degrees of freedom (`df`) must be 1.0 or higher.
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  module procedure f_dst_t_pdf
end interface

! t distribution cdf
interface fsml_t_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for student t distribution.
  !!
  !! The value for degrees of freedom (`df`) must be 1.0 or higher.
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  !! The tail option (`tail`) is an optional argument. If passed, it must be one of the following:
  !! *"left"*, *"right"*, *"two"*, or *"confidence"*. If not passed, it will default to "left".
  module procedure f_dst_t_cdf
end interface

! t distribution ppf
interface fsml_t_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for t distribution.
  !!
  !! Procedure uses bisection method.
  !! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
  !! will return large negative or positive numbers (highly dependent on the tolerance threshold).
  !!
  !! The value for degrees of freedom (`df`) must be 1.0 or higher.
  !! The location parameter (`mu`) is an optional argument and will default to 0.0 if not passed.
  !! The scale parameter (`sigma`) is an optional argument. If passed, it must be non-zero positive.
  !! It will default to 1.0 if not passed.
  module procedure f_dst_t_ppf
end interface

! gamma distribution pdf
interface fsml_gamma_pdf
  !! Probability density function for gamma distribution.
  !! Uses intrinsic exp function.
  !! $$ f(x) = \frac{\lambda^\alpha}{\Gamma(\alpha)} \cdot x^{\alpha - 1} \cdot e^{-\lambda \cdot x}, \quad x > 0, \ \alpha > 0, \ \lambda > 0 $$
  module procedure f_dst_gamma_pdf
end interface

! gamma distribution cdf
interface fsml_gamma_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for gamma distribution.
  module procedure f_dst_gamma_cdf
end interface

! gamma distribution ppf
interface fsml_gamma_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for gamma distribution.
  !! Procedure uses bisection method. p should be between 0.0 and 1.0.
  module procedure f_dst_gamma_ppf
end interface

! exponential distribution pdf
interface fsml_exp_pdf
  !! Probability density function for exponential distribution.
  !! Uses intrinsic exp function.
  !! $$ f(x) = \lambda \cdot e^{-\lambda \cdot x}, \quad x \geq 0, \ \lambda > 0 $$
  module procedure f_dst_exp_pdf
end interface

! exponential distribution cdf
interface fsml_exp_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for exponential distribution.
  module procedure f_dst_exp_cdf
end interface

! exponential distribution ppf
interface fsml_exp_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for exponential distribution.
  !! Procedure uses bisection method. p should be between 0.0 and 1.0.
  module procedure f_dst_exp_ppf
end interface

! chi-squared distribution pdf
interface fsml_chi2_pdf
  !! Probability density function for the chi-squared distribution.
  !! Uses intrinsic exp and gamma function.
  !! $$ f(x) = \frac{x^{\frac{k}{2} - 1} \cdot e^{ - \frac{x}{2} }}{2^{\frac{k}{2}} \cdot \Gamma\left(\frac{k}{2}\right)}, \quad x \geq 0, \ k > 0 $$
  !! where \(k\) = degrees of freedom (df) and \(\Gamma\) is the gamma function.
  module procedure f_dst_chi2_pdf
end interface

! chi-squared distribution cdf
interface fsml_chi2_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for the chi-squared distribution.
  module procedure f_dst_chi2_cdf
end interface

! chi-squared distribution ppf
interface fsml_chi2_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for the chi-squared distribution.
  !! Uses the bisection method for numerical inversion of the CDF.
  module procedure f_dst_chi2_ppf
end interface

! f distribution pdf
interface fsml_f_pdf
  !! Probability density function for the F distribution.
  !! $$ f(x) = \frac{1}{\mathrm{B}\left(\frac{d_1}{2}, \frac{d_2}{2}\right)} \cdot \left( \frac{d_1}{d_2} \right)^{ \frac{d_1}{2} } \cdot \frac{x^{ \frac{d_1}{2} - 1}}{\left(1 + \frac{d_1}{d_2} x \right)^{ \frac{d_1 + d_2}{2} }}, \quad x > 0 $$
  !! where \(d_1\) = numerator degrees of freedom, \(d_2\) = denominator degrees of freedom and \( B \) is the complete beta function.
  !! (Uses intrinsic gamma function for beta.)
  !!
  !! The F distribution is the distribution of \( X = \frac{U_1/d_1}{U_2/d_2} \), where \( U_1 \) and \( U_2 \) are are random variables with chi-square distributions with \( d_1 \) and \( d_2 \) degrees of freedom, respectively.
  module procedure f_dst_f_pdf
end interface

! f distribution cdf
interface fsml_f_cdf
  !! Cumulative density function \(F(x) = \mathbb{P}(X \leq x)\) for the F distribution.
  module procedure f_dst_f_cdf
end interface

! f distribution ppf
interface fsml_f_ppf
  !! Percent point function / quantile function \( Q(p) = F^{-1}(p) \) for the F distribution.
  !! Uses the bisection method to numerically invert the CDF.
  module procedure f_dst_f_ppf
end interface

! generalised pareto distribution pdf
interface fsml_gpd_pdf
  !! Probability density function for generalised pareto distribution.
  !! $$ f(x) = \frac{1}{\sigma} \cdot \left ( 1 + \frac{\xi \cdot (x - \mu)}{\sigma} \right)^{-\frac{1}{\xi} - 1}, \quad x \geq \mu, \ \sigma > 0, \ \xi \in \mathbb{R} $$
  !! where \(\xi\) is a shape parameter (xi), \(\sigma\) is the scale parameter (sigma), \(\mu\) (mu) is the location (not mean).
  module procedure f_dst_gpd_pdf
end interface

! generalised pareto distribution cdf
interface fsml_gpd_cdf
  !! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for generalised pareto distribution.
  module procedure f_dst_gpd_cdf
end interface

! generalised pareto distribution ppf
interface fsml_gpd_ppf
  !! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for generalised pareto distribution.
  !! Procedure uses bisection method. p must be between 0.0 and 1.0.
  module procedure f_dst_gpd_ppf
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Statistical Tests

! 1 sample t-test
interface fsml_ttest_1sample
  module procedure s_tst_ttest_1s
end interface

! paired sample t-test
interface fsml_ttest_paired
  module procedure s_tst_ttest_paired
end interface

! 2 sample t-test
interface fsml_ttest_2sample
  module procedure s_tst_ttest_2s
end interface

! Wilcoxon 1-sample signed rank test
interface fsml_signedrank_1sample
  module procedure s_tst_signedrank_1s
end interface

! Wilcoxon 2-sample (paired sample) signed rank test
interface fsml_signedrank_paired
  module procedure s_tst_signedrank_2s
end interface

! Wilcoxon rank-sum / Mann–Whitney U test
interface fsml_ranksum
  module procedure s_tst_ranksum
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Utilities

! ranks elements in an array
interface fsml_rank
  module procedure s_utl_rank
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Data Handling

! read csv file into dataframe
interface fsml_read_csv
  module procedure s_dat_read_csv
end interface

end module fsml

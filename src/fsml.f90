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

! TODO: copy module (to fsml_core), but expose pure core routines with same api. Allows devs to use either core routines or safe ones, depending on needs.

! FORD
!! FSML interface module.

  ! load modules
  use :: fsml_ini
  use :: fsml_typ
  use :: fsml_dat
  use :: fsml_sts
  use :: fsml_dst
  use :: fsml_tst
  use :: fsml_lin
  use :: fsml_nlp
  use :: fsml_utl

  ! basic options
  implicit none
  private

  ! public statistics procedures
  public :: fsml_mean, fsml_median
  public :: fsml_var, fsml_std
  public :: fsml_cov, fsml_trend, fsml_pcc, fsml_scc
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
  public :: fsml_anova_1way
  public :: fsml_signedrank_1sample, fsml_signedrank_paired, fsml_ranksum
  public :: fsml_kruskalwallis
  ! public linear (algebra) procedures
  public :: fsml_eof, fsml_pca, fsml_lda_2class, fsml_ols, fsml_ridge, fsml_mahalanobis
  ! public nonlinear procedures
  public :: fsml_hcluster, fsml_kmeans
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

! median
interface fsml_median
  !! Computes median of vector `x`. The procedures can handle tied ranks.
  module procedure f_sts_median
end interface

! variance
interface fsml_var
  !! Computes the population or sample variance (depending on passed arguments).
  !! $$ \operatorname{var}(x) = \frac{1}{n - \nu} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
  !! where \( n \) is the size of (or number of observations in) vector `x`,
  !! \( x_i \) are individual elements in `x`,
  !! \( \nu \) (`ddof`) is a degrees of freedom adjustment
  !! (`ddof = 0.0` for population variance, `ddof = 1.0` for sample variance), and
  !! \( \bar{x} \) is the arithmetic mean of `x`.
  module procedure f_sts_var
end interface

! standard deviation
interface fsml_std
  !! Computes the population or sample standard deviation (depending on passed arguments).
  !! $$ \sigma = \sqrt{\operatorname{var}(x)} $$
  !! where \( \operatorname{var}(x) \) is the variance of vector `x`.
  !! \( \nu \) (`ddof`) can also be passed and serves as a degrees of freedom adjustment
  !! when the variance is caulculated. (`ddof = 0.0` for population standard deviation,
  !! `ddof = 1.0` for sample standard deviation)
  module procedure f_sts_std
end interface

! covariance
interface fsml_cov
  !! Computes the population or sample covariance (depending on passed arguments).
  !! $$ \operatorname{cov}(x, y) = \frac{1}{n - \nu} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
  !! where \( n \) is the size of (or number of observations in) vectors `x` and `y`,
  !! \( x_i \) and \( y_i \) are individual elements in `x` and `y`,
  !! \( \nu \) (`ddof`) is a degrees of freedom adjustment
  !! (`ddof = 0.0` for population variance, `ddof = 1.0` for sample variance), and
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

! Spearman rank correlation coefficient
interface fsml_scc
  !! Computes the Spearman rank correlation coefficient (SCC).
  !! The procedure gets the ranks of cectors `x` and `y`, then
  !! calculates the Pearson correlation coefficient on these ranks.
  !!
  !! Vectors `x` and `y` must be the same size.
  module procedure f_sts_scc
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
  !! Procedure uses bisection method. `p` should be between 0.0 and 1.0.
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
  !! Procedure uses bisection method. `p` should be between 0.0 and 1.0.
  module procedure f_dst_exp_ppf
end interface

! chi-squared distribution pdf
interface fsml_chi2_pdf
  !! Probability density function for the chi-squared distribution.
  !! Uses intrinsic exp and gamma function.
  !! $$ f(x) = \frac{x^{\frac{k}{2} - 1} \cdot e^{ - \frac{x}{2} }}{2^{\frac{k}{2}} \cdot \Gamma\left(\frac{k}{2}\right)}, \quad x \geq 0, \ k > 0 $$
  !! where \(k\) = degrees of freedom (`df`) and \(\Gamma\) is the gamma function.
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
  !! Procedure uses bisection method. `p` must be between 0.0 and 1.0.
  module procedure f_dst_gpd_ppf
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Statistical Tests

! 1 sample t-test
interface fsml_ttest_1sample
  !! The 1-sample t-test determines if the sample mean has the value specified in the null hypothesis.
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
  !! \( H_{0} \): \( \bar{x}  =  \mu_0 \), and \( H_{1} \): \( \bar{x} \neq \mu_0 \)
  !!
  !! ### Procedure:
  !! The test statstic \( t \) is calculated as follows:
  !! $$ t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$
  !! where \( \bar{x} \) is the sample mean,
  !! \( s \) is the sample standard deviation,
  !! \( n \) is the sample size, and
  !! \( \mu_0 \) is the population mean.
  !!
  !! The degrees of freedom \( \nu \) is calculated as follows:
  !! $$ \nu = n -1 $$
  module procedure s_tst_ttest_1s
end interface

! paired sample t-test
interface fsml_ttest_paired
  !! The paired sample t-test (or  dependent sample t-test) determines if
  !! the mean difference between two sample sets are zero.
  !! It is mathematically equivalent to the 1-sample t-test conducted
  !! on the difference vector \( d \) with \( \mu_0 = 0 \).
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
  !! \( H_{0} \): \( \bar{d}  =  0 \), and \( H_{1} \): \( \bar{d} \neq 0 \)
  !!
  !! ### Procedure:
  !! The test statstic \( t \) is calculated as follows:
  !! $$ t = \frac{\bar{d} - 0}{s_d / \sqrt{n}}$$
  !! where \( \bar{d} \) is the mean of the differences between the sample sets,
  !! \( s_d \) is the standard deviation of the differences, and
  !! \( n \) is the number of paired samples.
  !!
  !! The degrees of freedom \( \nu \) is calculated as follows:
  !! $$ \nu = n -1 $$
  module procedure s_tst_ttest_paired
end interface

! 2 sample t-test
interface fsml_ttest_2sample
  !! The 2-sample t-test determines if two population means \( \mu_1 \) and \( \mu_2\) are the same.
  !! The procedure can handle 2-sample t-tests for equal variances and Welch's t-tests for unequal variances.
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
  !! \( H_{0} \): \( \mu_1 \ = \mu_2 \), and \( H_{1} \): \( \mu_1 \ \neq \mu_2\)
  !!
  !! ### Procedure:
  !! The procedure defaults to Welch's t-test for unequal variances if `eq_var` is not specified.
  !! In this case, the test statstic \( t \) is calculated as follows:
  !! $$t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{\frac{s^2_1}{n_1} + \frac{s^2_2}{n_2} }} $$
  !! where \( \bar{x}_1 \) and \( \bar{x}_2 \) are the sample means
  !! \( s^2_1 \) and \( s^2_2 \) are the sample standard deviations, and
  !! \( n_1 \) and \( n_1 \) are the sample sizes.
  !! The degrees of freedom \( \nu \) is approximated with the Welch–Satterthwaite equation:
  !! $$ \nu = \frac{\left( \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2} \right)^2} {\frac{\left( \frac{s_1^2}{n_1} \right)^2}{n_1 - 1} + \frac{\left( \frac{s_2^2}{n_2} \right)^2}{n_2 - 1}} $$
  !!
  !! If variances are assumed to be equal (`eq_var = .true.`),
  !! the procedure conducts a 2 sample t-test for equal variances, using the pooled standard
  !! deviation \( s_p \) to calculate the t-statistic:
  !! $$ t = \frac{\bar{x}_1 - \bar{x}_2}{s_p \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}} $$
  !!
  !! In case of assumed equal variances, the degrees of freedom is calculated as follows:
  !! $$ \nu = n_1 + n_2 - 2 $$
  module procedure s_tst_ttest_2s
end interface

! 1-way ANOVA
interface fsml_anova_1way
  !! The one-way ANOVA (Analysis of Variance) tests whether three or more population means
  !! \( \mu_1, \mu_2, \dots, \mu_k \) are equal.
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_0 \) and alternative hypothesis \( H_1 \) are defined as:
  !! \( H_0 \): \( \mu_1 = \mu_2 = \cdots = \mu_k \), and
  !! \( H_1 \): At least one \( \mu_j \) differs from the others.
  !!
  !! ### Procedure:
  !! The data is passed to the procedure as a rank-2 array `x`, where each column is a group of observations.
  !! The procedure partitions then the total variability in the data (  \( SS_{total} \) ) into
  !! variability between groups ( \( SS_{between} \); variability explained by groups ), and
  !! variability within groups ( \( SS_{within} \); unexplained or residual variability ), so that
  !! $$ SS_{total} =  SS_{between} + SS_{within} $$
  !!
  !! The F-statistic (`f`) is the ratio of the mean sum of squares between groups
  !! to the mean sum of squares within groups:
  !! $$ F = \frac{SS_{between} / (k - 1)}{SS_{within} / (n - k)} $$
  !! where \( k \) is the number of groups, \( n \) is the total number of observations,
  !! \( SS_{between} \) is the sum of squares between groups, and
  !! \( SS_{within} \) is the sum of squares within groups.
  !!
  !! The degrees of freedom are \( \nu_1 = k - 1 \) between groups (`df_b`)
  !! and \( \nu_2 = n - k \) within groups (`df_w`).
  !!
  !! The resulting p-value (`p`) is computed from the F-distribution:
  !! $$ p = P(F_{\nu_1, \nu_2} > F_{observed}) = 1 - \text{CDF}(F_{observed}) $$
  !!
  !! It is computed with the elemental procedure `f_dst_f_cdf_core`.
  !!
  !! The ANOVA makes the assumptions that
  !! a) the groups are independent,
  !! b) the observations within each group are normally distributed, and
  !! c) The variances within groups are equal.
  module procedure s_tst_anova_1w
end interface

! Wilcoxon 1-sample signed rank test
interface fsml_signedrank_1sample
  !! The 1-sample Wilcoxon signed rank test is a non-parametric test that
  !! determines if data comes from a symmetric population with centre \( mu_0 \).
  !! It can be regarded as a non-parametric version of the 1-sample t-test.
  !!
  !! ### Hypotheses:
  !! If the data consists of independent and similarly distributed samples
  !! from distribution \( D \), the null hypothesis \( H_0 \) can be expressed as:
  !!
  !! \( D \) is symmetric around \( \mu = \mu_0 \).
  !!
  !! The default alternative hypothesis \( H_1 \) is two-sided and also be
  !! set explicitly (`h1 = "two"`). It can be expressed as:
  !!
  !! \( D \) is symmetric around \( \mu \neq \mu_0 \)
  !!
  !! If the alternative hypothesis is set to "greater than" (`h1 = "gt"`), it is:
  !!
  !! \( D \) is symmetric around \( \mu > \mu_0 \)
  !!
  !! If the alternative hypothesis is set to "less than" (`h1 = "lt"`), it is:
  !!
  !! \( D \) is symmetric around \( \mu < \mu_0 \)
  !!
  !! ### Procedure:
  !! The test statistic \( W \) is the smaller of the sum of positive and
  !! negative signed ranks:
  !! $$ W = \min \left( \sum_{d_i > 0} R_i, \sum_{d_i < 0} R_i \right) $$
  !!
  !! The procedure takes into consideration tied ranks.
  module procedure s_tst_signedrank_1s
end interface

! Wilcoxon 2-sample (paired sample) signed rank test
interface fsml_signedrank_paired
  !! The Wilcoxon signed rank test is a non-parametric test that determines
  !! if two related paired samples come from the same distribution.
  !! It can be regarded as a non-parametric version of the paired t-test.
  !!
  !! ### Hypotheses:
  !! The Wilcoxon signed rank test is mathematically equivalent to the
  !! 1-sample Wilcoxon signed rank test conducted on the difference vector
  !! \( d = x_1 - x_2 \) with \( mu_0 \) set to zero. Consequently, the
  !! the null hypothesis \( H_0 \)  can be expressed as:
  !!
  !! Samples \( x_1 - x_2 \) are symmetric around \( \mu = 0 \).
  !!
  !! The default alternative hypothesis \( H_1 \) is two-sided and also be
  !! set explicitly (`h1 = "two"`). It can be expressed as:
  !!
  !! Samples \( x_1 - x_2 \) are symmetric around \( \mu \neq 0 \)
  !!
  !! If the alternative hypothesis is set to "greater than" (`h1 = "gt"`), it is:
  !!
  !! Samples \( x_1 - x_2 \) are symmetric around \( \mu > 0 \)
  !!
  !! If the alternative hypothesis is set to "less than" (`h1 = "lt"`), it is:
  !!
  !! Samples \( x_1 - x_2 \) are symmetric around \( \mu < 0 \)
  !!
  !! The procedure takes into consideration tied ranks.
  module procedure s_tst_signedrank_2s
end interface

! Wilcoxon rank-sum / Mann–Whitney U test
interface fsml_ranksum
  !! The ranks sum test (Wilcoxon rank-sum test or Mann–Whitney U test) is a
  !! non-parametric test to determine if two independent samples \( x_1 \) and
  !! \( x_2 \) are have the same distribution. It can be regarded as the non-parametric
  !! equivalent of the 2-sample t-test.
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
  !! \( H_0 \): the distributions of \( x_1 \) and \( x_2 \) are equal.
  !! \( H_1 \): the distributions of \( x_1 \) and \( x_2 \) are not equal.
  !!
  !! ### Procedure:
  !! The Mann–Whitney U statistic is calculated for each sample as follows:
  !! $$ U_i = R_i - \frac{n_i \cdot (n_i + 1)}{2} $$
  !! where \( R_i \) is the sum of ranks of sample set \( i \)
  !! and \( n_i \) is the sample size of sample set \( i \).
  !! The final U statistic is:
  !! $$ U = \min(U_1, U_2) $$
  !!
  !! The procedure takes into consideration tied ranks.
  module procedure s_tst_ranksum
end interface

! Kruskal Wallis H test
interface fsml_kruskalwallis
  !! The Kruskal-Wallis H-test is used to determine whether samples originate from the same
  !! distribution without assuming normality. It is therefore considered a nonparametric
  !! alternative to the one-way ANOVA (Analysis of Variance).
  !!
  !! ### Hypotheses:
  !! The null hypothesis \( H_0 \) and alternative hypothesis \( H_1 \) are defined as:
  !! \( H_0 \): The populations have the same distribution (medians are equal), and
  !! \( H_1 \): At least one population differs from the others.
  !!
  !! ### Procedure:
  !! The data is passed to the procedure as a rank-2 array `x`, where each column is a group of observations.
  !! All values are ranked across the entire dataset, with tied values assigned the average rank.
  !!
  !! The Kruskal-Wallis H-statistic (`h`) is computed as:
  !! $$
  !! H = \frac{12}{n(n+1)} \sum_{j=1}^{k} \frac{R_j^2}{n_j} - 3(n+1)
  !! $$
  !! where:
  !! - \( n \) is the total number of observations,
  !! - \( k \) is the number of groups,
  !! - \( n_j \) is the number of observations in group \( j \), and
  !! - \( R_j \) is the sum of ranks in group \( j \).
  !!
  !! The degrees of freedom are:
  !! $$ \nu = k - 1 $$
  !! and returned as `df`.
  !!
  !! The p-value (`p`) is computed from the chi-squared distribution:
  !! $$ p = P(\chi^2_{\nu} > H_{observed}) = 1 - \text{CDF}(H_{observed}) $$
  !!
  !! It is computed using the elemental procedure `f_dst_chi2_cdf_core`.
  !!
  !! The Kruskal-Wallis test assumes that:
  !! a) all groups are independent,
  !! b) the response variable is ordinal or continuous,
  !! c) the group distributions have the same shape,
  !! and d) observations are independent both within and between groups.
  module procedure s_tst_kruskalwallis
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Linear Procedures

! PCA
interface fsml_pca
  !! Principal Component Analysis (PCA) is a procedure to reduce the dimensionality
  !! of multivariate data by identifying a set of orthogonal vectors (eigenvectores)
  !! that represent directions of maximum variance in the dataset.
  !!
  !! The procedure `fsml_pca` is a wrapper for `fsml_eof` and offers a simpler,
  !! more familiar interface for non-geoscientists. The EOF interface allows for
  !! more options to be passed that are irrelevant to standard applications of PCA.
  !! The PCA procedure calls the EOF procedures with weights (`wt`) set to *1.0*,
  !! and matrix options set to `opt = 0` to force the use of the covariance matrix
  !! to be comparable to other common implementations of a PCA (e.g., sklearn).
  !!
  !! The covariance matrix \( \mathbf{C} \) is computed as:
  !! $$
  !! \mathbf{C} = \frac{1}{m - 1} \mathbf{X}^\top \mathbf{X}
  !! $$
  !! where \( \mathbf{X} \) is the preprocessed (centred and optionally standardised) data matrix,
  !! and \( m \) is the number of observations (rows in `x`).
  !!
  !! A symmetric eigen-decomposition is then performed:
  !! $$
  !! \mathbf{C} \mathbf{E} = \mathbf{E} \Lambda
  !! $$
  !! where \( \mathbf{E} \) contains the EOFs (`ev`), and \( \Lambda \) is a diagonal matrix
  !! of eigenvalues (`ew`).
  !!
  !! The principal components or scores (PCs, `pc`) are given by:
  !! $$
  !! \mathbf{PC} = \mathbf{X} \mathbf{E}
  !! $$
  !! The number of valid PC modes is determined by the number of non-zero eigenvalues.
  !! Arrays are initialised to zero and populated only where eigenvalues are strictly positive.
  !!
  !! The explained variance (`r2`) for each component is computed as a fraction:
  !! $$
  !! r^2_j = \frac{\lambda_j}{\sum_k \lambda_k}
  !! $$
  !! where \( j \) is the PC index, and \( k \) spans all retained eigenvalues,
  !! representing all principal components that explain variability in the data.
  !!
  !! **Note:** This subroutine uses `eigh` from the `stdlib_linalg` module to compute
  !! eigenvalues and eigenvectors of the symmetric covariance matrix.
  module procedure s_lin_pca
end interface

! EOF Analysis
interface fsml_eof
  !! Empirical Orthogonal Function (EOF) analysis is a procedure to reduce the dimensionality
  !! of multivariate data by identifying a set of orthogonal vectors (EOFs or eigenvectores)
  !! that represent directions of maximum variance in the dataset.
  !! The term *EOF analysis* is often used interchangably with the geographically weighted
  !! principal component analysis (PCA). The procedures are mathematically equivalent, but
  !! procedures for EOF analysis offer some additional options that are mostly relevant for
  !! geoscience. The procedure `fsml_pca` is a wrapper for `fsml_eof` that offers a simpler,
  !! more familiar interface for non-geoscientists.
  !!
  !! For a classic EOF analysis, the input matrix `x` holds data or observations that have been
  !! discretised in time and space. Rows (`m`) and columns (`n`) can therefore be interpreted
  !! as time and space dimensions, respectively. EOF analysis allows for geographical weighting,
  !! which translates to column-wise weighting prior to analysis in the procedure.
  !! Weights can be set by bassing the rank-1 array `wt` of dimension `n`. If this optional
  !! argument is not passed, the procedure will default to equal weights of value \( wt=1/n \).
  !! It is numerically more stable than *1.0*, which is the default for many implementations of a PCA.
  !!
  !! After the weighting is applied, the covariance or correlation matrix \( \mathbf{C} \) is computed:
  !! $$
  !! \mathbf{C} = \frac{1}{m - 1} \mathbf{X}^\top \mathbf{X}
  !! $$
  !! where \( \mathbf{X} \) is the preprocessed (centred and optionally standardised) data matrix,
  !! and \( m \) is the number of observations (rows in `x`).
  !! The value of the optional argument `opt` determines if the covariance matrix (`opt = 0`) or
  !! correlation matrix (`opt = 1`) is constructed. If the argument is not passed, the procedure will
  !! default to the use of the covariance matrix, as is the standard for a regular PCA.
  !!
  !! A symmetric eigen-decomposition is then performed:
  !! $$
  !! \mathbf{C} \mathbf{E} = \mathbf{E} \Lambda
  !! $$
  !! where \( \mathbf{E} \) contains the EOFs (`eof`), and \( \Lambda \) is a diagonal matrix
  !! of eigenvalues (`ew`).
  !!
  !! The principal components or scores (PCs, `pc`) are given by:
  !! $$
  !! \mathbf{PC} = \mathbf{X} \mathbf{E}
  !! $$
  !! The number of valid EOF/PC modes is determined by the number of non-zero eigenvalues.
  !! Arrays are initialised to zero and populated only where eigenvalues are strictly positive.
  !!
  !! The explained variance (`r2`) for each component is computed as a fraction:
  !! $$
  !! r^2_j = \frac{\lambda_j}{\sum_k \lambda_k}
  !! $$
  !! where \( j \) is the PC index, and \( k \) spans all retained eigenvalues,
  !! representing all principal components that explain variability in the data.
  !!
  !! EOFs may optionally be scaled (`eof_scaled`) for more convenient plotting:
  !! $$
  !! \text{EOF}_{\text{scaled}} = \text{EOF} \cdot \sqrt{\lambda_j}
  !! $$
  !!
  !! **Note:** This subroutine uses `eigh` from the `stdlib_linalg` module to compute
  !! eigenvalues and eigenvectors of the symmetric covariance matrix.
  module procedure s_lin_eof
end interface

interface fsml_lda_2class
! 2-Class LDA
  !! interface fsml_lda_2class
  !! The 2-class multivariate Linear Discriminant Analysis (LDA) is a statistical
  !! procedure for classification and the investigation and explanation of differences
  !! between two groups (or classes) with regard to their attribute variables.
  !! It quantifies the discriminability of the groups and the contribution of each of
  !! the attribute variables to this discriminability.
  !!
  !! The procedure finds a discriminant function that best separates the two groups.
  !! The function can be expressed as a linear combination of the attribute variables:
  !!
  !! $$
  !! Y = \nu_0 + \nu_1 X_1 + \nu_2 X_2 + \dots + \nu_m X_m + \dots + \nu_M X_M
  !! $$
  !!
  !! where \( Y \) is the discriminant function, \( X_m (m=1...M) \) are the attribute
  !! variables used in evaluating the differences between the groups, \( nu_m (m=1...M) \)
  !! are the discriminant coefficients associated with each variable, \( M \) (`nv`) is
  !! the number of variables, and \( \nu_0 \) is the y-intercept.
  !! (**Note**: Mathematically, it is analogous to a multivariate linear regression function.)
  !!
  !! Each attribute variable \( X_m \) contains elements \( x_{mn} (n=1…N) \) (`x`), where
  !! \( N \) (`nd`) is the number of elements in each group. Each element is associated with a
  !! discriminant value \( y_n \) described by:
  !!
  !! $$
  !! y_n = \nu_1 x_{1n} + \nu_2 x_{2n} + \dots + \nu_m x_{mn} + \dots + \nu_M x_{Mn}
  !! $$
  !!
  !! Geometrically, this can be visualised as elements \( y_n \) being projected on the
  !! discriminant axis \( Y \). The optimal discriminant function is then determined by
  !! finding an axis, on which the projected elements for the two groups are best separated.
  !! The best separation is given by maximising the discriminant criterion \( \Gamma \) (`g`),
  !! a signal to noise ratio, so that:
  !!
  !! $$
  !! \Gamma = \frac{\text{scatter between groups }}{\text{scatter within groups }}
  !! = \frac{(\bar{y}_{G1} - \bar{y}_{G2})^2}
  !! {\sum_{j=1}^{n_1} (y_{G1j} - \bar{y}_{G1})^2 + \sum_{j=1}^{n_2} (y_{G2j} - \bar{y}_{G2})^2}
  !! \rightarrow \max
  !! $$
  !!
  !! where \( n_1 \) and \( n_2 \) are the number of elements in groups \( G1 \) and \( G2 \),
  !! respectively. The procedure assumes that these are the same (`nd`) and only accepts 2 groups (`nc = 2`).
  !!
  !! The discriminant coefficients are then standardised (`sa`) using the standard deviations
  !! of respective variables. The discriminant function represents a model that best seperates
  !! the groups and can be used as a classification model. The skill of that model is determined
  !! by forgetting the association of each element with the groups and using the model to reclassify
  !! the elements. The score (`score`) is the fraction of correct classifications and can be
  !! interpreted as a measure of how well the function works as a classification model.
  !!
  !! The procedure optionally returns the Mahalanobis distance (`mh`) as a measure of distance
  !! between the groups.
  !!
  !! **Note:** This subroutine uses `eigh` from the `stdlib_linalg` module.
  module procedure s_lin_lda_2c
end interface

! OLS
interface fsml_ols
  !! The multiple linear Ordinary Least Squares (OLS) regression models the relationship
  !! or linear dependence between a dependent (predictand) variable and and one or more
  !! independent (predictor) variables. The procedure estimates the linear regression
  !! coefficients by minimising the sum of squared residuals.
  !!
  !! The estimated regression model is of the form:
  !!
  !! $$
  !! y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_m x_m + \dots + \beta_M x_M
  !! $$
  !!
  !! where \( y \) is the predictand variable, \( x_m \ (m = 1 \dots M) \) are the predictor variables (`x`)
  !! with `nd` observations, \( \beta_0 \) is the y-intercept (`b0`), \( \beta_m \ (m = 1 \dots M) \) (`b`)
  !! are the regression coefficients, and \( M \) (`nv`) is the number of predictors (excluding the intercept).
  !!
  !! The subroutine constructs a full matrix internally by prepending a column of ones to account for
  !! the intercept. The regression coefficients are estimated as:
  !!
  !! $$
  !! \hat{\beta} = (X^\top X)^{-1} X^\top y
  !! $$
  !!
  !! where \( X \) is the extended design matrix including the intercept term.
  !!
  !! The coefficient of determination \( R^2 \) (`r2`) which represents the proportion of the total
  !! variance of \( y \)  (`y`) explained by the predictors. The predicted values (`y_hat`),
  !! standard errors (`se`) of the coefficients, and the covariance matrix of the predictors (`cov_b`)
  !! can optionally be returned by the procedure, too.
  !!
  !! **Note:** This subroutine uses `eigh` from the `stdlib_linalg` module.
  !! **Note:** The intercept and predictor coefficients are computed separately and returned explicitly.
  module procedure s_lin_ols
end interface

! Ridge
interface fsml_ridge
  !! The multiple linear Ridge regression models the relationship
  !! or linear dependence between a dependent (predictand) variable and one or more
  !! independent (predictor) variables, incorporating a penalty term on the size of the
  !! regression coefficients to reduce multicollinearity and overfitting.
  !!
  !! The procedure estimates the linear regression coefficients by minimising the sum of squared
  !! residuals plus a penalty proportional to the square of the magnitude of coefficients:
  !!
  !! $$
  !! \hat{\beta} = (X^\top X + \lambda I)^{-1} X^\top y
  !! $$
  !!
  !! where \( \lambda \) (`lambda`) is the ridge penalty parameter, and \( I \) is the
  !! identity matrix with the first diagonal element corresponding to the intercept set to zero
  !! (no penalty on intercept).
  !!
  !! The estimated regression model is of the form:
  !!
  !! $$
  !! y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \dots + \beta_m x_m + \dots + \beta_M x_M
  !! $$
  !!
  !! where \( y \) is the predictand variable, \( x_m \ (m = 1 \dots M) \) are the predictor variables (`x`)
  !! with `nd` observations, \( \beta_0 \) is the y-intercept (`b0`), \( \beta_m \ (m = 1 \dots M) \) (`b`)
  !! are the ridge regression coefficients, and \( M \) (`nv`) is the number of predictors
  !! (excluding the intercept).
  !!
  !! The subroutine constructs a full matrix internally by prepending a column of ones to account for
  !! the intercept. The coefficient of determination \( R^2 \) (`r2`), predicted values (`y_hat`),
  !! ridge-adjusted standard errors (`se`) of the coefficients, and the ridge-adjusted covariance
  !! matrix of the predictors (`cov_b`) can optionally be returned. The covariance matrix and standard
  !! errors are adjusted for the ridge penalty as:
  !!
  !! $$
  !! \mathrm{cov}(\hat{\beta}) = \sigma^2 (X^\top X + \lambda I)^{-1} X^\top X (X^\top X + \lambda I)^{-1}
  !! $$
  !!
  !! where \( \sigma^2 \) is the residual variance estimate.
  !!
  !! **Note:** This subroutine uses `eigh` from the `stdlib_linalg` module.
  module procedure s_lin_ridge
end interface

! Mahalanobis distance
interface fsml_mahalanobis
  !! Compute Mahalanobis distance between the first two samples (rows) of data.
  !! using covariance estimated from data or the optional covariance matrix cov(m,m).
  !! data dims are (n, m) where each row is a sample (n samples), columns are features (m features).
  module procedure f_lin_mahalanobis
end interface


! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Nonlinear Procedures

! hierarchical clustering (agglomerative)
interface fsml_hcluster
!! Perform agglomerative hierarchical clustering using centroid linkage
!! and the Mahalanobis distance.
!! The resulting cluster centroids can be passed to a separate k-means
!! procedure for refinement.
  module procedure s_nlp_cluster_h
end interface

! kmeans clustering
interface fsml_kmeans
!! Perform k-means clustering using Mahalanobis distance.
!! Starts from initial centroids (e.g. from hierarchical clustering) and iteratively
!! reassigns samples until convergence or maximum iterations reached.
  module procedure s_nlp_cluster_kmeans
end interface


! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Utilities

! ranks elements in an array
interface fsml_rank
  !! Ranks all samples such that the smallest value obtains rank 1
  !! and the largest rank n. Handles tied ranks and assigns average
  !! rank to tied elements within one group of tied elements.
  module procedure s_utl_rank
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- Data Handling

! read csv file into dataframe
interface fsml_read_csv
  !! Read CSV file directly into dataframe.
  module procedure s_dat_read_csv
end interface

end module fsml

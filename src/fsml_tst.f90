module fsml_tst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for common statistical tests.                               |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for common statistical tests.

  ! load modules
  use :: fsml_typ
  use :: fsml_sts
  use :: fsml_dst

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_tst_t1s

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_t1s(x, mu0, t, df, p, h1)

! ==== Description
!! 1 sample t-test, testing if the sample mean has the value specified in the null hypothesis.
!! The test statstic \( t \) is calculated as follows:
!! $$ t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$
!! where \( \bar{x} \) is the sample mean,
!! \( s \) is the sample standard deviation,
!! \( n \) is the sample size, and
!! \( \mu_0 \) is the population mean.
!!
!! The degrees of freedom \( \nu \) are:
!! $$ \nu = n -1 $$
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_{0} \): \( \bar{x}  =  \mu_0 \), and \( H_{1} \): \( \bar{x} \neq \mu_0 \)

! ==== Declarations
  real(wp)         , intent(in)           :: x(:) !! x vector (samples)
  real(wp)         , intent(in)           :: mu0  !! population mean (null hypothesis expected value)
  character(len=*) , intent(in), optional :: h1   !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out)          :: t    !! test statistic
  integer(i4)      , intent(out)          :: df   !! degrees of freedom
  real(wp)         , intent(out)          :: p    !! p-value
  character(len=16)                       :: w_h1 !! final value for h1
  real(wp)                                :: xbar !! sample mean
  real(wp)                                :: s    !! sample standard deviation

! ==== Instructions

! ---- handle input

  ! assume two-sided if option not passed
  if (present(h1)) then
     w_h1 = h1
  else
     w_h1 = "two"
  endif

! ---- conduct test

  ! get sample standard deviation (using n-1)
  xbar = f_sts_mean(x)
  s = sqrt( sum( (x - xbar) * (x - xbar) ) / real( (size(x)-1), kind=wp ) )

  ! get test statistic
  t = f_tst_t1s_t( f_sts_mean(x), s, size(x), mu0 )

  ! get degrees of freedom
  df = size(x) - 1

  ! get p-value
  select case (w_h1)
     case("lt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "left")
     case("gt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "right")
     case("two")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "two")
  end select

  contains

  ! --------------------------------------------------------------- !
  elemental function f_tst_t1s_t(x, s, n, mu0) result(t)

  ! ==== Description
  !! Calculates the test statstic \( t \) for 1 sample t-test.
  ! NOTE: Alternatively offer as a separate public procedure.

  ! ==== Declarations
  real(wp)   , intent(in) :: x   !! sample mean
  real(wp)   , intent(in) :: s   !! sample standard deviation
  integer(i4), intent(in) :: n   !! sample size
  real(wp)   , intent(in) :: mu0 !! population mean
  real(wp)                :: t   !! test statistic

  ! ==== Instructions
  t = (x - mu0) / ( s / sqrt( real(n, kind=wp) ) )

  end function f_tst_t1s_t

end subroutine s_tst_t1s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_t2s(x1, x2, t, df, p, h1)

! ==== Description
!! Welch's t-test (the generalised version for the 2 sample t-test),
!! testing if two population means \( \mu_1 \) and \( \mu_2\) are the same.
!! The test does not assume equal variances.
!! The test statstic \( t \) is calculated as follows:
!! $$t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{\frac{s^2_1}{n_1} + \frac{s^2_2}{n_2} }} $$
!! where \( \bar{x}_1 \) and \( \bar{x}_2 \) are the sample means
!! \( s^2_1 \) and \( s^2_2 \) are the sample standard deviations, and
!! \( n_1 \) and \( n_1 \) are the sample sizes.
!!
!! The degrees of freedom \( \nu \) is approximated as follows:
!! $$ \nu = \frac{\left( \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2} \right)^2} {\frac{\left( \frac{s_1^2}{n_1} \right)^2}{n_1 - 1} + \frac{\left( \frac{s_2^2}{n_2} \right)^2}{n_2 - 1}} $$
!!
!! If the variances are equal, it is the equivalent of the 2 sample t-test for equal variances:
!! $$ t = \frac{\bar{x}_1 - \bar{x}_2}{s_p \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}} $$
!!
!! With \( \nu \) degrees of freedom:
!! $$ \nu = 2 \cdot n - 2 $$
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_{0} \): \( \mu_1 \ = \mu_2 \), and \( H_{1} \): \( \mu_1 \ \neq \mu_2\)

! ==== Declarations
  real(wp)        , intent(in)           :: x1(:) !! x vector 1 (of samples)
  real(wp)        , intent(in)           :: x2(:) !! x vector 2 (of samples)
  character(len=*), intent(in), optional :: h1    !! \( H_{1} \) option: two (default), le, ge
  real(wp)        , intent(out)          :: t     !! test statistic
  integer(i4)     , intent(out)          :: df    !! degrees of freedom
  real(wp)        , intent(out)          :: p     !! p-value
  real(wp)                               :: s1    !! sample 1 standard deviation
  real(wp)                               :: s2    !! sample 2 standard deviation

! ==== Instructions
  t = 0.0_wp
  p = 0.0_wp

end subroutine s_tst_t2s




end module fsml_tst

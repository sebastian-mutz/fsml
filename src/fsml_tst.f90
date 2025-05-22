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
  use :: fsml_ini
  use :: fsml_sts
  use :: fsml_dst

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_tst_ttest_1s, s_tst_ttest_paired, s_tst_ttest_2s

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_1s(x, mu0, t, df, p, h1)

! ==== Description
!! The 1-sample t-test determines if the sample mean has the value specified in the null hypothesis.
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_{0} \): \( \bar{x}  =  \mu_0 \), and \( H_{1} \): \( \bar{x} \neq \mu_0 \)
!!
!! The test statstic \( t \) is calculated as follows:
!! $$ t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$
!! where \( \bar{x} \) is the sample mean,
!! \( s \) is the sample standard deviation,
!! \( n \) is the sample size, and
!! \( \mu_0 \) is the population mean.
!!
!! The degrees of freedom \( \nu \) is calculated as follows:
!! $$ \nu = n -1 $$

! ==== Declarations
  real(wp)         , intent(in)           :: x(:) !! x vector (samples)
  real(wp)         , intent(in)           :: mu0  !! population mean (null hypothesis expected value)
  character(len=*) , intent(in), optional :: h1   !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out)          :: t    !! test statistic
  real(wp)         , intent(out)          :: df   !! degrees of freedom
  real(wp)         , intent(out)          :: p    !! p-value
  character(len=16)                       :: h1_w !! final value for h1
  real(wp)                                :: xbar !! sample mean
  real(wp)                                :: s    !! sample standard deviation
  integer(i4)                             :: n    !! sample size

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

! ---- conduct test

  ! get mean, sample size, and sample standard deviation (using n-1)
  xbar = f_sts_mean(x)
  n = size(x)
  s = sqrt( dot_product( (x - xbar), (x - xbar) ) / real( (n-1), kind=wp ) )

  ! get test statistic
  t = f_tst_ttest_1s_t(xbar, s, n, mu0)

  ! get degrees of freedom
  df = real(n, kind=wp) - 1.0_wp

  ! get p-value
  select case(h1_w)
     ! less than
     case("lt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "two")
     ! invalid option
     case default
        p = -1.0_wp
  end select

  contains

  ! --------------------------------------------------------------- !
  pure function f_tst_ttest_1s_t(xbar, s, n, mu0) result(t)

     ! ==== Description
     !! Calculates the test statstic \( t \) for 1 sample t-test.
     ! TODO: Think about making elemental and public for batch processing

     ! ==== Declarations
     real(wp)   , intent(in) :: xbar !! sample mean
     real(wp)   , intent(in) :: s    !! sample standard deviation
     integer(i4), intent(in) :: n    !! sample size
     real(wp)   , intent(in) :: mu0  !! population mean
     real(wp)                :: t    !! test statistic

     ! ==== Instructions
     t = (xbar - mu0) / ( s / sqrt( real(n, kind=wp) ) )

  end function f_tst_ttest_1s_t

end subroutine s_tst_ttest_1s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_paired(x1, x2, t, df, p, h1)

! ==== Description
!! The paired sample t-test (or  dependent sample t-test) determines if
!! the mean difference between two sample sets are zero.
!! It is mathematically equivalent to the 1-sample t-test conducted
!! on the difference vector \( d \) with \( \mu_0 = 0 \).
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_{0} \): \( \bar{d}  =  0 \), and \( H_{1} \): \( \bar{d} \neq 0 \)
!!
!! The test statstic \( t \) is calculated as follows:
!! $$ t = \frac{\bar{d} - 0}{s_d / \sqrt{n}}$$
!! where \( \bar{d} \) is the mean of the differences between the sample sets,
!! \( s_d \) is the standard deviation of the differences, and
!! \( n \) is the number of paired samples.
!!
!! The degrees of freedom \( \nu \) is calculated as follows:
!! $$ \nu = n -1 $$

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:) !! x1 vector (samples)
  real(wp)         , intent(in)           :: x2(:) !! x2 vector (samples); must be same length as x1
  character(len=*) , intent(in), optional :: h1    !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out)          :: t     !! test statistic
  real(wp)         , intent(out)          :: df    !! degrees of freedom
  real(wp)         , intent(out)          :: p     !! p-value
  character(len=16)                       :: h1_w  !! final value for h1

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

! ---- conduct test on difference vector

  ! use procedure for 1-sample t-test and set mu0 to 0
  call s_tst_ttest_1s( (x1 - x2), 0.0_wp, t, df, p, h1)

end subroutine s_tst_ttest_paired




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_2s(x1, x2, t, df, p, eq_var, h1)

! ==== Description
!! The 2-sample t-test determines if two population means \( \mu_1 \) and \( \mu_2\) are the same.
!! The procedure can handle 2-sample t-tests for equal variances and Welch's t-tests for unequal variances.
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_{0} \): \( \mu_1 \ = \mu_2 \), and \( H_{1} \): \( \mu_1 \ \neq \mu_2\)
!!
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

! ==== Declarations
  real(wp)        , intent(in)           :: x1(:)    !! x1 vector (samples)
  real(wp)        , intent(in)           :: x2(:)    !! x2 vector (samples)
  character(len=*), intent(in), optional :: h1       !! \( H_{1} \) option: two (default), le, ge
  logical         , intent(in), optional :: eq_var   !! true if equal variances assumed
  real(wp)        , intent(out)          :: t        !! test statistic
  real(wp)        , intent(out)          :: df       !! degrees of freedom
  real(wp)        , intent(out)          :: p        !! p-value
  character(len=16)                      :: h1_w     !! final value for h1
  logical                                :: eq_var_w !! final value for eq_var
  real(wp)                               :: x1bar    !! sample mean for x1
  real(wp)                               :: x2bar    !! sample mean for x2
  real(wp)                               :: s1       !! sample 1 standard deviation
  real(wp)                               :: s2       !! sample 2 standard deviation
  integer(i4)                            :: n1       !! sample size for x1
  integer(i4)                            :: n2       !! sample size for x2
  real(wp)                               :: s1n, s2n !! equation terms

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! assume unequal variances, overwrite if option is passed
  eq_var_w = .false.
  if (present(eq_var)) eq_var_w = eq_var

! ---- conduct test

  ! get means, sample sizes, and sample standard deviations (using n-1)
  x1bar = f_sts_mean(x1)
  x2bar = f_sts_mean(x2)
  n1 = size(x1)
  n2 = size(x2)
  s1 = sqrt( dot_product( (x1 - x1bar), (x1 - x1bar) ) / &
     & real( (n1-1), kind=wp ) )
  s2 = sqrt( dot_product( (x2 - x2bar), (x2 - x2bar) ) / &
     & real( (n2-1), kind=wp ) )

  ! get test statistic
  if (eq_var_w) then
     t = f_tst_ttest_2s_pooled_t(x1bar, x2bar, s1, s2, n1, n2)
  else
     t = f_tst_ttest_2s_welch_t(x1bar, x2bar, s1, s2, n1, n2)
  endif

  ! get degrees of freedom
  if (eq_var_w) then
     df = real( (n1 + n2 - 2), kind=wp )
  else
     s1n = (s1 * s1) / real( (n1), kind=wp )
     s2n = (s2 * s2) / real( (n2), kind=wp )
     df = ( (s1n + s2n) * (s1n + s2n) ) / ( &
        & (s1n * s1n) / real( (n1-1), kind=wp ) + &
        & (s2n * s2n) / real( (n2-1), kind=wp ) )
  endif

  ! get p-value
  select case(h1_w)
     ! less than
     case("lt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf (t, df, 0.0_wp, 1.0_wp, "two")
     ! invalid option
     case default
        p = -1.0_wp
  end select

  contains

  ! --------------------------------------------------------------- !
  pure function f_tst_ttest_2s_welch_t(x1bar, x2bar, s1, s2, n1, n2) result(t)

     ! ==== Description
     !! Calculates the test statstic \( t \) for 2 sample t-test for unequal variances.
     ! TODO: Think about making elemental and public for batch processing

     ! ==== Declarations
     real(wp)   , intent(in) :: x1bar !! sample mean for x1
     real(wp)   , intent(in) :: x2bar !! sample mean for x2
     real(wp)   , intent(in) :: s1    !! sample 1 standard deviation
     real(wp)   , intent(in) :: s2    !! sample 2 standard deviation
     integer(i4), intent(in) :: n1    !! sample size for x1
     integer(i4), intent(in) :: n2    !! sample size for x2
     real(wp)                :: t     !! test statistic

     ! ==== Instructions
     t = (x1bar - x2bar) / ( sqrt( &
       & ( (s1 * s1) / real(n1, kind=wp)) + &
       & ( (s2 * s2) / real(n2, kind=wp)) ) )

  end function f_tst_ttest_2s_welch_t

  ! --------------------------------------------------------------- !
  pure function f_tst_ttest_2s_pooled_t(x1bar, x2bar, s1, s2, n1, n2) result(t)

     ! ==== Description
     !! Calculates the test statistic \( t \) for a two-sample t-test for equal variances.
     !! The function uses the pooled standard deviation.
     ! TODO: Think about making elemental and public for batch processing

     ! ==== Declarations
     real(wp)   , intent(in) :: x1bar !! sample mean for x1
     real(wp)   , intent(in) :: x2bar !! sample mean for x2
     real(wp)   , intent(in) :: s1    !! sample 1 standard deviation
     real(wp)   , intent(in) :: s2    !! sample 2 standard deviation
     integer(i4), intent(in) :: n1    !! sample size for x1
     integer(i4), intent(in) :: n2    !! sample size for x2
     real(wp)                :: t     !! test statistic
     real(wp)                :: sp    !! pooled variance

     ! ==== Instructions

     ! pooled standard deviation
     sp = sqrt( ( &
        & real(n1 - 1, kind=wp) * (s1 * s1) + &
        & real(n2 - 1, kind=wp) * (s2 * s2) ) / real(n1 + n2 - 2, kind=wp) )

     ! test statistic
     t = (x1bar - x2bar) / ( sp * sqrt( &
       & 1.0_wp / real(n1, kind=wp) + 1.0_wp / real(n2, kind=wp) ) )

  end function f_tst_ttest_2s_pooled_t

end subroutine s_tst_ttest_2s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ranksum(x1, x2, u, p, h1)

! ==== Description
!! The ranks sum test (Wilcoxon rank-sum test or Mann–Whitney U test) is a
!! nonparametric test to determine if two independent samples \( x_1 \) and
!! \( x_2 \) are have the same distribution.
!!
!! The null hypothesis \( H_{0} \) and alternative hypothesis \( H_{1} \) can be written as:
!! \( H_0 \): the distributions of \( x_1 \) and \( x_2 \) are equal.
!! \( H_1 \): the distributions of \( x_1 \) and \( x_2 \) are not equal.

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:)    !! x1 vector (samples)
  real(wp)         , intent(in)           :: x2(:)    !! x2 vector (samples)
  real(wp)         , intent(out)          :: u        !! U statistic
  real(wp)         , intent(out)          :: p        !! p-value
  character(len=*) , intent(in), optional :: h1       !! \( H_{1} \) option: "two" (default), "lt", or "gt"
  character(len=16)                       :: h1_w     !! final value for h1
  integer(i4)                             :: n1       !! sample size for x1
  integer(i4)                             :: n2       !! sample size for x2
  real(wp)         , allocatable          :: x(:)     !! combined x1 and x2
  integer(i4)      , allocatable          :: ranks(:) !! stores ranks
  real(wp)                                :: rx1      !! rank sum for x1
  real(wp)                                :: rx2      !! rank sum for x2
  real(wp)                                :: u1       !! U statistic 1
  real(wp)                                :: u2       !! U statistic 2
  real(wp)                                :: mu       !! mean of U (under H0)
  real(wp)                                :: s        !! standard deviation of U (under H0)
  real(wp)                                :: z        !! z statistic
  integer(i4)                             :: i

! ==== Instructions

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! combine and rank all values
  n1 = size(x1)
  n2 = size(x2)
  allocate(x(n1+n2))
  allocate(ranks(n1+n2))
  x(1:n1)  = x1
  x(n1+1:) = x2

  ! assigns ranks, use stdlib procedure
  call s_rank(x, ranks)

  ! sum ranks of sample x
  rx1 = sum(real(ranks(1:n1),  kind=wp))
  rx2 = sum(real(ranks(n1+1:), kind=wp))

  ! deallocate
  deallocate(x)
  deallocate(ranks)

  ! compute U statistics
  u1 = real(n1, kind=wp) * real(n2, kind=wp) + &
     & (real(n1, kind=wp) * (real(n1, kind=wp) + 1.0_wp) / 2.0_wp) - rx1
  u2 = real(n1, kind=wp) * real(n2, kind=wp) + &
     & (real(n2, kind=wp) * (real(n2, kind=wp) + 1.0_wp) / 2.0_wp) - rx2
  u  = min(u1, u2)

  ! mean and standard deviation of U under H0
  mu = n1 * n2 / 2.0_wp
  s  = sqrt(n1 * n2 * (n1 + n2 + 1.0_wp) / 12.0_wp)

  ! z statistic (U is approximately normal for large samples)
  z = (u - mu) / s

  ! get p-value
  select case(h1_w)
     ! less than
     case("lt")
        p = f_dst_norm_cdf(z, 0.0_wp, 1.0_wp, "left")
     ! greater than
     case("gt")
        p = f_dst_norm_cdf(z, 0.0_wp, 1.0_wp, "right")
     ! two-sided
     case("two")
        p = f_dst_norm_cdf(z, 0.0_wp, 1.0_wp, "two")
     ! invalid option
     case default
        p = -1.0_wp
  end select

  contains

  ! --------------------------------------------------------------- !
  pure subroutine s_rank(x, ranks)

     ! ==== Description
     !! Ranks all samples such that the smallest value obtains rank 1
     !! and the largest rank n.

     ! ==== Declarations
     real(wp)                , intent(in)  :: x(:)
     integer(i4), allocatable, intent(out) :: ranks(:)
     integer(i4), allocatable              :: idx(:)
     integer(i4), allocatable              :: sorted(:)
     real(wp)                              :: rank_sum
     integer(i4)                           :: count
     integer(i4)                           :: n, i, j, k

     ! ==== Instructions

     ! allocate
     n = size(x)
     allocate(idx(n))
     allocate(ranks(n))

     ! create index vector
     do i = 1, n
        idx(i) = i
     enddo

     ! sort index based on x
     do i = 2, n
        j = i
        do while (j .gt. 1 .and. x( idx(j) ) .lt. x( idx(j - 1) ))
           k = idx(j)
           idx(j) = idx(j - 1)
           idx(j - 1) = k
           j = j - 1
        enddo
     enddo

     ! assign ranks (with tie averaging)
     i = 1
     do while (i .le. n)
        rank_sum = real(i, kind=wp)
        count = 1
        do j = i + 1, n
           if (x( idx(j) ) .eq. x( idx(i) )) then
              rank_sum = rank_sum + real(j, kind=wp)
              count = count + 1
           else
              exit
           endif
        enddo
        rank_sum = rank_sum / real(count, kind=wp)
        do k = i, i + count - 1
           ranks(idx(k)) = rank_sum
        enddo
        i = i + count
     enddo

     ! deallocate
     deallocate(idx)

  end subroutine s_rank

end subroutine s_tst_ranksum

end module fsml_tst

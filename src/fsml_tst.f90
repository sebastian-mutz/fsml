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
  use :: fsml_utl
  use :: fsml_sts
  use :: fsml_dst

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_tst_ttest_1s, s_tst_ttest_paired, s_tst_ttest_2s
  public :: s_tst_ranksum, s_tst_signedrank_1s, s_tst_signedrank_2s
  ! TODO: wrapper procedures to handle errors, invalid args, etc.

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_1s(x, mu0, t, df, p, h1)

! ==== Description
!! The 1-sample t-test.

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
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="two")
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
!! The paired sample t-test (or  dependent sample t-test).

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
!! The 2-sample t-test.

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
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="two")
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
pure subroutine s_tst_anova_1w(x, f, p)

! ==== Description
!! Performs one-way ANOVA for `k` groups stored in a rank-2 array (`x`).
!! Each column in `x` is a group of observations.

! ==== Declarations
  real(wp), intent(in)  :: x(:,:)    !! 2D array, each column is a group
  real(wp), intent(out) :: f         !! F-statistic
  real(wp), intent(out) :: p         !! p-value from F distribution
  real(wp)              :: df_b      !! between-group degrees of freedom
  real(wp)              :: df_w      !! within-group degrees of freedom
  integer(i4)           :: ni        !! number of elements in groups
  integer(i4)           :: nt        !! number of elements in total
  integer(i4)           :: ng        !! number of groups
  real(wp)              :: mu_t      !! grand/total mean
  real(wp)              :: mu_g      !! group mean
  real(wp)              :: ss_b      !! ss between groups
  real(wp)              :: ss_w      !! ss within groups
  real(wp), allocatable :: x_flat(:) !! flattened x
  integer(i4)           :: i, j      !! flexible integers

! ==== Instructions

  ! flatten all elements to compute grand mean and total n
  allocate( x_flat( product( shape(x) ) ) )
  nt = size(x)
  x_flat = reshape(x, [nt])

  ! get grand mean
  mu_t = f_sts_mean(x_flat(:nt))

  ! initialise sums
  ss_b = 0.0_wp
  ss_w = 0.0_wp

  ! get number of groups
  ng = size(x, 2)

  do j = 1, ng
     ni = size(x, 1)
     mu_g = f_sts_mean(x(:,j))
     ss_b = ss_b + real(ni, kind=wp) * (mu_g - mu_t) ** 2
     ss_w = ss_w + sum( (x(:,j) - mu_g) ** 2 )
  enddo

  ! degrees of freedom
  df_b = real(ng - 1, kind=wp)
  df_w = real(nt - ng, kind=wp)

  ! calculate f statistics
  f = (ss_b / df_b) / (ss_w / df_w)

  ! get right tail p value with F distribution CDF procedure; use default loc, scale
  p = f_dst_f_cdf_core(F, df_b, df_w, 0.0_wp, 1.0_wp, "right")

end subroutine s_tst_anova_1w




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_signedrank_1s(x, mu0, w, p, h1)

! ==== Description
!! The 1-sample Wilcoxon signed rank test.

! ==== Declarations
  real(wp)         , intent(in)           :: x(:)      !! x vector (samples)
  real(wp)         , intent(in)           :: mu0       !! population mean (null hypothesis expected value)
  real(wp)         , intent(out)          :: w         !! W statistic (sum of signed ranks)
  real(wp)         , intent(out)          :: p         !! p-value
  character(len=*) , intent(in), optional :: h1        !! \( H_{1} \): "two" (default), "lt", "gt"
  character(len=16)                       :: h1_w      !! final value for h1
  integer(i4)                             :: n         !! sample size
  integer(i4)                             :: m         !! non-zero sample size
  real(wp)                                :: mr        !! float sample size for x
  real(wp)         , allocatable          :: d(:)      !! differences
  real(wp)         , allocatable          :: dm(:)     !! absolute (modulus of) differences
  integer(i4)      , allocatable          :: ranks(:)  !! ranks of absolute differences
  integer(i4)      , allocatable          :: idx(:)    !! indeces for x > 0
  real(wp)                                :: rpos      !! sum of positive ranks
  real(wp)                                :: rneg      !! sum of negative ranks
  real(wp)                                :: mu        !! expected W under H0
  real(wp)                                :: s         !! standard deviation under H0
  real(wp)                                :: z         !! z statistic
  integer(i4)                             :: i, j

! ==== Instructions

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! get dims and allocate
  n = size(x)
  allocate(d(n))

  ! compute differences and absolute values
  d  = x - mu0

  ! filter out zero differences
  m = count(abs(d) .gt. 0.0_wp)
  if (m .eq. 0) then
     w = 0.0_wp
     p = 1.0_wp
     return
  endif

  ! allocation based on new vector
  allocate(dm(m))
  allocate(ranks(m))
  allocate(idx(m))
  mr = real(m, kind=wp)

  ! extract non-zero entries and their signs
  i = 0
  do j = 1, n
     if (abs(d(j)) .gt. 0.0_wp) then
        i = i + 1
        dm(i)  = abs(d(j))
        idx(i) = j
     endif
  enddo

  ! rank non-zero differences
  call s_utl_rank(dm, ranks)

  ! compute rank sums
  rpos = 0.0_wp
  rneg = 0.0_wp
  do i = 1, m
     if (d(idx(i)) .gt. 0.0_wp) then
        rpos = rpos + real(ranks(i), wp)
     else
        rneg = rneg + real(ranks(i), wp)
     endif
  enddo

  ! W is the smaller of the rank sums
  w = min(rpos, rneg)

  ! expected value and standard deviation under H0
  mu = mr * (mr + 1.0_wp) / 4.0_wp
  s  = sqrt(mr * (mr + 1.0_wp) * (2.0_wp * mr + 1.0_wp) / 24.0_wp)

  ! z statistic
  z = (w - mu) / s

  ! get p-value
  select case(h1_w)
     ! less than
     case("lt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="two")
     ! invalid option
     case default
        p = -1.0_wp
  end select

  ! deallocate
  deallocate(d)
  deallocate(dm)
  deallocate(ranks)
  deallocate(idx)

end subroutine s_tst_signedrank_1s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_signedrank_2s(x1, x2, w, p, h1)

! ==== Description
!! The Wilcoxon signed rank test.

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:)     !! sample 1 (paired data)
  real(wp)         , intent(in)           :: x2(:)     !! sample 2 (paired data)
  real(wp)         , intent(out)          :: w         !! W statistic (sum of signed ranks)
  real(wp)         , intent(out)          :: p         !! p-value
  character(len=*) , intent(in), optional :: h1        !! \( H_{1} \): "two" (default), "lt", "gt"
  character(len=16)                       :: h1_w      !! final value for h1
  real(wp)         , allocatable          :: d(:)      !! differences

! ==== Instructions

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! get dims and allocate
  allocate(d( size(x1) ))

  ! compute differences and absolute values
  d  = x1 - x2

  ! use 1-sample signed-rank test with mu0 = 0
  call s_tst_signedrank_1s(d, 0.0_wp, w, p, h1_w)

  ! deallocate
  deallocate(d)

end subroutine s_tst_signedrank_2s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ranksum(x1, x2, u, p, h1)

! ==== Description
!! The ranks sum test (Wilcoxon rank-sum test or Mann–Whitney U test).

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:)    !! x1 vector (samples)
  real(wp)         , intent(in)           :: x2(:)    !! x2 vector (samples)
  real(wp)         , intent(out)          :: u        !! U statistic
  real(wp)         , intent(out)          :: p        !! p-value
  character(len=*) , intent(in), optional :: h1       !! \( H_{1} \) option: "two" (default), "lt", or "gt"
  character(len=16)                       :: h1_w     !! final value for h1
  integer(i4)                             :: n1       !! sample size for x1
  integer(i4)                             :: n2       !! sample size for x2
  real(wp)                                :: n1r      !! float sample size for x1
  real(wp)                                :: n2r      !! float sample size for x2
  real(wp)         , allocatable          :: x(:)     !! combined x1 and x2
  integer(i4)      , allocatable          :: ranks(:) !! stores ranks
  real(wp)                                :: rx1      !! rank sum for x1
  real(wp)                                :: rx2      !! rank sum for x2
  real(wp)                                :: u1       !! U statistic 1
  real(wp)                                :: u2       !! U statistic 2
  real(wp)                                :: mu       !! expected value of U (under H0)
  real(wp)                                :: s        !! standard deviation of U (under H0)
  real(wp)                                :: z        !! z statistic
  integer(i4)                             :: i

! ==== Instructions

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! get dims and combine and rank all values
  n1 = size(x1)
  n2 = size(x2)
  n1r = real(n1, kind=wp)
  n2r = real(n2, kind=wp)
  allocate(x(n1+n2))
  allocate(ranks(n1+n2))
  x(1:n1)  = x1
  x(n1+1:) = x2

  ! assigns ranks, use stdlib procedure
  call s_utl_rank(x, ranks)

  ! sum ranks of sample x
  rx1 = sum(real(ranks(1:n1),  kind=wp))
  rx2 = sum(real(ranks(n1+1:), kind=wp))

  ! compute U statistics
  u1 = n1r * n2r + (n1r * (n1r + 1.0_wp) / 2.0_wp) - rx1
  u2 = n1r * n2r + (n2r * (n2r + 1.0_wp) / 2.0_wp) - rx2
  u  = min(u1, u2)

  ! expected value and standard deviation of U under H0
  mu = n1r * n2r / 2.0_wp
  s  = sqrt(n1r * n2r * (n1r + n2r + 1.0_wp) / 12.0_wp)

  ! z statistic (U is approximately normal for large samples)
  z = (u - mu) / s

  ! get p-value
  select case(h1_w)
     ! less than
     case("lt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="two")
     ! invalid option
     case default
        p = -1.0_wp
  end select

  ! deallocate
  deallocate(x)
  deallocate(ranks)

end subroutine s_tst_ranksum

end module fsml_tst

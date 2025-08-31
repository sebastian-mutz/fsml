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
  use :: fsml_con
  use :: fsml_err
  use :: fsml_sts
  use :: fsml_dst

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_tst_ttest_1s, s_tst_ttest_paired, s_tst_ttest_2s
  public :: s_tst_anova_1w
  public :: s_tst_ranksum, s_tst_signedrank_1s, s_tst_signedrank_2s
  public :: s_tst_kruskalwallis

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_ttest_1s(x, mu0, t, df, p, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_ttest_1s_core`.

! ==== Declarations
  real(wp)         , intent(in)           :: x(:) !! x vector (samples)
  real(wp)         , intent(in)           :: mu0  !! population mean (null hypothesis expected value)
  character(len=*) , intent(in), optional :: h1   !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out)          :: t    !! test statistic
  real(wp)         , intent(out)          :: df   !! degrees of freedom
  real(wp)         , intent(out)          :: p    !! p-value
  character(len=16)                       :: h1_w !! final value for h1

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure function
  call s_tst_ttest_1s_core(x, mu0, t, df, p, h1_w)

end subroutine s_tst_ttest_1s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_1s_core(x, mu0, t, df, p, h1)

! ==== Description
!! The 1-sample t-test.

! ==== Declarations
  real(wp)         , intent(in)  :: x(:) !! x vector (samples)
  real(wp)         , intent(in)  :: mu0  !! population mean (null hypothesis expected value)
  character(len=*) , intent(in)  :: h1   !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out) :: t    !! test statistic
  real(wp)         , intent(out) :: df   !! degrees of freedom
  real(wp)         , intent(out) :: p    !! p-value
  real(wp)                       :: xbar !! sample mean
  real(wp)                       :: s    !! sample standard deviation
  integer(i4)                    :: n    !! sample size

! ==== Instructions

  ! get mean, sample size, and sample standard deviation (using n-1)
  xbar = f_sts_mean_core(x)
  n = size(x)
  s = sqrt( dot_product( (x - xbar), (x - xbar) ) / real( (n-1), kind=wp ) )

  ! get test statistic
  t = f_tst_ttest_1s_t(xbar, s, n, mu0)

  ! get degrees of freedom
  df = real(n, kind=wp) - 1.0_wp

  ! get p-value
  select case(h1)
     ! less than
     case("lt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="two")
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

end subroutine s_tst_ttest_1s_core






! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_ttest_paired(x1, x2, t, df, p, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_ttest_paired_core`.

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

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x1) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if x1 and x2 have same size
  if (size(x1) .ne. size(x2)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_ttest_paired_core(x1, x2, t, df, p, h1_w)

end subroutine s_tst_ttest_paired




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_paired_core(x1, x2, t, df, p, h1)

! ==== Description
!! The paired sample t-test (or dependent sample t-test).
!! It is a special case of `s_tst_ttest_1s` and uses
!! the same pure procedure (`s_tst_ttest_1s_core`).

! ==== Declarations
  real(wp)         , intent(in)  :: x1(:) !! x1 vector (samples)
  real(wp)         , intent(in)  :: x2(:) !! x2 vector (samples); must be same length as x1
  character(len=*) , intent(in)  :: h1    !! \( H_{1} \) option: two (default), le, ge
  real(wp)         , intent(out) :: t     !! test statistic
  real(wp)         , intent(out) :: df    !! degrees of freedom
  real(wp)         , intent(out) :: p     !! p-value

! ==== Instructions

  ! use procedure for 1-sample t-test on difference vector and set mu0 to 0
  call s_tst_ttest_1s_core( (x1 - x2), 0.0_wp, t, df, p, h1)

end subroutine s_tst_ttest_paired_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_ttest_2s(x1, x2, t, df, p, eq_var, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_ttest_2s_core`.

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

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! assume unequal variances, overwrite if option is passed
  eq_var_w = .false.
  if (present(eq_var)) eq_var_w = eq_var

  ! check if size is valid
  if (size(x1) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     t  = c_sentinel_r
     df = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_ttest_2s_core(x1, x2, t, df, p, eq_var_w, h1_w)

end subroutine s_tst_ttest_2s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ttest_2s_core(x1, x2, t, df, p, eq_var, h1)

! ==== Description
!! The 2-sample t-test.

! ==== Declarations
  real(wp)        , intent(in)  :: x1(:)    !! x1 vector (samples)
  real(wp)        , intent(in)  :: x2(:)    !! x2 vector (samples)
  character(len=*), intent(in)  :: h1       !! \( H_{1} \) option: two (default), le, ge
  logical         , intent(in)  :: eq_var   !! true if equal variances assumed
  real(wp)        , intent(out) :: t        !! test statistic
  real(wp)        , intent(out) :: df       !! degrees of freedom
  real(wp)        , intent(out) :: p        !! p-value
  real(wp)                      :: x1bar    !! sample mean for x1
  real(wp)                      :: x2bar    !! sample mean for x2
  real(wp)                      :: s1       !! sample 1 standard deviation
  real(wp)                      :: s2       !! sample 2 standard deviation
  integer(i4)                   :: n1       !! sample size for x1
  integer(i4)                   :: n2       !! sample size for x2
  real(wp)                      :: s1n, s2n !! equation terms

! ==== Instructions

! ---- conduct test

  ! get means, sample sizes, and sample standard deviations (using n-1)
  x1bar = f_sts_mean_core(x1)
  x2bar = f_sts_mean_core(x2)
  n1 = size(x1)
  n2 = size(x2)
  s1 = sqrt( dot_product( (x1 - x1bar), (x1 - x1bar) ) / &
     & real( (n1-1), kind=wp ) )
  s2 = sqrt( dot_product( (x2 - x2bar), (x2 - x2bar) ) / &
     & real( (n2-1), kind=wp ) )

  ! get test statistic
  if (eq_var) then
     t = f_tst_ttest_2s_pooled_t(x1bar, x2bar, s1, s2, n1, n2)
  else
     t = f_tst_ttest_2s_welch_t(x1bar, x2bar, s1, s2, n1, n2)
  endif

  ! get degrees of freedom
  if (eq_var) then
     df = real( (n1 + n2 - 2), kind=wp )
  else
     s1n = (s1 * s1) / real( (n1), kind=wp )
     s2n = (s2 * s2) / real( (n2), kind=wp )
     df = ( (s1n + s2n) * (s1n + s2n) ) / ( &
        & (s1n * s1n) / real( (n1-1), kind=wp ) + &
        & (s2n * s2n) / real( (n2-1), kind=wp ) )
  endif

  ! get p-value
  select case(h1)
     ! less than
     case("lt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_t_cdf_core(t, df, mu=0.0_wp, sigma=1.0_wp, tail="two")
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

end subroutine s_tst_ttest_2s_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_anova_1w(x, f, df_b, df_w, p)

! ==== Description
!! Impure wrapper procedure for `s_tst_anova_1w_core`.

! ==== Declarations
  real(wp), intent(in)  :: x(:,:) !! 2D array, each column is a group
  real(wp), intent(out) :: f      !! F-statistic
  real(wp), intent(out) :: p      !! p-value from F distribution
  real(wp), intent(out) :: df_b   !! degrees of freedom between groups
  real(wp), intent(out) :: df_w   !! degrees of freedom within groups

! ==== Instructions

! ---- handle input

  ! check that no. of elements and groups (must both be 2 or higher)
  if (size(x, 1) .le. 1 .or. size(x, 2) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     f    = c_sentinel_r
     p    = c_sentinel_r
     df_b = c_sentinel_r
     df_w = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_anova_1w_core(x, f, df_b, df_w, p)

end subroutine s_tst_anova_1w




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_anova_1w_core(x, f, df_b, df_w, p)

! ==== Description
!! One-way ANOVA.

! ==== Declarations
  real(wp), intent(in)  :: x(:,:)    !! 2D array, each column is a group
  real(wp), intent(out) :: f         !! F-statistic
  real(wp), intent(out) :: p         !! p-value from F distribution
  real(wp), intent(out) :: df_b      !! degrees of freedom between groups
  real(wp), intent(out) :: df_w      !! degrees of freedom within groups
  integer(i4)           :: ni        !! number of elements in groups
  integer(i4)           :: n         !! number of elements in total
  integer(i4)           :: k         !! number of groups
  real(wp)              :: mu_t      !! grand/total mean
  real(wp)              :: mu_g      !! group mean
  real(wp)              :: ss_b      !! ss between groups
  real(wp)              :: ss_w      !! ss within groups
  real(wp), allocatable :: x_flat(:) !! flattened x
  integer(i4)           :: i, j      !! flexible integers

! ==== Instructions

  ! flatten all elements to compute grand mean and total n
  n = size(x)
  allocate(x_flat(n))
  x_flat = reshape(x, [n])

  ! get grand mean
  mu_t = f_sts_mean_core(x_flat)

  ! initialise sums
  ss_b = 0.0_wp
  ss_w = 0.0_wp

  ! get number of groups
  k = size(x, 2)

  ni = size(x, 1)
  do j = 1, k
     mu_g = f_sts_mean_core(x(:,j))
     ss_b = ss_b + real(ni, kind=wp) * (mu_g - mu_t) ** 2
     ss_w = ss_w + sum( (x(:,j) - mu_g) ** 2 )
  enddo

  ! degrees of freedom
  df_b = real(k - 1, kind=wp)
  df_w = real(n - k, kind=wp)

  ! calculate f statistics
  f = (ss_b / df_b) / (ss_w / df_w)

  ! get right tail p value with F distribution CDF procedure; use default loc, scale
  p = f_dst_f_cdf_core(f, df_b, df_w, 0.0_wp, 1.0_wp, "right")

end subroutine s_tst_anova_1w_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_signedrank_1s(x, mu0, w, p, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_signedrank_1s_core`.

! ==== Declarations
  real(wp)         , intent(in)           :: x(:) !! x vector (samples)
  real(wp)         , intent(in)           :: mu0  !! population mean (null hypothesis expected value)
  character(len=*) , intent(in), optional :: h1   !! \( H_{1} \): "two" (default), "lt", "gt"
  real(wp)         , intent(out)          :: w    !! W statistic (sum of signed ranks)
  real(wp)         , intent(out)          :: p    !! p-value
  character(len=16)                       :: h1_w !! final value for h1

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     w  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     w  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_signedrank_1s_core(x, mu0, w, p, h1_w)

end subroutine s_tst_signedrank_1s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_signedrank_1s_core(x, mu0, w, p, h1)

! ==== Description
!! The 1-sample Wilcoxon signed rank test.

! ==== Declarations
  real(wp)         , intent(in)  :: x(:)     !! x vector (samples)
  real(wp)         , intent(in)  :: mu0      !! population mean (null hypothesis expected value)
  real(wp)         , intent(out) :: w        !! W statistic (sum of signed ranks)
  real(wp)         , intent(out) :: p        !! p-value
  character(len=*) , intent(in)  :: h1       !! \( H_{1} \): "two" (default), "lt", "gt"
  integer(i4)                    :: n        !! sample size
  integer(i4)                    :: m        !! non-zero sample size
  real(wp)                       :: mr       !! float sample size for x
  real(wp)         , allocatable :: d(:)     !! differences
  real(wp)         , allocatable :: dm(:)    !! absolute (modulus of) differences
  real(wp)         , allocatable :: ranks(:) !! ranks of absolute differences
  integer(i4)      , allocatable :: idx(:)   !! indeces for x > 0
  real(wp)                       :: rpos     !! sum of positive ranks
  real(wp)                       :: rneg     !! sum of negative ranks
  real(wp)                       :: mu       !! expected W under H0
  real(wp)                       :: s        !! standard deviation under H0
  real(wp)                       :: z        !! z statistic
  integer(i4)                    :: i, j

! ==== Instructions

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
        rpos = rpos + ranks(i)
     else
        rneg = rneg + ranks(i)
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
  select case(h1)
     ! less than
     case("lt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="two")
  end select

  ! deallocate
  deallocate(d)
  deallocate(dm)
  deallocate(ranks)
  deallocate(idx)

end subroutine s_tst_signedrank_1s_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_signedrank_2s(x1, x2, w, p, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_signedrank_2s_core`.

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:) !! sample 1 (paired data)
  real(wp)         , intent(in)           :: x2(:) !! sample 2 (paired data)
  character(len=*) , intent(in), optional :: h1    !! \( H_{1} \): "two" (default), "lt", "gt"
  real(wp)         , intent(out)          :: w     !! W statistic (sum of signed ranks)
  real(wp)         , intent(out)          :: p     !! p-value
  character(len=16)                       :: h1_w  !! final value for h1

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     w  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x1) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     w  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if x1 and x2 have same size
  if (size(x1) .ne. size(x2)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     w  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_signedrank_2s_core(x1, x2, w, p, h1_w)

end subroutine s_tst_signedrank_2s




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_signedrank_2s_core(x1, x2, w, p, h1)

! ==== Description
!! The Wilcoxon signed rank test.

! ==== Declarations
  real(wp)        , intent(in)  :: x1(:) !! sample 1 (paired data)
  real(wp)        , intent(in)  :: x2(:) !! sample 2 (paired data)
  character(len=*), intent(in)  :: h1    !! \( H_{1} \): "two" (default), "lt", "gt"
  real(wp)        , intent(out) :: w     !! W statistic (sum of signed ranks)
  real(wp)        , intent(out) :: p     !! p-value
  real(wp)        , allocatable :: d(:)  !! differences

! ==== Instructions

  ! get dims and allocate
  allocate(d( size(x1) ))

  ! compute differences and absolute values
  d  = x1 - x2

  ! use 1-sample signed-rank test with mu0 = 0
  call s_tst_signedrank_1s_core(d, 0.0_wp, w, p, h1)

  ! deallocate
  deallocate(d)

end subroutine s_tst_signedrank_2s_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_ranksum(x1, x2, u, p, h1)

! ==== Description
!! Impure wrapper procedure for `s_tst_ranksum_core`.

! ==== Declarations
  real(wp)         , intent(in)           :: x1(:) !! x1 vector (samples)
  real(wp)         , intent(in)           :: x2(:) !! x2 vector (samples)
  real(wp)         , intent(out)          :: u     !! U statistic
  real(wp)         , intent(out)          :: p     !! p-value
  character(len=*) , intent(in), optional :: h1    !! \( H_{1} \) option: "two" (default), "lt", or "gt"
  character(len=16)                       :: h1_w  !! final value for h1

! ==== Instructions

! ---- handle input

  ! assume two-sided, overwrite if option is passed
  h1_w = "two"
  if (present(h1)) h1_w = h1

  ! check if h1 (tail) options are valid
  if (h1_w .ne. "lt" .and. h1_w .ne. "gt" .and. h1_w .ne. "two") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     u  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x1) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     u  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

  ! check if x1 and x2 have same size
  if (size(x1) .ne. size(x2)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     u  = c_sentinel_r
     p  = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_ranksum_core(x1, x2, u, p, h1_w)

end subroutine s_tst_ranksum




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_ranksum_core(x1, x2, u, p, h1)

! ==== Description
!! The ranks sum test (Wilcoxon rank-sum test or Mann–Whitney U test).

! ==== Declarations
  real(wp)         , intent(in)  :: x1(:)    !! x1 vector (samples)
  real(wp)         , intent(in)  :: x2(:)    !! x2 vector (samples)
  real(wp)         , intent(out) :: u        !! U statistic
  real(wp)         , intent(out) :: p        !! p-value
  character(len=*) , intent(in)  :: h1       !! \( H_{1} \) option: "two" (default), "lt", or "gt"
  integer(i4)                    :: n1       !! sample size for x1
  integer(i4)                    :: n2       !! sample size for x2
  real(wp)                       :: n1r      !! float sample size for x1
  real(wp)                       :: n2r      !! float sample size for x2
  real(wp)         , allocatable :: x(:)     !! combined x1 and x2
  real(wp)         , allocatable :: ranks(:) !! stores ranks
  real(wp)                       :: rx1      !! rank sum for x1
  real(wp)                       :: rx2      !! rank sum for x2
  real(wp)                       :: u1       !! U statistic 1
  real(wp)                       :: u2       !! U statistic 2
  real(wp)                       :: mu       !! expected value of U (under H0)
  real(wp)                       :: s        !! standard deviation of U (under H0)
  real(wp)                       :: z        !! z statistic
  integer(i4)                    :: i

! ==== Instructions

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
  rx1 = sum( ranks(1:n1) )
  rx2 = sum( ranks(n1+1:) )

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
  select case(h1)
     ! less than
     case("lt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="left")
     ! greater than
     case("gt")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="right")
     ! two-sided
     case("two")
        p = f_dst_norm_cdf_core(z, mu=0.0_wp, sigma=1.0_wp, tail="two")
  end select

  ! deallocate
  deallocate(x)
  deallocate(ranks)

end subroutine s_tst_ranksum_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_tst_kruskalwallis(x, h, df, p)

! ==== Description
!! Impure wrapper procedure for `s_tst_kruskalwallis_core`.

! ==== Declarations
  real(wp)   , intent(in)  :: x(:,:) !! 2D array, each column is a group
  real(wp)   , intent(out) :: h      !! Kruskal-Wallis H-statistic
  real(wp)   , intent(out) :: df     !! degrees of freedom (k - 1)
  real(wp)   , intent(out) :: p      !! p-value from chi-squared distribution

! ==== Instructions

! ---- handle input

  ! check that no. of elements and groups (must both be 2 or higher)
  if (size(x, 1) .le. 1 .or. size(x, 2) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     h  = c_sentinel_r
     p  = c_sentinel_r
     df = c_sentinel_r
     return
  endif

! ---- conduct test

  ! call pure procedure
  call s_tst_kruskalwallis_core(x, h, df, p)


end subroutine s_tst_kruskalwallis




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_kruskalwallis_core(x, h, df, p)

! ==== Description
!! Kruskal-Wallis H-test for independent samples. No tie correction.

! ==== Declarations
  real(wp)   , intent(in)  :: x(:,:)     !! 2D array, each column is a group
  real(wp)   , intent(out) :: h          !! Kruskal-Wallis H-statistic
  real(wp)   , intent(out) :: df         !! degrees of freedom (k - 1)
  real(wp)   , intent(out) :: p          !! p-value from chi-squared distribution
  integer(i4)              :: k          !! number of groups
  integer(i4)              :: n          !! total number of samples
  integer(i4)              :: ni         !! group size (assumes balanced)
  real(wp)   , allocatable :: x_flat(:)  !! flattened x array
  real(wp)   , allocatable :: ranks(:)   !! real ranks
  real(wp)   , allocatable :: r_sum(:)   !! sum of ranks per group
  integer(i4)              :: idx        !! group indexing
  integer(i4)              :: i, j       !! flexible integers

! ==== Instructions

  ! get sizes
  k = size(x, 2)
  ni = size(x, 1)
  n = k * ni

  ! allocate ranks
  allocate(x_flat(n))
  allocate(ranks(n))
  allocate(r_sum(k))

  ! flatten x
  x_flat = reshape(x, [n])

  ! compute ranks
  call s_utl_rank(x_flat, ranks)

  ! allocate and compute rank sums for each group
  r_sum = 0.0_wp
  idx = 1
  do j = 1, k
     r_sum(j) = sum( ranks(idx:idx + ni - 1) )
     idx = idx + ni
  end do

  ! compute H-statistic
  h = 0.0_wp
  do j = 1, k
     h = h + ( r_sum(j) ** 2 ) / real(ni, kind=wp)
  end do
  h = 12.0_wp * h / real(n * (n + 1), wp) - 3.0_wp * real(n + 1, kind=wp)

  ! degrees of freedom and p-value
  df = real(k - 1, kind=wp)
  p = f_dst_chi2_cdf_core(h, df, 0.0_wp, 1.0_wp, "right")

  ! deallocate
  deallocate(x_flat)
  deallocate(ranks)
  deallocate(r_sum)

end subroutine s_tst_kruskalwallis_core


end module fsml_tst

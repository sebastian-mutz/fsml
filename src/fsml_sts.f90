module fsml_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for basic (sample) statistics.                              |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic sample statistics.
! TODO: change ddof to ddf for consistency with df

  ! load modules
  use :: fsml_ini
  use :: fsml_con
  use :: fsml_err
  use :: fsml_utl

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: f_sts_mean, f_sts_mean_core
  public :: f_sts_median, f_sts_median_core
  public :: f_sts_var, f_sts_var_core
  public :: f_sts_std, f_sts_std_core
  public :: f_sts_cov, f_sts_cov_core
  public :: f_sts_trend, f_sts_trend_core
  public :: f_sts_pcc, f_sts_pcc_core
  public :: f_sts_scc, f_sts_scc_core

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_mean(x) result(mean)

! ==== Description
!! Impure wrapper function for `f_sts_mean_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: mean   !! arithmetic mean

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     mean = c_sentinel_r
     return
  endif

! ---- compute mean

  ! call pure function
  mean = f_sts_mean_core(x)

end function f_sts_mean


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_mean_core(x) result(mean)

! ==== Description
!! Computes arithmetic mean.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: mean   !! arithmetic mean

! ==== Instructions
  mean = sum(x) / real(size(x), kind=wp)

end function f_sts_mean_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_median(x) result(median)

! ==== Description
!! Impure wrapper function for `f_sts_median_core`.

! ==== Declarations
  real(wp), intent(in)  :: x(:)   !! x vector (assumed size array)
  real(wp)              :: median !! median

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     median = c_sentinel_r
     return
  endif

! ---- compute median

  ! call pure function
  median = f_sts_median_core(x)

end function f_sts_median


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_median_core(x) result(median)

! ==== Description
!! Computes median using s_utl_rank for tie-aware ranking

! ==== Declarations
  real(wp), intent(in)  :: x(:)   !! x vector (assumed size array)
  real(wp)              :: median !! median
  real(wp), allocatable :: rx(:)  !! ranks of x
  integer(i4)           :: n      !! dimension of x
  integer(i4)           :: i1, i2

! ==== Instructions

  ! get array dimension
  n = size(x)

  ! get ranks for x; rank arrays allocated in ranking
  call s_utl_rank(x, rx)

  if (mod(n, 2) .eq. 1) then
     ! If n is odd, the middle rank is (n+1)/2
     i1 = maxloc(rx, mask = (rx .eq. real((n+1)/2, wp)), dim=1)
     median = x(i1)
  else
     ! If n is even, average elements with ranks n/2 and n/2+1
     i1 = maxloc(rx, mask = (rx .eq. real(n, wp)/2.0_wp)         , dim=1)
     i2 = maxloc(rx, mask = (rx .eq. 1.0_wp + real(n, wp)/2.0_wp), dim=1)
     median = 0.5_wp * (x(i1) + x(i2))
  endif

  deallocate(rx)

end function f_sts_median_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_var(x, ddf) result(var)

! ==== Description
!! Impure wrapper function for `f_sts_var_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in), optional :: ddf    !! delta degrees of freedom
  real(wp)                       :: ddf_w  !! final value for ddf
  real(wp)                       :: xbar   !! mean of x
  real(wp)                       :: var    !! variance

! ==== Instructions

! ---- handle input

  ! assume ddf = 0 (population, not sample statistics)
  ddf_w = 0.0_wp
  if (present(ddf)) ddf_w = ddf

  ! check if value is valid
  if (ddf_w .ne. 1.0_wp .and. ddf_w .ne. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     var = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     var = c_sentinel_r
     return
  endif

! ---- compute variance

  ! call pure function
  var  = f_sts_var_core(x, ddf_w)

end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var_core(x, ddf) result(var)

! ==== Description
!! Computes (sample) variance.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: ddf  !! delta degrees of freedom
  real(wp)             :: xbar !! mean of x
  real(wp)             :: var  !! variance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  var  = dot_product( (x - xbar), (x - xbar) ) / &
       & (real(size(x), kind=wp) - ddf)

end function f_sts_var_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_std(x, ddf) result(std)

! ==== Description
!! Impure wrapper function for `f_sts_std_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)  !! x vector (assumed size array)
  real(wp), intent(in), optional :: ddf   !! delta degrees of freedom
  real(wp)                       :: ddf_w !! final value for ddf
  real(wp)                       :: std   !! standard deviation

! ==== Instructions

! ---- handle input

  ! assume ddf = 0 (population, not sample statistics)
  ddf_w = 0.0_wp
  if (present(ddf)) ddf_w = ddf

  ! check if value is valid
  if (ddf_w .ne. 1.0_wp .and. ddf_w .ne. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     std = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     std = c_sentinel_r
     return
  endif

! ---- compute standard deviation

  ! call pure function
  std = sqrt( f_sts_var_core(x, ddf_w) )

end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std_core(x, ddf) result(std)

! ==== Description
!! Computes standard deviation.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: ddf  !! delta degrees of freedom
  real(wp)             :: std  !! standard deviation

! ==== Instructions

  ! call pure function
  std = sqrt( f_sts_var_core(x, 0.0_wp) )

end function f_sts_std_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_cov(x, y, ddf) result(cov)

! ==== Description
!! Impure wrapper function for `f_sts_cov_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in)           :: y(:)   !! y vector (assumed size array)
  real(wp), intent(in), optional :: ddf    !! delta degrees of freedom
  real(wp)                       :: ddf_w  !! final value for ddf
  real(wp)                       :: xbar   !! mean of x
  real(wp)                       :: ybar   !! mean of y
  real(wp)                       :: cov    !! covariance

! ==== Instructions

! ---- handle input

  ! assume ddf = 0 (population, not sample statistics)
  ddf_w = 0.0_wp
  if (present(ddf)) ddf_w = ddf

  ! check if value is valid
  if (ddf_w .ne. 1.0_wp .and. ddf_w .ne. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     cov = c_sentinel_r
     return
  endif

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     cov = c_sentinel_r
     return
  endif

  ! check if x and y have same size
  if (size(x) .ne. size(y)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     cov = c_sentinel_r
     return
  endif

! ---- compute covariance

  ! call pure function
  cov = f_sts_cov_core(x, y, ddf_w)

end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov_core(x, y, ddf) result(cov)

! ==== Description
!! Computes covariance.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: y(:) !! y vector (assumed size array)
  real(wp), intent(in) :: ddf  !! delta degrees of freedom
  real(wp)             :: xbar !! mean of x
  real(wp)             :: ybar !! mean of y
  real(wp)             :: cov  !! covariance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  ybar = f_sts_mean_core(y)
  cov = dot_product( (x - xbar), (y - ybar) ) / &
      & (real(size(x), kind=wp) - ddf)

end function f_sts_cov_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_trend(x, y) result(trend)

! ==== Description
!! Impure wrapper function for `f_sts_trend_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)  !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)  !! y vector (assumed size array)
  real(wp)             :: trend !! trend/regression slope

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     trend = c_sentinel_r
     return
  endif

  ! check if x and y have same size
  if (size(x) .ne. size(y)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     trend = c_sentinel_r
     return
  endif

! ---- compute trend

  ! call pure function
  trend = f_sts_trend_core(x, y)

end function f_sts_trend


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_trend_core(x, y) result(trend)

! ==== Description
!! Computes regression coefficient/trend.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: trend  !! trend/regression slope

! ==== Instructions
  trend = f_sts_cov_core(x, y, 0.0_wp) / f_sts_var_core(x, 0.0_wp)

end function f_sts_trend_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_pcc(x, y) result(corr)

! ==== Description
!! Impure wrapper function for `f_sts_trend_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: y(:) !! y vector (assumed size array)
  real(wp)             :: corr !! Pearson correlation coefficient

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     corr = c_sentinel_r
     return
  endif

  ! check if x and y have same size
  if (size(x) .ne. size(y)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     corr = c_sentinel_r
     return
  endif

! ---- compute Pearson correlation coefficient

  ! call pure function
  corr = f_sts_pcc_core(x,y)

end function f_sts_pcc


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_pcc_core(x, y) result(corr)

! ==== Description
!! Computes Pearson correlation coefficient.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: y(:) !! y vector (assumed size array)
  real(wp)             :: corr !! Pearson correlation coefficient

! ==== Instructions
  corr = f_sts_cov_core(x, y, 0.0_wp) / &
       & sqrt( f_sts_var_core(x, 0.0_wp) * f_sts_var_core(y, 0.0_wp) )

end function f_sts_pcc_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_scc(x, y) result(corr)

! ==== Description
!! Impure wrapper for `f_sts_scc_core`.

! ==== Declarations
  real(wp), intent(in)  :: x(:)  !! x vector (assumed size array)
  real(wp), intent(in)  :: y(:)  !! y vector (assumed size array)
  real(wp)              :: corr  !! Spearman correlation coefficient

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     corr = c_sentinel_r
     return
  endif

  ! check if x and y have same size
  if (size(x) .ne. size(y)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     corr = c_sentinel_r
     return
  endif

! ---- compute Spearman rank correlation coefficient

  ! call pure function
  corr = f_sts_scc_core(x,y)

end function f_sts_scc


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_scc_core(x, y) result(corr)

! ==== Description
!! Computes Spearman rank correlation coefficient between x and y.
!! Uses `f_sts_pcc_core` on ranks.

! ==== Declarations
  real(wp), intent(in)  :: x(:)  !! x vector (assumed size array)
  real(wp), intent(in)  :: y(:)  !! y vector (assumed size array)
  real(wp)              :: corr  !! Spearman correlation coefficient
  real(wp), allocatable :: rx(:) !! ranks of x
  real(wp), allocatable :: ry(:) !! ranks of y

! ==== Instructions

  ! rank both arrays (uses your existing s_utl_rank)
  ! rank arrays allocated in ranking
  call s_utl_rank(x, rx)
  call s_utl_rank(y, ry)

  ! Pearson correlation on ranks
  corr = f_sts_pcc_core(rx, ry)

  ! deallocate
  deallocate(rx)
  deallocate(ry)

end function f_sts_scc_core

end module fsml_sts

module fsml_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for basic statistics.                                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic sample statistics.

  ! load modules
  use :: fsml_ini
  use :: fsml_con
  use :: fsml_err

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: f_sts_mean, f_sts_mean_core
  public :: f_sts_var, f_sts_var_core
  public :: f_sts_std, f_sts_std_core
  public :: f_sts_cov, f_sts_cov_core
  public :: f_sts_trend, f_sts_trend_core
  public :: f_sts_pcc, f_sts_pcc_core

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
impure function f_sts_var(x) result(var)

! ==== Description
!! Impure wrapper function for `f_sts_var_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: xbar   !! mean of x
  real(wp)             :: var    !! variance

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     var = c_sentinel_r
     return
  endif

! ---- compute variance
  var  = f_sts_var_core(x)

end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var_core(x) result(var)

! ==== Description
!! Computes variance.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: xbar   !! mean of x
  real(wp)             :: var    !! variance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  var  = dot_product( (x - xbar), (x - xbar) ) / real(size(x), kind=wp)

end function f_sts_var_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_std(x) result(std)

! ==== Description
!! Impure wrapper function for `f_sts_std_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: std    !! standard deviation

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     std = c_sentinel_r
     return
  endif

! ---- compute standard deviation
  std = sqrt( f_sts_var_core(x) )

end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std_core(x) result(std)

! ==== Description
!! Computes standard deviation.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: std    !! standard deviation

! ==== Instructions
  std = sqrt( f_sts_var_core(x) )

end function f_sts_std_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_cov(x,y) result(cov)

! ==== Description
!! Impure wrapper function for `f_sts_cov_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: xbar   !! mean of x
  real(wp)             :: ybar   !! mean of y
  real(wp)             :: cov    !! covariance

! ==== Instructions

! ---- handle input

  ! check if size is valid
  if (size(x) .le. 1) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     cov = c_sentinel_r
     return
  endif

  ! check if x and y have same size
  if (size(x) .ne. size(x)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     cov = c_sentinel_r
     return
  endif

! ---- compute covariance
  cov = f_sts_cov_core(x,y)

end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov_core(x,y) result(cov)

! ==== Description
!! Computes covariance.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: xbar   !! mean of x
  real(wp)             :: ybar   !! mean of y
  real(wp)             :: cov    !! covariance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  ybar = f_sts_mean_core(y)
  cov = dot_product( (x - xbar), (y - ybar) ) / real(size(x), kind=wp)

end function f_sts_cov_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_trend(x,y) result(trend)

! ==== Description
!! Impure wrapper function for `f_sts_trend_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: trend  !! trend/regression slope

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
  if (size(x) .ne. size(x)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     trend = c_sentinel_r
     return
  endif

! ---- compute trend
  trend = f_sts_trend_core(x,y)

end function f_sts_trend


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_trend_core(x,y) result(trend)

! ==== Description
!! Computes regression coefficient/trend.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: trend  !! trend/regression slope

! ==== Instructions
  trend = f_sts_cov_core(x,y) / f_sts_var_core(x)

end function f_sts_trend_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_pcc(x,y) result(corr)

! ==== Description
!! Impure wrapper function for `f_sts_trend_core`.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: corr   !! Pearson correlation coefficient

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
  if (size(x) .ne. size(x)) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(4))
     corr = c_sentinel_r
     return
  endif

! ---- compute Pearson correlation coefficient
  corr = f_sts_pcc_core(x,y)

end function f_sts_pcc


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_pcc_core(x,y) result(corr)

! ==== Description
!! Computes Pearson correlation coefficient.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: corr   !! Pearson correlation coefficient

! ==== Instructions
  corr = f_sts_cov_core(x,y) / sqrt( f_sts_var_core(x) * f_sts_var_core(y) )

end function f_sts_pcc_core

end module fsml_sts

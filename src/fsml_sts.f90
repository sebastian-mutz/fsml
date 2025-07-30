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
impure function f_sts_var(x, ddof) result(var)

! ==== Description
!! Impure wrapper function for `f_sts_var_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in), optional :: ddof   !! delta degrees of freedom
  real(wp)                       :: ddof_w !! final value for ddof
  real(wp)                       :: xbar   !! mean of x
  real(wp)                       :: var    !! variance

! ==== Instructions

! ---- handle input

  ! assume ddof = 0 (population, not sample statistics)
  ddof_w = 0.0_wp
  if (present(ddof)) ddof_w = ddof

  ! check if value is valid
  if (ddof_w .ne. 1.0_wp .and. ddof_w .ne. 0.0_wp) then
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
  var  = f_sts_var_core(x, ddof_w)

end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var_core(x, ddof) result(var)

! ==== Description
!! Computes (sample) variance.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: ddof !! delta degrees of freedom
  real(wp)             :: xbar !! mean of x
  real(wp)             :: var  !! variance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  var  = dot_product( (x - xbar), (x - xbar) ) / &
       & (real(size(x), kind=wp) - ddof)

end function f_sts_var_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_std(x, ddof) result(std)

! ==== Description
!! Impure wrapper function for `f_sts_std_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in), optional :: ddof   !! delta degrees of freedom
  real(wp)                       :: ddof_w !! final value for ddof
  real(wp)                       :: std    !! standard deviation

! ==== Instructions

! ---- handle input

  ! assume ddof = 0 (population, not sample statistics)
  ddof_w = 0.0_wp
  if (present(ddof)) ddof_w = ddof

  ! check if value is valid
  if (ddof_w .ne. 1.0_wp .and. ddof_w .ne. 0.0_wp) then
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
  std = sqrt( f_sts_var_core(x, ddof_w) )

end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std_core(x, ddof) result(std)

! ==== Description
!! Computes standard deviation.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: ddof !! delta degrees of freedom
  real(wp)             :: std  !! standard deviation

! ==== Instructions

  ! call pure function
  std = sqrt( f_sts_var_core(x, 0.0_wp) )

end function f_sts_std_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_sts_cov(x, y, ddof) result(cov)

! ==== Description
!! Impure wrapper function for `f_sts_cov_core`.

! ==== Declarations
  real(wp), intent(in)           :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in)           :: y(:)   !! y vector (assumed size array)
  real(wp), intent(in), optional :: ddof   !! delta degrees of freedom
  real(wp)                       :: ddof_w !! final value for ddof
  real(wp)                       :: xbar   !! mean of x
  real(wp)                       :: ybar   !! mean of y
  real(wp)                       :: cov    !! covariance

! ==== Instructions

! ---- handle input

  ! assume ddof = 0 (population, not sample statistics)
  ddof_w = 0.0_wp
  if (present(ddof)) ddof_w = ddof

  ! check if value is valid
  if (ddof_w .ne. 1.0_wp .and. ddof_w .ne. 0.0_wp) then
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
  cov = f_sts_cov_core(x, y, ddof_w)

end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov_core(x, y, ddof) result(cov)

! ==== Description
!! Computes covariance.

! ==== Declarations
  real(wp), intent(in) :: x(:) !! x vector (assumed size array)
  real(wp), intent(in) :: y(:) !! y vector (assumed size array)
  real(wp), intent(in) :: ddof !! delta degrees of freedom
  real(wp)             :: xbar !! mean of x
  real(wp)             :: ybar !! mean of y
  real(wp)             :: cov  !! covariance

! ==== Instructions
  xbar = f_sts_mean_core(x)
  ybar = f_sts_mean_core(y)
  cov = dot_product( (x - xbar), (y - ybar) ) / &
      & (real(size(x), kind=wp) - ddof)

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
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: corr   !! Pearson correlation coefficient

! ==== Instructions
  corr = f_sts_cov_core(x, y, 0.0_wp) / &
       & sqrt( f_sts_var_core(x, 0.0_wp) * f_sts_var_core(y, 0.0_wp) )

end function f_sts_pcc_core

end module fsml_sts

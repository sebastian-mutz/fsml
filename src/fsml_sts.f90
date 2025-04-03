module fsml_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Test application for fsml lib.                                     |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic statistics.

! load modules
  use :: fsml_typ

! basic options
  implicit none
  private

! declare public procedures
  public :: f_sts_mean, f_sts_var, f_sts_std, f_sts_cov, f_sts_trend, f_sts_corr

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_mean(x) result(mean)
  !! Computes arithmetic mean.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)                :: mean   !! arithmetic mean
  real(wp)                :: sum
  integer(kind=4)         :: i
  sum = 0.0_wp
  do i = 1, size(x)
     sum = sum + x(i)
  enddo
  mean = sum / size(x)
end function f_sts_mean


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var(x) result(var)
  !! Computes variance.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)                :: var    !! variance
  real(wp)                :: sum
  integer(i4)             :: i
  sum = 0.0_wp
  do i = 1, size(x)
     sum = sum + ( x(i) - f_sts_mean(x) ) * (x(i) - f_sts_mean(x) )
  enddo
  var = sum / size(x)
end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std(x) result(std)
  !! Computes standard deviation.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)                :: std    !! standard deviation
  std = sqrt( f_sts_var(x) )
end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov(x,y) result(cov)
  !! Computes covariance. x and y must be the same size.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)   , intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)                :: cov    !! covariance
  real(wp)                :: sum
  integer(i4)             :: i
  sum = 0.0_wp
  do i = 1, size(x)
     sum = sum + ( x(i) - f_sts_mean(x) ) * (y(i) - f_sts_mean(y) )
  enddo
  cov = sum / size(x)
end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_trend(x,y) result(trend)
  !! Computes regression coefficient/trend. x and y must be the same size.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)   , intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)                :: trend  !! trend/regression slope
  trend = f_sts_cov(x,y) / f_sts_var(x)
end function f_sts_trend


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_corr(x,y) result(corr)
  !! Computes Pearson correlation coefficient. x and y must be the same size.
  real(wp)   , intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)   , intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)                :: corr   !! Pearson correlation coefficient
  corr = f_sts_cov(x,y) / sqrt( f_sts_var(x) * f_sts_var(y) )
end function f_sts_corr


end module fsml_sts

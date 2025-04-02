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
pure function f_sts_mean(x,n) result(mean)
  !! Computes arithmetic mean.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)                :: mean
  real(wp)                :: sum
  integer(kind=4)         :: i
  sum = 0.0_wp
  do i = 1, n
     sum = sum + x(i)
  enddo
  mean = sum / n
end function f_sts_mean


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var(x,n) result(var)
  !! Computes variance.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)                :: var
  real(wp)                :: sum
  integer(i4)             :: i
  sum = 0.0_wp
  do i = 1, n
     sum = sum + ( x(i) - f_sts_mean(x,n) ) * (x(i) - f_sts_mean(x,n) )
  enddo
  var = sum / n
end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std(x,n) result(std)
  !! Computes standard deviation.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)                :: std
  std = sqrt( f_sts_var(x,n) )
end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov(x,y,n) result(cov)
  !! Computes covariance.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)   , intent(in) :: y(n)   !! y vector
  real(wp)                :: cov
  real(wp)                :: sum
  integer(i4)             :: i
  sum = 0.0_wp
  do i = 1, n
     sum = sum + ( x(i) - f_sts_mean(x,n) ) * (y(i) - f_sts_mean(y,n) )
  enddo
  cov = sum / n
end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_trend(x,y,n) result(trend)
  !! Computes regression coefficient/trend.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)   , intent(in) :: y(n)   !! y vector
  real(wp)                :: trend
  trend = f_sts_cov(x,y,n) / f_sts_var(x,n)
end function f_sts_trend


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_corr(x,y,n) result(corr)
  !! Computes Pearson correlation coefficient.
  integer(i4), intent(in) :: n      !! vector length
  real(wp)   , intent(in) :: x(n)   !! x vector
  real(wp)   , intent(in) :: y(n)   !! y vector
  real(wp)                :: corr
  corr = f_sts_cov(x,y,n) / sqrt( f_sts_var(x,n) * f_sts_var(y,n) )
end function f_sts_corr


end module fsml_sts

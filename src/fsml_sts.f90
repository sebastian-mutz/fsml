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
  public :: f_sts_mean, f_sts_var, f_sts_cov, f_sts_reg, f_sts_corr

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
function f_sts_mean(x,n) result(mean)
  integer(i4), intent(in) :: n
  real(wp)   , intent(in) :: x(n)
  real(wp)                :: mean
  real(wp)                :: sum
  integer(kind=4)             :: i
  sum = 0.0_wp
  do i = 1, n
     sum = sum + x(i)
  enddo
  mean = sum / n
end function f_sts_mean


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_sts_var(x,n) result(var)
  integer(i4), intent(in) :: n
  real(wp)   , intent(in) :: x(n)
  real(wp)                :: var
  real(wp)                :: sum, mn
  integer(i4)             :: i
  mn = f_sts_mean(x,n)
  sum = 0.0_wp
  do i = 1, n
     sum = sum + (x(i) - mn) * (x(i) - mn)
  enddo
  var = sum / n
end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_sts_cov(x,y,n) result(cov)
  integer(i4), intent(in) :: n
  real(wp)   , intent(in) :: x(n), y(n)
  real(wp)                :: cov
  real(wp)                :: sum, mnx, mny
  integer(i4)             :: i
  mnx = f_sts_mean(x,n)
  mny = f_sts_mean(y,n)
  sum = 0.0_wp
  do i = 1, n
     sum = sum + (x(i) - mnx) * (y(i) - mny)
  enddo
  cov = sum / n
end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_sts_reg(x,y,n) result(reg)
  integer(i4), intent(in) :: n
  real(wp)   , intent(in) :: x(n), y(n)
  real(wp)                :: reg
  reg = f_sts_cov(x,y,n) / f_sts_var(x,n)
end function f_sts_reg


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_sts_corr(x,y,n) result(corr)
  integer(i4), intent(in) :: n
  real(wp)   , intent(in) :: x(n), y(n)
  real(wp)                :: corr
  corr = f_sts_cov(x,y,n) / sqrt( f_sts_var(x,n) * f_sts_var(y,n) )
end function f_sts_corr


end module fsml_sts

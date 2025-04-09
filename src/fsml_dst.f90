module fsml_dst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for statistical distributions.                              |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Statistical distribution module.

! load modules
  use :: fsml_typ
  use :: fsml_con
  use :: stdlib_quadrature, only: gauss_legendre

! basic options
  implicit none
  private

! declare public procedures
  public :: f_dst_pdf_norm, f_dst_cdf_norm, f_dst_pdf_t

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_pdf_norm(x, loc, scale) result(fx)

! ==== Description
!! Probability density function for for normal distribution.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: loc     !! distribution location (mean)
  real(wp), intent(in), optional :: scale   !! distribution dispersion (standard deviation)
  real(wp)                       :: w_loc   !! final value of loc
  real(wp)                       :: w_scale !! final value of scale
  real(wp)                       :: fx

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not passed
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

  ! calculate probability/fx
  fx = (1.0_wp / (w_scale * sqrt(2.0_wp * c_pi))) * &
     & exp(-0.5_wp * ((x - w_loc) / w_scale)**2.0_wp)

end function


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_pdf_t(x, df, loc, scale) result(fx)

! ==== Description
!! Probability density function for for student t distribution.
!! Uses intrinsic gamma function (Fortran 2008 and later)

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp), intent(in), optional :: loc     !! distribution location (~mean)
  real(wp), intent(in), optional :: scale   !! distribution dispersion (~standard deviation)
  real(wp)                       :: w_loc   !! final value of loc
  real(wp)                       :: w_scale !! final value of scale
  real(wp)                       :: fx

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not passed
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

  ! calculate probability/fx
  fx = gamma((df + 1.0_wp) / 2.0_wp) / &
     & (w_scale * sqrt(df * c_pi) * gamma(df / 2.0_wp)) * &
     & (1.0_wp + ( ( (x - w_loc) / w_scale )**2 ) / df)** &
     & (-(df + 1.0_wp) / 2.0_wp)

end function


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_cdf_norm(x, loc, scale) result(p)

! ==== Description
!! Cumulative distribution function for normal distribution.
!! Uses Gauss-Legendre quadrature (stdlib).

! ==== Declarations
  real(wp), intent(in)           :: x           !! sample position
  real(wp), intent(in), optional :: loc         !! distribution location (mean)
  real(wp), intent(in), optional :: scale       !! distribution dispersion (standard deviation)
  real(wp)                       :: w_loc       !! final value of loc
  real(wp)                       :: w_scale     !! final value of scale
  real(wp)                       :: x_arr(c_qn) !! x values array
  real(wp)                       :: w_arr(c_qn) !! weights array
  real(wp)                       :: p           !! returned probability integral
  integer(i4)                    :: i, j

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not passed
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

  ! TODO: determine integral by passed options (x value + 2-tailed or 1-tailed),
  ! then call quadrature subroutine.

  ! gaussian quadrature
  call gauss_legendre(x_arr, w_arr)

  ! find index of closest x value
  j = 1
  do i = 2, c_qn
     if ( abs( x - x_arr(i) ) .le. abs( x - x_arr(j) ) ) j = i
  enddo

  ! compute integral
  p = 0.0_wp
  do i = 1, j
     p = p + f_dst_pdf_norm( x_arr(i), w_loc, w_scale ) * w_arr(i)
  enddo

end function


end module fsml_dst

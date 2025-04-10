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

! basic options
  implicit none
  private

! declare public procedures
  public :: f_dst_pdf_norm, f_dst_cdf_norm, f_dst_pdf_t

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_pdf_norm(x, mu, sigma) result(fx)

! ==== Description
!! Probability density function for for normal distribution.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: mu      !! distribution location (mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)                       :: w_mu    !! final value of mu
  real(wp)                       :: w_sigma !! final value of sigma
  real(wp)                       :: fx

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     w_mu = mu
  else
     w_mu = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     w_sigma = sigma
  else
     w_sigma = 1.0_wp
  endif

  ! calculate probability/fx
  fx = (1.0_wp / (w_sigma * sqrt(2.0_wp * c_pi))) * &
     & exp(-0.5_wp * ((x - w_mu) / w_sigma)**2.0_wp)

end function


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_pdf_t(x, df, mu, sigma) result(fx)

! ==== Description
!! Probability density function for for student t distribution.
!! Uses intrinsic gamma function (Fortran 2008 and later)

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp), intent(in), optional :: mu      !! distribution location (~mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (~standard deviation)
  real(wp)                       :: w_mu    !! final value of mu
  real(wp)                       :: w_sigma !! final value of sigma
  real(wp)                       :: fx

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     w_mu = mu
  else
     w_mu = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     w_sigma = sigma
  else
     w_sigma = 1.0_wp
  endif

  ! calculate probability/fx
  fx = gamma((df + 1.0_wp) / 2.0_wp) / &
     & (w_sigma * sqrt(df * c_pi) * gamma(df / 2.0_wp)) * &
     & (1.0_wp + ( ( (x - w_mu) / w_sigma )**2 ) / df)** &
     & (-(df + 1.0_wp) / 2.0_wp)

end function


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_cdf_norm(x, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for normal distribution.
!! Uses Gauss-Legendre quadrature (stdlib).

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: w_mu    !! final value of mu
  real(wp)                               :: w_sigma !! final value of sigma
  character(len=16)                      :: w_tail  !! final tail option
  real(wp)                               :: z       !! z-score
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     w_mu = mu
  else
     w_mu = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     w_sigma = sigma
  else
     w_sigma = 1.0_wp
  endif

  ! assume two-tailed if not specified
  if (present(tail)) then
     w_tail = tail
  else
     w_tail = "two"
  endif

  ! compute z-score
  z = (x - w_mu) / (w_sigma * sqrt(2.0_wp))

  ! compute integral (left tailed)
  p = 0.5_wp * (1.0_wp + erf(z))

  ! tail options
  select case(w_tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. w_mu) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. w_mu) then
           p = 2.0_wp * p
        endif
     ! confidence level
     case("confidence")
        if (x .gt. w_mu) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. w_mu) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function


end module fsml_dst

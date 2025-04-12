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
  public :: f_dst_pdf_norm, f_dst_cdf_norm, f_dst_pdf_t, f_dst_cdf_t

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

! ---- handle input

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

  ! assume left-tailed if not specified
  if (present(tail)) then
     w_tail = tail
  else
     w_tail = "left"
  endif

! ---- compute CDF

  ! compute z-score
  z = (x - w_mu) / (w_sigma * sqrt(2.0_wp))

  ! compute integral (left tailed)
  p = 0.5_wp * (1.0_wp + erf(z))

  ! tail options
  ! NOTE: alternatively, compare z to 0.0 instead of x to mu
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
     ! confidence interval
     case("confidence")
        if (x .gt. w_mu) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. w_mu) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_cdf_t(t, df, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for student t distribution.
!! Degrees of freedom must be positive.

! ==== Declarations
  real(wp)        , intent(in)           :: t           !! sample position
  integer(i4)     , intent(in)           :: df          !! degrees of freedom
  real(wp)        , intent(in), optional :: mu          !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma       !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail        !! tail options
  real(wp)                               :: w_mu        !! final value of mu
  real(wp)                               :: w_sigma     !! final value of sigma
  character(len=16)                      :: w_tail      !! final tail option
  real(wp)                               :: z           !! z-score
  real(wp)                               :: xbeta, a, b !! parameters for beta function
  real(wp)                               :: p           !! returned probability integral

! ==== Instructions

! ---- handle input

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

  ! assume left-tailed if not specified
  if (present(tail)) then
     w_tail = tail
  else
     w_tail = "left"
  endif

! ---- compute CDF

  ! compute z-score
  z = (t - w_mu) / w_sigma

  ! shape parameters for beta function
  a = 0.5_wp * df
  b = 0.5_wp
  xbeta = df / (df + z**2)

  ! compute integral (left tailed)
  if (z .ge. 0.0_wp) then
    p = 1.0_wp - 0.5_wp * beta_inc(xbeta, a, b)
  else
    p = 0.5_wp * beta_inc(xbeta, a, b)
  endif

  ! tail options
  ! NOTE: alternatively, compare z to 0.0 instead of t to mu
  select case(w_tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (t .gt. w_mu) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (t .le. w_mu) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (t .gt. w_mu) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (t .le. w_mu) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

  contains

  !! NOTE: beta_inc and beta_cf algorithms are thrown together based on public domain Fortran and C code and numerical recipes, and modified to use 2008+ intrinsics. May be improved.

     ! --------------------------------------------------------------- !
     pure function beta_inc(x, a, b) result(betai)

     ! ==== Description
     !! computes the regularised incomplete beta function.

     ! ==== Declarations
     real(wp), intent(in) :: x      !! upper limit of integral
     real(wp), intent(in) :: a, b   !! shape parameters for beta dist.
     real(wp)             :: betai  !! regularised incomplete beta function
     real(wp)             :: lnbeta !! log of beta function
     real(wp)             :: bt     !! pre-multiplier
     real(wp)             :: cf     !! continued fraction

     ! ==== Instructions

     ! handle edge cases
     if (x .le. 0.0_wp) then
        betai = 0.0_wp
        return
     elseif (x .ge. 1.0_wp) then
        betai = 1.0_wp
        return
     elseif (a .le. 0.0_wp .or. b .le. 0.0_wp) then
        betai = 0.0_wp  ! alternative: signal an error
        return
     endif

     ! compute log(beta(a, b)) for numerical stability
     lnbeta = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
     bt = exp(a * log(x) + b * log(1.0_wp - x) - lnbeta)

     ! use symmetry transformation if x > (a+1)/(a+b+2)
     if (x .lt. (a + 1.0_wp) / (a + b + 2.0_wp)) then
       cf = beta_cf(x, a, b)
       betai = bt * cf / a
     else
       cf = beta_cf(1.0_wp - x, b, a)
       betai = 1.0_wp - bt * cf / b
     endif

     end function

     ! --------------------------------------------------------------- !
     pure function beta_cf(x, a, b) result(cf)

     ! ==== Description
     !! computes the continued fraction expansion of incomplete beta function.
     !! Based on Lentz's algorithm (1976)

     ! ==== Declarations
     real(wp), intent(in)   :: x                  !! upper limit of integral
     real(wp), intent(in)   :: a, b               !! shape parameters for beta dist.
     real(wp)               :: cf                 !! continued fraction
     real(wp)               :: c, d
     real(wp)               :: aa, del, qab, qam, qap
     real(wp)   , parameter :: eps   = 1.0e-12_wp !! Convergence threshold ( how close to 1 the fractional delta must be to stop iterating)
     real(wp)   , parameter :: fpmin = 1.0e-30_wp !! small number to prevent division by zero
     integer(i4), parameter :: max_i = 200        !! max. iteration numbers
     integer(i4)            :: m

     ! ==== Instructions

    ! starting conditions
     qab = a + b
     qap = a + 1.0_wp
     qam = a - 1.0_wp

     c = 1.0_wp
     d = 1.0_wp - qab * x / qap
     if (abs(d) .lt. fpmin) d = fpmin
     d = 1.0_wp / d
     cf = d

     ! iterative approx.
     do m = 1, max_i
        ! even term
        aa = m * (b - m) * x / ((qam + 2.0_wp * m) * (a + 2.0_wp * m))
        d = 1.0_wp + aa * d
        if (abs(d) .lt. fpmin) d = fpmin
        c = 1.0_wp + aa / c
        if (abs(c) .lt. fpmin) c = fpmin
        d = 1.0_wp / d
        cf = cf * d * c

        ! odd term
        aa = -(a + m) * (qab + m) * x / ((a + 2.0_wp * m) * (qap + 2.0_wp * m))
        d = 1.0_wp + aa * d
        if (abs(d) .lt. fpmin) d = fpmin
        c = 1.0_wp + aa / c
        if (abs(c) .lt. fpmin) c = fpmin
        d = 1.0_wp / d
        del = d * c
        cf = cf * del

        if (abs(del - 1.0_wp) .lt. eps) exit
     enddo

   end function

end function

end module fsml_dst

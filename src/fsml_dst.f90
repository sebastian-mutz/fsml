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
!! Statistical distribution module. The arguments for functions use common
!! mathematical notation (e.g., mu for the normal distribution mean).

! load modules
  use :: fsml_typ
  use :: fsml_con

! basic options
  implicit none
  private

! declare public procedures
  public :: f_dst_norm_pdf, f_dst_norm_cdf, f_dst_norm_ppf
  public :: f_dst_t_pdf, f_dst_t_cdf, f_dst_t_ppf
  public :: f_dst_gamma_pdf, f_dst_gamma_cdf, f_dst_gamma_ppf
  public :: f_dst_exp_pdf, f_dst_exp_cdf, f_dst_exp_ppf
  public :: f_dst_chi2_pdf, f_dst_chi2_cdf, f_dst_chi2_ppf
  public :: f_dst_gpd_pdf, f_dst_gpd_cdf, f_dst_gpd_ppf

contains

! TODO: handle invalid arguments (e.g., sigma = negative)
! TODO: specify optional arguments passed in ppf procedures to cdf procedures
!       to not rely on guessing by order

! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_norm_pdf(x, mu, sigma) result(fx)

! ==== Description
!! Probability density function for normal distribution.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: mu      !! distribution location (mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)                       :: w_mu    !! final value of mu
  real(wp)                       :: w_sigma !! final value of sigma
  real(wp)                       :: fx

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

! ---- compute PDF

  ! calculate probability/fx
  fx = (1.0_wp / (w_sigma * sqrt(2.0_wp * c_pi))) * &
     & exp(-0.5_wp * ((x - w_mu) / w_sigma)**2.0_wp)

end function f_dst_norm_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_norm_cdf(x, mu, sigma, tail) result(p)

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

end function f_dst_norm_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_norm_ppf(p, mu, sigma) result(x)

! ==== Description
!! Percent point function(PPF) (quantile function or inverse CDF) for normal distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.
!! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
!! will return large negative or positive numbers (highly dependent on the tolerance threshold).

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: mu               !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma            !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: w_mu             !! final value of mu
  real(wp)                          :: w_sigma          !! final value of sigma
  integer(i4), parameter            :: i_max = 200      !! max. iteration numbers
  real(wp)   , parameter            :: tol = 1.0e-12_wp !! p deviation tolerance
  real(wp)                          :: a, b             !! section bounds for bisection algorithm
  real(wp)                          :: x_mid, p_mid     !! x and p mid points in bisection algorithm
  integer(i4)                       :: i                !! for iteration
  real(wp)                          :: x                !! sample position

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

! ---- compute PPF

  ! set initial section
  a = -20.0_wp
  b = 20.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_norm_cdf(x_mid, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value and adjust for mu and sigma
        x = w_mu + w_sigma * x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if p not within valid range or x not found in iterations, set x to 0
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp .or. i .eq. i_max) x = 0.0_wp

end function f_dst_norm_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_t_pdf(x, df, mu, sigma) result(fx)

! ==== Description
!! Probability density function for student t distribution.
!! Uses intrinsic gamma function (Fortran 2008 and later)

! ==== Declarations
  real(wp)   , intent(in)           :: x       !! sample position
  integer(i4), intent(in)           :: df      !! degrees of freedom
  real(wp)   , intent(in), optional :: mu      !! distribution location (~mean)
  real(wp)   , intent(in), optional :: sigma   !! distribution dispersion/scale (~standard deviation)
  real(wp)                          :: w_mu    !! final value of mu
  real(wp)                          :: w_sigma !! final value of sigma
  real(wp)                          :: fx

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

! ---- compute PDF

  ! calculate probability/fx
  fx = gamma((df + 1.0_wp) / 2.0_wp) / &
     & (w_sigma * sqrt(df * c_pi) * gamma(df / 2.0_wp)) * &
     & (1.0_wp + ( ( (x - w_mu) / w_sigma )**2 ) / df)** &
     & (-(df + 1.0_wp) / 2.0_wp)

end function f_dst_t_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_t_cdf(t, df, mu, sigma, tail) result(p)

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
    p = 1.0_wp - 0.5_wp * f_dst_beta_inc(xbeta, a, b)
  else
    p = 0.5_wp * f_dst_beta_inc(xbeta, a, b)
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

end function f_dst_t_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_t_ppf(p, df, mu, sigma) result(x)

! ==== Description
!! Percent point function(PPF) (quantile function or inverse CDF) for t distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.
!! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
!! will return large negative or positive numbers (highly dependent on the tolerance threshold).

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  integer(i4), intent(in)           :: df               !! degrees of freedom
  real(wp)   , intent(in), optional :: mu               !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma            !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: w_mu             !! final value of mu
  real(wp)                          :: w_sigma          !! final value of sigma
  integer(i4), parameter            :: i_max = 200      !! max. iteration numbers
  real(wp)   , parameter            :: tol = 1.0e-12_wp !! p deviation tolerance
  real(wp)                          :: a, b             !! section bounds for bisection algorithm
  real(wp)                          :: x_mid, p_mid     !! x and p mid points in bisection algorithm
  integer(i4)                       :: i                !! for iteration
  real(wp)                          :: x                !! sample position

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

! ---- compute PPF

  ! set initial section bounds
  a = -20.0_wp
  b = 20.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_t_cdf(x_mid, df, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value and adjust for mu and sigma
        x = w_mu + w_sigma * x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if p not within valid range or x not found in iterations, set x to 0
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp .or. i .eq. i_max) x = 0.0_wp

end function f_dst_t_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gamma_pdf(x, alpha, beta, loc) result(fx)

! ==== Description
!! Probability density function for gamma distribution.
!! Uses intrinsic exp function.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: alpha   !! shape  parameter
  real(wp), intent(in), optional :: beta    !! scale parameter
  real(wp), intent(in), optional :: loc     !! location parameter
  real(wp)                       :: w_alpha !! final value for alpha
  real(wp)                       :: w_beta  !! final value for beta
  real(wp)                       :: w_loc   !! final value for loc
  real(wp)                       :: z       !! shifted and scaled variable
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume alpha = 1 if not specified
  if (present(alpha)) then
     w_alpha = alpha
  else
     w_alpha = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     w_beta = beta
  else
     w_beta = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif


! ---- compute PDF

  ! calculate probability/fx
  if (x .le. 0.0_wp .or. w_alpha .le. 0.0_wp .or. w_beta .le. 0.0_wp) then
     fx = 0.0_wp
  else
     z = (x - w_loc) / w_beta
     fx = ( x ** (w_alpha - 1.0_wp) ) * exp( -z) / &
        & ( w_beta ** w_alpha * gamma(w_alpha))
  endif

end function f_dst_gamma_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gamma_cdf(x, alpha, beta, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function for gamma distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in), optional :: alpha   !! shape  parameter
  real(wp)        , intent(in), optional :: beta    !! scale parameter
  real(wp)        , intent(in), optional :: loc     !! location parameter
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: w_alpha !! final value for alpha
  real(wp)                               :: w_beta  !! final value for beta
  real(wp)                               :: w_loc   !! final value for loc
  character(len=16)                      :: w_tail  !! final tail option
  real(wp)                               :: z       !! standardised variable
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume alpha = 1 if not specified
  if (present(alpha)) then
     w_alpha = alpha
  else
     w_alpha = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     w_beta = beta
  else
     w_beta = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     w_tail = tail
  else
     w_tail = "left"
  endif

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .le. w_loc .or. alpha .le. 0.0_wp .or. beta .le. 0.0_wp) then
     p = 0.0_wp
  else
     z = (x - w_loc) / w_beta
     p = f_dst_gamma_inc(alpha, z)
  endif

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
        p = 2.0_wp * (1.0_wp - p)
     ! confidence interval
     case("confidence")
        p = 1.0_wp - 2.0_wp * (1.0_wp - p)
   end select

end function f_dst_gamma_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gamma_ppf(p, alpha, beta, loc) result(x)

! ==== Description
!! Percent point function(PPF) (quantile function or inverse CDF) for gamma distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: alpha            !! shape  parameter
  real(wp)   , intent(in), optional :: beta             !! scale parameter
  real(wp)   , intent(in), optional :: loc              !! location parameter
  real(wp)                          :: w_alpha          !! final value for alpha
  real(wp)                          :: w_beta           !! final value for beta
  real(wp)                          :: w_loc            !! final value for loc
  integer(i4), parameter            :: i_max = 200      !! max. iteration numbers
  real(wp)   , parameter            :: tol = 1.0e-12_wp !! p deviation tolerance
  real(wp)                          :: a, b             !! section bounds for bisection algorithm
  real(wp)                          :: x_mid, p_mid     !! x and p mid points in bisection algorithm
  integer(i4)                       :: i                !! for iteration
  real(wp)                          :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume alpha = 1 if not specified
  if (present(alpha)) then
     w_alpha = alpha
  else
     w_alpha = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     w_beta = beta
  else
     w_beta = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

! ---- compute PPF

  ! set initial section
  a = w_loc
  b = w_loc + w_alpha * w_beta * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_gamma_cdf(x_mid, w_alpha, w_beta, w_loc, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value
        x = x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if p not within valid range or x not found in iterations, set x to 0
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp .or. i .eq. i_max) x = 0.0_wp

end function f_dst_gamma_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_exp_pdf(x, lambda, loc) result(fx)

! ==== Description
!! Probability density function for exponential distribution.
!! Uses intrinsic exp function.

! ==== Declarations
  real(wp), intent(in)           :: x        !! sample position
  real(wp), intent(in), optional :: loc      !! location parameter
  real(wp), intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                       :: w_loc    !! final value for loc
  real(wp)                       :: w_lambda !! final value for lambda
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     w_lambda = lambda
  else
     w_lambda = 1.0_wp
  endif

! ---- compute PDF

  ! calculate probability/fx
  if (x .lt. w_loc) then
     fx = 0.0_wp
  elseif (w_lambda .lt. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = w_lambda * exp( -w_lambda * (x - w_loc) )
  endif

end function f_dst_exp_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_exp_cdf(x, lambda, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function for exponential distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x        !! sample position
  real(wp)        , intent(in), optional :: loc      !! location parameter
  real(wp)        , intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  character(len=*), intent(in), optional :: tail     !! tail options
  real(wp)                               :: w_loc    !! final value for loc
  real(wp)                               :: w_lambda !! final value for lambda
  character(len=16)                      :: w_tail   !! final tail option
  real(wp)                               :: p        !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     w_lambda = lambda
  else
     w_lambda = 1.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     w_tail = tail
  else
     w_tail = "left"
  endif

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .lt. w_loc) then
     p = 0.0_wp
  elseif (w_lambda .lt. 0.0_wp) then
     p = 0.0_wp
  else
     p = 1.0_wp - exp( -lambda * (x - w_loc) )
  endif

  ! tail options
  select case(w_tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
   end select


end function f_dst_exp_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_exp_ppf(p, lambda, loc) result(x)

! ==== Description
!! Percent point function(PPF) (quantile function or inverse CDF) for exponential distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: loc              !! location parameter
  real(wp)   , intent(in), optional :: lambda           !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                          :: w_loc            !! final value for loc
  real(wp)                          :: w_lambda         !! final value for lambda
  integer(i4), parameter            :: i_max = 200      !! max. iteration numbers
  real(wp)   , parameter            :: tol = 1.0e-12_wp !! p deviation tolerance
  real(wp)                          :: a, b             !! section bounds for bisection algorithm
  real(wp)                          :: x_mid, p_mid     !! x and p mid points in bisection algorithm
  integer(i4)                       :: i                !! for iteration
  real(wp)                          :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     w_lambda = lambda
  else
     w_lambda = 1.0_wp
  endif

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = w_loc
  b = w_loc - log(1.0_wp - p) / w_lambda * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_exp_cdf(x_mid, lambda=w_lambda, loc=w_loc, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value
        x = x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if p not within valid range or x not found in iterations, set x to 0
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp .or. i .eq. i_max) x = 0.0_wp

end function f_dst_exp_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_chi2_pdf(x, df, loc, scale) result(fx)

! ==== Description
!! Probability density function for the chi-squared distribution.
!! Uses intrinsic exp and gamma function.

! ==== Declarations
  real(wp)   , intent(in)           :: x       !! sample position
  integer(i4), intent(in)           :: df      !! degrees of freedom
  real(wp)   , intent(in), optional :: loc     !! location parameter
  real(wp)   , intent(in), optional :: scale   !! scale parameter
  real(wp)                          :: w_loc   !! final value for loc
  real(wp)                          :: w_scale !! final value for scale
  real(wp)                          :: fx      !! resulting PDF value
  real(wp)                          :: z       !! standardised variable

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

! ----compute PDF
  z = (x - w_loc) / w_scale
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. w_scale .le. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = (1.0_wp / (2.0_wp ** (df / 2.0_wp) * &
        & gamma(df / 2.0_wp) * w_scale)) * &
        & z ** (df / 2.0_wp - 1.0_wp) * exp(-z / 2.0_wp)
  endif

end function f_dst_chi2_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_chi2_cdf(x, df, loc, scale, tail) result(p)

! ==== Description
!! Cumulative distribution function for the chi-squared distribution.

  ! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  integer(i4)     , intent(in)           :: df      !! degrees of freedom
  real(wp)        , intent(in), optional :: loc     !! location parameter
  real(wp)        , intent(in), optional :: scale   !! scale parameter
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: w_loc   !! final value for loc
  real(wp)                               :: w_scale !! final value for scale
  character(len=16)                      :: w_tail  !! final tail option
  real(wp)                               :: p       !! resulting CDF value
  real(wp)                               :: z       !! standardised variable

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

! ----compute CDF

  ! compute integral (left tailed)
  z = (x - w_loc) / w_scale
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. w_scale .le. 0.0_wp) then
     p = 0.0_wp
  else
     p = f_dst_gamma_inc(df / 2.0_wp, z / 2.0_wp)
  endif

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
        if (x .gt. w_loc) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. w_loc) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. w_loc) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. w_loc) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function f_dst_chi2_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_chi2_ppf(p, df, loc, scale) result(x)

! ==== Description
!! Percent point function (PPF) for the chi-squared distribution.
!! Uses the bisection method for numerical inversion of the CDF.

! ==== Declarations
  real(wp), intent(in)           :: p                !! probability between 0.0 and 1.0
  integer(i4), intent(in)        :: df               !! degrees of freedom
  real(wp), intent(in), optional :: loc              !! location parameter
  real(wp), intent(in), optional :: scale            !! scale parameter
  real(wp)                       :: w_loc            !! final value for loc
  real(wp)                       :: w_scale          !! final value for scale
  integer(i4), parameter         :: i_max = 200      !! maximum iterations
  real(wp), parameter            :: tol = 1.0e-12_wp !! tolerance for convergence
  real(wp)                       :: a, b             !! interval bounds
  real(wp)                       :: x_mid, p_mid     !! midpoint and corresponding CDF value
  integer(i4)                    :: i                !! iteration counter
  real(wp)                       :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     w_loc = loc
  else
     w_loc = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     w_scale = scale
  else
     w_scale = 1.0_wp
  endif

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = w_loc
  b = w_loc + df * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     p_mid = f_dst_chi2_cdf(x_mid, df, loc=w_loc, scale=w_scale) - p
     if (abs(p_mid) .lt. tol) then
        x = x_mid
        return
     else if (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     end if
  end do

  ! if p not within valid range or x not found in iterations, set x to 0
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp .or. i .eq. i_max) x = 0.0_wp

end function f_dst_chi2_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gpd_pdf(x, xi, mu, sigma) result(fx)

! ==== Description
!! Probability density function for generalised pareto distribution.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: xi      !! distribution shape parameter
  real(wp), intent(in), optional :: mu      !! distribution location (mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)                       :: w_mu    !! final value of mu
  real(wp)                       :: w_sigma !! final value of sigma
  real(wp)                       :: z       !! z-score
  real(wp)                       :: fx

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

! ---- compute PDF

  ! compute z-score
  z = (x - w_mu) / w_sigma

  ! calculate probability/fx
  if (xi .eq. 0.0_wp) then
     if (z .lt. 0.0_wp) then
        fx = 0.0_wp
     else
        fx = exp(-z) / sigma
     endif
  else
     if (z .lt. 0.0_wp .or. &
     &  (xi .lt. 0.0_wp .and. z .gt. -1.0_wp / xi)) then
        fx = 0.0_wp
     else
        fx = (1.0_wp + xi * z) ** (-1.0_wp / xi - 1.0_wp) / sigma
     endif
  endif

end function f_dst_gpd_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gpd_cdf(x, xi, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for generalised pareto distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in)           :: xi      !! distribution shape parameter
  real(wp)        , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
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
  z = (x - w_mu) / w_sigma

  ! compute integral (left tailed)
  if (xi .eq. 0.0_wp) then
     if (z .lt. 0.0_wp) then
        p = 0.0_wp
     else
        p = 1.0_wp - exp(-z)
     end if
  else
     if (z .lt. 0.0_wp .or. &
     &  (xi .lt. 0.0_wp .and. z .gt. -1.0_wp / xi)) then
        p = 0.0_wp
     else
        p = 1.0_wp - (1.0_wp + xi * z) ** (-1.0_wp / xi)
     endif
  endif

  ! tail options
  select case(w_tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
   end select

end function f_dst_gpd_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gpd_ppf(p, xi, mu, sigma) result(x)

! ==== Description
!! Percent point function for generalised pareto distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)        , intent(in)           :: p       !! probability between 0.0 - 1.0
  real(wp)        , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)        , intent(in), optional :: xi      !! distribution shape parameter
  real(wp)                               :: w_mu    !! final value of mu
  real(wp)                               :: w_sigma !! final value of sigma
  real(wp)                               :: x       !! sample position

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

! ---- compute PPF

  ! compute inverse cdf based on xi
    if (abs(xi) .lt. 1.0e-12_wp) then
      ! if xi is approximately zero, use exponential distribution
      x = w_mu - w_sigma * log(1.0_wp - p)
    else
      ! if xi is not zero, use general formula
      x = w_mu + (w_sigma / xi) * ( (1.0_wp - p) ** (-xi) - 1.0_wp )
    endif

end function f_dst_gpd_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_gamma_inc(a, x) result(p)

! ==== Description
!! Incomplete gamma function. Needed by gamma and chi-squared cdf.
!! Uses Fortran 2008+ intrinsics.

! ==== Declarations
  real(wp), intent(in)   :: a, x
  real(wp)               :: p
  real(wp)               :: sum, term, lngamma_a
  real(wp)               :: ap, del, b, c, d, h
  real(wp)   , parameter :: eps = 1.0e-12_wp     !! convergence threshold
  real(wp)   , parameter :: fpmin = 1.0e-30_wp   !! small number to prevent division by zero
  integer(i4), parameter :: i_max = 200          !! max. iteration numbers
  integer(i4)            :: i

! ==== Instructions
  if (x .lt. 0.0_wp .or. a .le. 0.0_wp) then
     p = 0.0_wp
     return
  endif

  if (x .eq. 0.0_wp) then
     p = 0.0_wp
     return
  endif

  lngamma_a = log(gamma(a))

  if (x .lt. a + 1.0_wp) then
     ! use series expansion
     sum = 1.0_wp / a
     term = sum
     ap = a
     do
        ap = ap + 1.0_wp
        term = term * x / ap
        sum = sum + term
        if (abs(term) .lt. abs(sum) * eps) exit
     enddo
     p = sum * exp(-x + a * log(x) - lngamma_a)
  else
     ! use continued fraction expansion
     b = x + 1.0_wp - a
     c = 1.0_wp / fpmin
     d = 1.0_wp / b
     h = d
     i = 1
     do
       i = i + 1
       if (mod(i, 2) .eq. 0) then
          del = (i / 2) * (a - (i / 2))
       else
          del = -((i - 1) / 2) * ((i - 1) / 2)
       endif
       b = b + 2.0_wp
       d = del * d + b
       if (abs(d) .lt. fpmin) d = fpmin
       c = b + del / c
       if (abs(c) .lt. fpmin) c = fpmin
       d = 1.0_wp / d
       h = h * d * c
       if (abs(d * c - 1.0_wp) .lt. eps) exit
    enddo
    p = 1.0_wp - exp(-x + a * log(x) - lngamma_a) * h
  endif

end function f_dst_gamma_inc


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_dst_beta_inc(x, a, b) result(betai)

! ==== Description
!! Computes the regularised incomplete beta function. beta_inc and beta_cf
!!  algorithms are based on several public domain Fortran and C code,
!!  Lentz's algorithm (1976), and modified to use 2008+ intrinsics.

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

  ! avoid problems from large a/b
  bt = exp(a * log(x) + b * log(1.0_wp - x) - lnbeta)

  ! switch x with 1-x for num. stability if x > (a+1)/(a+b+2)
  if (x .lt. (a + 1.0_wp) / (a + b + 2.0_wp)) then
    cf = beta_cf(x, a, b)
    betai = bt * cf / a
  else
    cf = beta_cf(1.0_wp - x, b, a)
    betai = 1.0_wp - bt * cf / b
  endif

  contains

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
     real(wp)   , parameter :: eps   = 1.0e-12_wp !! Convergence threshold (how close to 1 the fractional delta must be to stop iterating)
     real(wp)   , parameter :: fpmin = 1.0e-30_wp !! small number to prevent division by zero
     integer(i4), parameter :: i_max = 200        !! max. iteration numbers
     integer(i4)            :: i

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
     do i = 1, i_max
        ! even term
        aa = i * (b - i) * x / ((qam + 2.0_wp * i) * (a + 2.0_wp * i))
        d = 1.0_wp + aa * d
        if (abs(d) .lt. fpmin) d = fpmin
        c = 1.0_wp + aa / c
        if (abs(c) .lt. fpmin) c = fpmin
        d = 1.0_wp / d
        cf = cf * d * c

        ! odd term
        aa = -(a + i) * (qab + i) * x / ((a + 2.0_wp * i) * (qap + 2.0_wp * i))
        d = 1.0_wp + aa * d
        if (abs(d) .lt. fpmin) d = fpmin
        c = 1.0_wp + aa / c
        if (abs(c) .lt. fpmin) c = fpmin
        d = 1.0_wp / d
        del = d * c
        cf = cf * del

        ! check if fractional delta against convergence threshold
        if (abs(del - 1.0_wp) .lt. eps) exit
     enddo

   end function beta_cf

end function f_dst_beta_inc

end module fsml_dst

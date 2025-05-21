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
!! Statistical distribution module. Common parameterisation and mathematical
!! notation is used. Additional parameters are sometimes added to enhance functionality.

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
!       to not rely on guessing by order; also check where loc/scale are passed

! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_pdf(x, mu, sigma) result(fx)

! ==== Description
!! Probability density function for normal distribution.
!! $$ f(x) = \frac{1}{\sigma \cdot \sqrt{2 \cdot \pi}} e^{ -\frac{1}{2} \left( \frac{x - \mu}{\sigma} \right)^2 } $$

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: mu      !! distribution location (mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)                       :: mu_w    !! final value of mu
  real(wp)                       :: sigma_w !! final value of sigma
  real(wp)                       :: z       !! z-score
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

! ---- compute PDF

  ! compute z-score
  z = (x - mu_w) / sigma_w

  ! calculate probability/fx
  fx = (1.0_wp / (sigma_w * sqrt(2.0_wp * c_pi))) * exp( -0.5_wp * (z * z) )

end function f_dst_norm_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_cdf(x, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for normal distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: mu_w    !! final value of mu
  real(wp)                               :: sigma_w !! final value of sigma
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: z       !! z-score
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     tail_w = tail
  else
     tail_w = "left"
  endif

! ---- compute CDF

  ! compute z-score
  z = (x - mu_w) / (sigma_w * sqrt(2.0_wp))

  ! compute integral (left tailed)
  p = 0.5_wp * (1.0_wp + erf(z))

  ! tail options
  ! NOTE: alternatively, compare z to 0.0 instead of x to mu
  select case(tail_w)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. mu_w) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu_w) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. mu_w) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu_w) then
           p = 1.0_wp - 2.0_wp * p
        endif
     ! invalid option
     case default
        p = -1.0_wp
   end select

end function f_dst_norm_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_ppf(p, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for normal distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.
!! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
!! will return large negative or positive numbers (highly dependent on the tolerance threshold).

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: mu               !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma            !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: mu_w             !! final value of mu
  real(wp)                          :: sigma_w          !! final value of sigma
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
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
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
        x = mu_w + sigma_w * x_mid
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
elemental function f_dst_t_pdf(x, df, mu, sigma) result(fx)

! ==== Description
!! Probability density function for student t distribution.
!! Uses intrinsic gamma function (Fortran 2008 and later).
!! $$ f(x) = \frac{\Gamma\left(\frac{\nu + 1}{2}\right)}{\sqrt{\nu \cdot \pi}\, \cdot \Gamma\left(\frac{\nu}{2}\right)} \left(1 + \frac{x^2}{\nu}\right)^{-\frac{\nu + 1}{2}} $$
!! where  \(v\) = degrees of freedom (df) and \(\Gamma\) is the gamma function.

! ==== Declarations
  real(wp)   , intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp)   , intent(in), optional :: mu      !! distribution location (~mean)
  real(wp)   , intent(in), optional :: sigma   !! distribution dispersion/scale (~standard deviation)
  real(wp)                          :: mu_w    !! final value of mu
  real(wp)                          :: sigma_w !! final value of sigma
  real(wp)                          :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

! ---- compute PDF

  ! calculate probability/fx
  fx = gamma((df + 1.0_wp) / 2.0_wp) / &
     & (sigma_w * sqrt(df * c_pi) * gamma(df / 2.0_wp)) * &
     & (1.0_wp + ( ( (x - mu_w) / sigma_w )**2 ) / df)** &
     & (-(df + 1.0_wp) / 2.0_wp)

end function f_dst_t_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_t_cdf(x, df, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for student t distribution.
!! Degrees of freedom must be positive.

! ==== Declarations
  real(wp)        , intent(in)           :: x           !! sample position
  real(wp)     , intent(in)           :: df          !! degrees of freedom
  real(wp)        , intent(in), optional :: mu          !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma       !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail        !! tail options
  real(wp)                               :: mu_w        !! final value of mu
  real(wp)                               :: sigma_w     !! final value of sigma
  character(len=16)                      :: tail_w      !! final tail option
  real(wp)                               :: z           !! z-score
  real(wp)                               :: xbeta, a, b !! parameters for beta function
  real(wp)                               :: p           !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     tail_w = tail
  else
     tail_w = "left"
  endif

! ---- compute CDF

  ! compute z-score
  z = (x - mu_w) / sigma_w

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
  ! NOTE: alternatively, compare z to 0.0 instead of x to mu
  select case(tail_w)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. mu_w) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu_w) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. mu_w) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu_w) then
           p = 1.0_wp - 2.0_wp * p
        endif
     ! invalid option
     case default
        p = -1.0_wp
   end select

end function f_dst_t_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_t_ppf(p, df, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for t distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.
!! Conditions p=0.0 and p=1.0 cannot return negative and positive infinity;
!! will return large negative or positive numbers (highly dependent on the tolerance threshold).

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp), intent(in)           :: df               !! degrees of freedom
  real(wp)   , intent(in), optional :: mu               !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma            !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: mu_w             !! final value of mu
  real(wp)                          :: sigma_w          !! final value of sigma
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
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
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
        x = mu_w + sigma_w * x_mid
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
elemental function f_dst_gamma_pdf(x, alpha, beta, loc) result(fx)

! ==== Description
!! Probability density function for gamma distribution.
!! Uses intrinsic exp function.
!! $$ f(x) = \frac{\lambda^\alpha}{\Gamma(\alpha)} x^{\alpha - 1} e^{-\lambda x}, \quad x > 0, \ \alpha > 0, \ \lambda > 0 $$

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: alpha   !! shape  parameter
  real(wp), intent(in), optional :: beta    !! scale parameter
  real(wp), intent(in), optional :: loc     !! location parameter
  real(wp)                       :: alpha_w !! final value for alpha
  real(wp)                       :: beta_w  !! final value for beta
  real(wp)                       :: loc_w   !! final value for loc
  real(wp)                       :: z       !! shifted and scaled variable
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume alpha = 1 if not specified
  if (present(alpha)) then
     alpha_w = alpha
  else
     alpha_w = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     beta_w = beta
  else
     beta_w = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif


! ---- compute PDF

  ! calculate probability/fx
  if (x .le. 0.0_wp .or. alpha_w .le. 0.0_wp .or. beta_w .le. 0.0_wp) then
     fx = 0.0_wp
  else
     z = (x - loc_w) / beta_w
     fx = ( x ** (alpha_w - 1.0_wp) ) * exp( -z) / &
        & ( beta_w ** alpha_w * gamma(alpha_w))
  endif

end function f_dst_gamma_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_cdf(x, alpha, beta, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for gamma distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in), optional :: alpha   !! shape  parameter
  real(wp)        , intent(in), optional :: beta    !! scale parameter
  real(wp)        , intent(in), optional :: loc     !! location parameter
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: alpha_w !! final value for alpha
  real(wp)                               :: beta_w  !! final value for beta
  real(wp)                               :: loc_w   !! final value for loc
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: z       !! standardised variable
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume alpha = 1 if not specified
  if (present(alpha)) then
     alpha_w = alpha
  else
     alpha_w = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     beta_w = beta
  else
     beta_w = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     tail_w = tail
  else
     tail_w = "left"
  endif

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .le. loc_w .or. alpha .le. 0.0_wp .or. beta .le. 0.0_wp) then
     p = 0.0_wp
  else
     z = (x - loc_w) / beta_w
     p = f_dst_gamma_inc(alpha, z)
  endif

  ! tail options
  select case(tail_w)
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
     ! invalid option
     case default
        p = -1.0_wp
   end select

end function f_dst_gamma_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_ppf(p, alpha, beta, loc) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for gamma distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: alpha            !! shape  parameter
  real(wp)   , intent(in), optional :: beta             !! scale parameter
  real(wp)   , intent(in), optional :: loc              !! location parameter
  real(wp)                          :: alpha_w          !! final value for alpha
  real(wp)                          :: beta_w           !! final value for beta
  real(wp)                          :: loc_w            !! final value for loc
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
     alpha_w = alpha
  else
     alpha_w = 1.0_wp
  endif

  ! assume beta = 1 if not specified
  if (present(beta)) then
     beta_w = beta
  else
     beta_w = 1.0_wp
  endif

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

! ---- compute PPF

  ! set initial section
  a = loc_w
  b = loc_w + alpha_w * beta_w * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_gamma_cdf(x_mid, alpha=alpha_w, beta=beta_w&
          &, loc=loc_w, tail="left") - p
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
elemental function f_dst_exp_pdf(x, lambda, loc) result(fx)

! ==== Description
!! Probability density function for exponential distribution.
!! Uses intrinsic exp function.
!! $$ f(x) = \lambda e^{-\lambda x}, \quad x \geq 0, \ \lambda > 0 $$

! ==== Declarations
  real(wp), intent(in)           :: x        !! sample position
  real(wp), intent(in), optional :: loc      !! location parameter
  real(wp), intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                       :: loc_w    !! final value for loc
  real(wp)                       :: lambda_w !! final value for lambda
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     lambda_w = lambda
  else
     lambda_w = 1.0_wp
  endif

! ---- compute PDF

  ! calculate probability/fx
  if (x .lt. loc_w) then
     fx = 0.0_wp
  elseif (lambda_w .lt. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = lambda_w * exp( -lambda_w * (x - loc_w) )
  endif

end function f_dst_exp_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_exp_cdf(x, lambda, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for exponential distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x        !! sample position
  real(wp)        , intent(in), optional :: loc      !! location parameter
  real(wp)        , intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  character(len=*), intent(in), optional :: tail     !! tail options
  real(wp)                               :: loc_w    !! final value for loc
  real(wp)                               :: lambda_w !! final value for lambda
  character(len=16)                      :: tail_w   !! final tail option
  real(wp)                               :: p        !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     lambda_w = lambda
  else
     lambda_w = 1.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     tail_w = tail
  else
     tail_w = "left"
  endif

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .lt. loc_w) then
     p = 0.0_wp
  elseif (lambda_w .lt. 0.0_wp) then
     p = 0.0_wp
  else
     p = 1.0_wp - exp( -lambda * (x - loc_w) )
  endif

  ! tail options
  select case(tail_w)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! invalid option
     case default
        p = -1.0_wp
   end select


end function f_dst_exp_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_exp_ppf(p, lambda, loc) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for exponential distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: loc              !! location parameter
  real(wp)   , intent(in), optional :: lambda           !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                          :: loc_w            !! final value for loc
  real(wp)                          :: lambda_w         !! final value for lambda
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
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume lambda = 1 if not specified
  if (present(lambda)) then
     lambda_w = lambda
  else
     lambda_w = 1.0_wp
  endif

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = loc_w
  b = loc_w - log(1.0_wp - p) / lambda_w * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_exp_cdf(x_mid, lambda=lambda_w, loc=loc_w, tail="left") - p
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
elemental function f_dst_chi2_pdf(x, df, loc, scale) result(fx)

! ==== Description
!! Probability density function for the chi-squared distribution.
!! Uses intrinsic exp and gamma function.
!! $$ f(x) = \frac{x^{\frac{k}{2} - 1} \cdot e^{-x/2}}{2^{\frac{k}{2}} \cdot \Gamma\left(\frac{k}{2}\right)}, \quad x \geq 0, \ k > 0 $$
!! where \(k\) = degrees of freedom (df) and \(\Gamma\) is the gamma function.

! ==== Declarations
  real(wp)   , intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp)   , intent(in), optional :: loc     !! location parameter
  real(wp)   , intent(in), optional :: scale   !! scale parameter
  real(wp)                          :: loc_w   !! final value for loc
  real(wp)                          :: scale_w !! final value for scale
  real(wp)                          :: fx      !! resulting PDF value
  real(wp)                          :: z       !! standardised variable

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     scale_w = scale
  else
     scale_w = 1.0_wp
  endif

! ----compute PDF
  z = (x - loc_w) / scale_w
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. scale_w .le. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = (1.0_wp / (2.0_wp ** (df / 2.0_wp) * &
        & gamma(df / 2.0_wp) * scale_w)) * &
        & z ** (df / 2.0_wp - 1.0_wp) * exp(-z / 2.0_wp)
  endif

end function f_dst_chi2_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_chi2_cdf(x, df, loc, scale, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for the chi-squared distribution.

  ! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)     , intent(in)           :: df      !! degrees of freedom
  real(wp)        , intent(in), optional :: loc     !! location parameter
  real(wp)        , intent(in), optional :: scale   !! scale parameter
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: loc_w   !! final value for loc
  real(wp)                               :: scale_w !! final value for scale
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: p       !! resulting CDF value
  real(wp)                               :: z       !! standardised variable

! ==== Instructions

! ---- handle input

  ! assume loc = 0 if not specified
  if (present(loc)) then
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     scale_w = scale
  else
     scale_w = 1.0_wp
  endif

! ----compute CDF

  ! compute integral (left tailed)
  z = (x - loc_w) / scale_w
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. scale_w .le. 0.0_wp) then
     p = 0.0_wp
  else
     p = f_dst_gamma_inc(df / 2.0_wp, z / 2.0_wp)
  endif

  ! tail options
  select case(tail_w)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. loc_w) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. loc_w) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. loc_w) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. loc_w) then
           p = 1.0_wp - 2.0_wp * p
        endif
     ! invalid option
     case default
        p = -1.0_wp
   end select

end function f_dst_chi2_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_chi2_ppf(p, df, loc, scale) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for the chi-squared distribution.
!! Uses the bisection method for numerical inversion of the CDF.

! ==== Declarations
  real(wp), intent(in)           :: p                !! probability between 0.0 and 1.0
  real(wp), intent(in)        :: df               !! degrees of freedom
  real(wp), intent(in), optional :: loc              !! location parameter
  real(wp), intent(in), optional :: scale            !! scale parameter
  real(wp)                       :: loc_w            !! final value for loc
  real(wp)                       :: scale_w          !! final value for scale
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
     loc_w = loc
  else
     loc_w = 0.0_wp
  endif

  ! assume scale = 1 if not specified
  if (present(scale)) then
     scale_w = scale
  else
     scale_w = 1.0_wp
  endif

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = loc_w
  b = loc_w + df * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     p_mid = f_dst_chi2_cdf(x_mid, df, loc=loc_w, scale=scale_w, tail="left") - p
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
elemental function f_dst_gpd_pdf(x, xi, mu, sigma) result(fx)

! ==== Description
!! Probability density function for generalised pareto distribution.
!! $$ f(x) = \frac{1}{\sigma} \left( 1 + \frac{\xi (x - \mu)}{\sigma} \right)^{-\frac{1}{\xi} - 1}, \quad x \geq \mu, \ \sigma > 0, \ \xi \in \mathbb{R} $$
!! where \(\xi\) is a shape parameter (xi), \(\sigma\) is the scale parameter (sigma), \(\mu\) (mu) is the location (not mean).

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: xi      !! distribution shape parameter
  real(wp), intent(in), optional :: mu      !! distribution location
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)                       :: mu_w    !! final value of mu
  real(wp)                       :: sigma_w !! final value of sigma
  real(wp)                       :: z       !! z-score
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

! ---- compute PDF

  ! compute z-score
  z = (x - mu_w) / sigma_w

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
elemental function f_dst_gpd_cdf(x, xi, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function \(F(x) = \mathbb{P}(X \leq x)\) for generalised pareto distribution.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in)           :: xi      !! distribution shape parameter
  real(wp)        , intent(in), optional :: mu      !! distribution location
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: mu_w    !! final value of mu
  real(wp)                               :: sigma_w !! final value of sigma
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: z       !! z-score
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

  ! assume left-tailed if not specified
  if (present(tail)) then
     tail_w = tail
  else
     tail_w = "left"
  endif

! ---- compute CDF

  ! compute z-score
  z = (x - mu_w) / sigma_w

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
  select case(tail_w)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! invalid option
     case default
        p = -1.0_wp
   end select

end function f_dst_gpd_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gpd_ppf(p, xi, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function \(Q(p) = {F}_{x}^{-1}(p)\) for generalised pareto distribution.
!! Procedure uses bisection method. p should be between 0.0 and 1.0.

! ==== Declarations
  real(wp)        , intent(in)           :: p       !! probability between 0.0 - 1.0
  real(wp)        , intent(in), optional :: mu      !! distribution location
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)        , intent(in), optional :: xi      !! distribution shape parameter
  real(wp)                               :: mu_w    !! final value of mu
  real(wp)                               :: sigma_w !! final value of sigma
  real(wp)                               :: x       !! sample position

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0 if not passed
  if (present(mu)) then
     mu_w = mu
  else
     mu_w = 0.0_wp
  endif

  ! assume sigma = 1 if not passed
  if (present(sigma)) then
     sigma_w = sigma
  else
     sigma_w = 1.0_wp
  endif

! ---- compute PPF

  ! compute inverse cdf based on xi
  if (abs(xi) .lt. 1.0e-12_wp) then
     ! if xi is approximately zero, use exponential distribution
     x = mu_w - sigma_w * log(1.0_wp - p)
  else
     ! if xi is not zero, use general formula
     x = mu_w + (sigma_w / xi) * ( (1.0_wp - p) ** (-xi) - 1.0_wp )
  endif

end function f_dst_gpd_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_inc(a, x) result(p)

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
elemental function f_dst_beta_inc(x, a, b) result(betai)

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
  ! NOTE: no intrinsic log_beta function; could improve
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
     elemental function beta_cf(x, a, b) result(cf)

     ! ==== Description
     !! computes the continued fraction expansion of incomplete beta function.
     !! Based on Lentz's algorithm (1976)

     ! ==== Declarations
     real(wp), intent(in)   :: x                  !! upper limit of integral
     real(wp), intent(in)   :: a, b               !! shape parameters for beta dist.
     real(wp)               :: cf                 !! continued fraction
     real(wp)               :: c, d
     real(wp)               :: aa, del, qab, qam, qap
     real(wp)   , parameter :: eps   = 1.0e-14_wp !! Convergence threshold (how close to 1 the fractional delta must be to stop iterating)
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

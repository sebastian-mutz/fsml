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
  use :: fsml_ini
  use :: fsml_con
  use :: fsml_err

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: f_dst_norm_pdf, f_dst_norm_pdf_core
  public :: f_dst_norm_cdf, f_dst_norm_cdf_core
  public :: f_dst_norm_ppf, f_dst_norm_ppf_core
  public :: f_dst_t_pdf, f_dst_t_pdf_core
  public :: f_dst_t_cdf, f_dst_t_cdf_core
  public :: f_dst_t_ppf, f_dst_t_ppf_core
  public :: f_dst_gamma_pdf, f_dst_gamma_pdf_core
  public :: f_dst_gamma_cdf, f_dst_gamma_cdf_core
  public :: f_dst_gamma_ppf, f_dst_gamma_ppf_core
  public :: f_dst_exp_pdf, f_dst_exp_pdf_core
  public :: f_dst_exp_cdf, f_dst_exp_cdf_core
  public :: f_dst_exp_ppf, f_dst_exp_ppf_core
  public :: f_dst_chi2_pdf, f_dst_chi2_pdf_core
  public :: f_dst_chi2_cdf, f_dst_chi2_cdf_core
  public :: f_dst_chi2_ppf, f_dst_chi2_ppf_core
  public :: f_dst_f_pdf, f_dst_f_pdf_core
  public :: f_dst_f_cdf, f_dst_f_cdf_core
  public :: f_dst_f_ppf, f_dst_f_ppf_core
  public :: f_dst_gpd_pdf, f_dst_gpd_pdf_core
  public :: f_dst_gpd_cdf, f_dst_gpd_cdf_core
  public :: f_dst_gpd_ppf, f_dst_gpd_ppf_core

  ! additional core features; only exposed in fsml_dst
  public :: f_dst_gammai_core, f_dst_betai_core

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_norm_pdf(x, mu, sigma) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_norm_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in), optional :: mu      !! distribution location (mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)                       :: mu_w    !! final value of mu
  real(wp)                       :: sigma_w !! final value of sigma
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ---- compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_norm_pdf_core(x, mu_w, sigma_w)

end function f_dst_norm_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_pdf_core(x, mu, sigma) result(fx)

! ==== Description
!! Probability density function for normal distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: mu      !! distribution location (mean)
  real(wp), intent(in) :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)             :: z       !! z-score
  real(wp)             :: fx

! ==== Instructions

! ---- compute PDF

  ! compute z-score
  z = (x - mu) / sigma

  ! calculate probability/fx
  fx = (1.0_wp / (sigma * sqrt(2.0_wp * c_pi))) * exp( -0.5_wp * (z * z) )

end function f_dst_norm_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_norm_cdf(x, mu, sigma, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_norm_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: mu_w    !! final value of mu
  real(wp)                               :: sigma_w !! final value of sigma
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ---- compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_norm_cdf_core(x, mu_w, sigma_w, tail_w)

end function f_dst_norm_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_cdf_core(x, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for normal distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x       !! sample position
  real(wp)        , intent(in) :: mu      !! distribution location (mean)
  real(wp)        , intent(in) :: sigma   !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in) :: tail    !! tail options
  real(wp)                     :: z       !! z-score
  real(wp)                     :: p       !! returned probability integral

! ==== Instructions

! ---- compute CDF

  ! compute z-score
  z = (x - mu) / (sigma * sqrt(2.0_wp))

  ! compute integral (left tailed)
  p = 0.5_wp * (1.0_wp + erf(z))

  ! tail options
  ! NOTE: alternatively, compare z to 0.0 instead of x to mu
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. mu) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. mu) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function f_dst_norm_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_norm_ppf(p, mu, sigma) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_norm_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: mu               !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma            !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: mu_w             !! final value of mu
  real(wp)                          :: sigma_w          !! final value of sigma
  real(wp)                          :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_norm_ppf_core(p, mu_w, sigma_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_norm_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_norm_ppf_core(p, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function for normal distribution.

! ==== Declarations
  real(wp)   , intent(in) :: p                  !! probability between 0.0 - 1.0
  real(wp)   , intent(in) :: mu                 !! distribution location (mean)
  real(wp)   , intent(in) :: sigma              !! distribution dispersion/scale (standard deviation)
  integer(i4), parameter  :: i_max = c_bisect_i !! max. number of iterations
  real(wp)   , parameter  :: tol = c_bisect_tol !! p deviation tolerance
  real(wp)                :: a, b               !! section bounds for bisection algorithm
  real(wp)                :: x_mid, p_mid       !! x and p mid points in bisection algorithm
  integer(i4)             :: i                  !! for iteration
  real(wp)                :: x                  !! sample position

! ==== Instructions

! ---- compute PPF

  ! set initial section
  a = -20.0_wp
  b = 20.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_norm_cdf_core(x_mid, mu=0.0_wp, sigma=1.0_wp, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value and adjust for mu and sigma
        x = mu + sigma * x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_norm_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_t_pdf(x, df, mu, sigma) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_t_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp), intent(in), optional :: mu      !! distribution location (~mean)
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (~standard deviation)
  real(wp)                       :: mu_w    !! final value of mu
  real(wp)                       :: sigma_w !! final value of sigma
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ---- compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_t_pdf_core(x, df, mu_w, sigma_w)

end function f_dst_t_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_t_pdf_core(x, df, mu, sigma) result(fx)

! ==== Description
!! Probability density function for student t distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: df      !! degrees of freedom
  real(wp), intent(in) :: mu      !! distribution location (~mean)
  real(wp), intent(in) :: sigma   !! distribution dispersion/scale (~standard deviation)
  real(wp)             :: fx

! ==== Instructions

! ---- compute PDF

  ! calculate probability/fx
  fx = gamma((df + 1.0_wp) / 2.0_wp) / &
     & (sigma * sqrt(df * c_pi) * gamma(df / 2.0_wp)) * &
     & (1.0_wp + ( ( (x - mu) / sigma )**2 ) / df)** &
     & (-(df + 1.0_wp) / 2.0_wp)

end function f_dst_t_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_t_cdf(x, df, mu, sigma, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_t_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)        , intent(in)           :: x           !! sample position
  real(wp)        , intent(in)           :: df          !! degrees of freedom
  real(wp)        , intent(in), optional :: mu          !! distribution location (mean)
  real(wp)        , intent(in), optional :: sigma       !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in), optional :: tail        !! tail options
  real(wp)                               :: mu_w        !! final value of mu
  real(wp)                               :: sigma_w     !! final value of sigma
  character(len=16)                      :: tail_w      !! final tail option
  real(wp)                               :: p           !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

! ---- compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_t_cdf_core(x, df, mu_w, sigma_w, tail_w)

end function f_dst_t_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_t_cdf_core(x, df, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for student t distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x           !! sample position
  real(wp)        , intent(in) :: df          !! degrees of freedom
  real(wp)        , intent(in) :: mu          !! distribution location (mean)
  real(wp)        , intent(in) :: sigma       !! distribution dispersion/scale (standard deviation)
  character(len=*), intent(in) :: tail        !! tail options
  real(wp)                     :: z           !! z-score
  real(wp)                     :: xbeta, a, b !! parameters for beta function
  real(wp)                     :: p           !! returned probability integral

! ==== Instructions

! ---- compute CDF

  ! compute z-score
  z = (x - mu) / sigma

  ! shape parameters for beta function
  a = 0.5_wp * df
  b = 0.5_wp
  xbeta = df / (df + z**2)

  ! compute integral (left tailed)
  if (z .ge. 0.0_wp) then
    p = 1.0_wp - 0.5_wp * f_dst_betai_core(xbeta, a, b)
  else
    p = 0.5_wp * f_dst_betai_core(xbeta, a, b)
  endif

  ! tail options
  ! NOTE: alternatively, compare z to 0.0 instead of x to mu
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. mu) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. mu) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. mu) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function f_dst_t_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_t_ppf(p, df, mu, sigma) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_t_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)   , intent(in)           :: p       !! probability between 0.0 - 1.0
  real(wp)   , intent(in)           :: df      !! degrees of freedom
  real(wp)   , intent(in), optional :: mu      !! distribution location (mean)
  real(wp)   , intent(in), optional :: sigma   !! distribution dispersion/scale (standard deviation)
  real(wp)                          :: mu_w    !! final value of mu
  real(wp)                          :: sigma_w !! final value of sigma
  real(wp)                          :: x       !! sample position

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_t_ppf_core(p, df, mu_w, sigma_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_t_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_t_ppf_core(p, df, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function for t distribution.

! ==== Declarations
  real(wp)   , intent(in) :: p                  !! probability between 0.0 - 1.0
  real(wp)   , intent(in) :: df                 !! degrees of freedom
  real(wp)   , intent(in) :: mu                 !! distribution location (mean)
  real(wp)   , intent(in) :: sigma              !! distribution dispersion/scale (standard deviation)
  integer(i4), parameter  :: i_max = c_bisect_i !! max. number of iterations
  real(wp)   , parameter  :: tol = c_bisect_tol !! p deviation tolerance
  real(wp)                :: a, b               !! section bounds for bisection algorithm
  real(wp)                :: x_mid, p_mid       !! x and p mid points in bisection algorithm
  integer(i4)             :: i                  !! for iteration
  real(wp)                :: x                  !! sample position

! ==== Instructions

! ---- compute PPF

  ! set initial section bounds
  a = -20.0_wp
  b = 20.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_t_cdf_core(x_mid, df, mu=0.0_wp, sigma=1.0_wp, tail="left") - p
     ! check if difference is acceptable, update section if not
     if (abs(p_mid) .lt. tol) then
        ! pass final x value and adjust for mu and sigma
        x = mu + sigma * x_mid
        return
     elseif (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     endif
  enddo

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_t_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gamma_pdf(x, alpha, beta, loc) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_gamma_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

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

  ! assume alpha = 1, overwrite if specified
  alpha_w = 1.0_wp
  if (present(alpha)) alpha_w = alpha

  ! assume beta = 1, overwrite if specified
  beta_w = 1.0_wp
  if (present(beta)) beta_w = beta

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! check if alpha value is valid
  if (alpha_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

  ! check if beta value is valid
  if (beta_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ---- compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_gamma_pdf_core(x, alpha_w, beta_w, loc_w)

end function f_dst_gamma_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_pdf_core(x, alpha, beta, loc) result(fx)

! ==== Description
!! Probability density function for gamma distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: alpha   !! shape  parameter
  real(wp), intent(in) :: beta    !! scale parameter
  real(wp), intent(in) :: loc     !! location parameter
  real(wp)             :: z       !! shifted and scaled variable
  real(wp)             :: fx

! ==== Instructions

! ---- compute PDF

  ! calculate probability/fx
  if (x .le. 0.0_wp .or. alpha .le. 0.0_wp .or. beta .le. 0.0_wp) then
     fx = 0.0_wp
  else
     z = (x - loc) / beta
     fx = ( x ** (alpha - 1.0_wp) ) * exp( -z) / &
        & ( beta ** alpha * gamma(alpha))
  endif

end function f_dst_gamma_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gamma_cdf(x, alpha, beta, loc, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_gamma_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

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

  ! assume alpha = 1, overwrite if specified
  alpha_w = 1.0_wp
  if (present(alpha)) alpha_w = alpha

  ! assume beta = 1, overwrite if specified
  beta_w = 1.0_wp
  if (present(beta)) beta_w = beta

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if alpha value is valid
  if (alpha_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if beta value is valid
  if (beta_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ---- compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_gamma_cdf_core(x, alpha_w, beta_w, loc_w, tail_w)


end function f_dst_gamma_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_cdf_core(x, alpha, beta, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function for gamma distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x       !! sample position
  real(wp)        , intent(in) :: alpha   !! shape  parameter
  real(wp)        , intent(in) :: beta    !! scale parameter
  real(wp)        , intent(in) :: loc     !! location parameter
  character(len=*), intent(in) :: tail    !! tail options
  real(wp)                     :: z       !! standardised variable
  real(wp)                     :: p       !! returned probability integral

! ==== Instructions

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .le. loc .or. alpha .le. 0.0_wp .or. beta .le. 0.0_wp) then
     p = 0.0_wp
  else
     z = (x - loc) / beta
     p = f_dst_gammai_core(z, alpha)
  endif

  ! tail options
  select case(tail)
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

end function f_dst_gamma_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gamma_ppf(p, alpha, beta, loc) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_gamma_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)   , intent(in)           :: p                !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: alpha            !! shape  parameter
  real(wp)   , intent(in), optional :: beta             !! scale parameter
  real(wp)   , intent(in), optional :: loc              !! location parameter
  real(wp)                          :: alpha_w          !! final value for alpha
  real(wp)                          :: beta_w           !! final value for beta
  real(wp)                          :: loc_w            !! final value for loc
  real(wp)                          :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume alpha = 1, overwrite if specified
  alpha_w = 1.0_wp
  if (present(alpha)) alpha_w = alpha

  ! assume beta = 1, overwrite if specified
  beta_w = 1.0_wp
  if (present(beta)) beta_w = beta

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! check if alpha value is valid
  if (alpha_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if beta value is valid
  if (beta_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_gamma_ppf_core(p, alpha_w, beta_w, loc_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_gamma_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gamma_ppf_core(p, alpha, beta, loc) result(x)

! ==== Description
!! Percent point function/quantile function for gamma distribution.

! ==== Declarations
  real(wp)   , intent(in) :: p                  !! probability between 0.0 - 1.0
  real(wp)   , intent(in) :: alpha              !! shape  parameter
  real(wp)   , intent(in) :: beta               !! scale parameter
  real(wp)   , intent(in) :: loc                !! location parameter
  integer(i4), parameter  :: i_max = c_bisect_i !! max. number of iterations
  real(wp)   , parameter  :: tol = c_bisect_tol !! p deviation tolerance
  real(wp)                :: a, b               !! section bounds for bisection algorithm
  real(wp)                :: x_mid, p_mid       !! x and p mid points in bisection algorithm
  integer(i4)             :: i                  !! for iteration
  real(wp)                :: x                  !! sample position

! ==== Instructions

! ---- compute PPF

  ! set initial section
  a = loc
  b = loc + alpha * beta * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_gamma_cdf_core(x_mid, alpha=alpha&
          &, beta=beta, loc=loc, tail="left") - p
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

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_gamma_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_exp_pdf(x, lambda, loc) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_exp_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x        !! sample position
  real(wp), intent(in), optional :: loc      !! location parameter
  real(wp), intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                       :: loc_w    !! final value for loc
  real(wp)                       :: lambda_w !! final value for lambda
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume lambda = 1, overwrite if specified
  lambda_w = 1.0_wp
  if (present(lambda)) lambda_w = lambda

  ! check if lambda value is valid
  if (lambda_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ---- compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_exp_pdf_core(x, lambda_w, loc_w)

end function f_dst_exp_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_exp_pdf_core(x, lambda, loc) result(fx)

! ==== Description
!! Probability density function for exponential distribution.

! ==== Declarations
  real(wp), intent(in) :: x        !! sample position
  real(wp), intent(in) :: loc      !! location parameter
  real(wp), intent(in) :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)             :: fx

! ==== Instructions

! ---- compute PDF

  ! calculate probability/fx
  if (x .lt. loc) then
     fx = 0.0_wp
  elseif (lambda .lt. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = lambda * exp( -lambda * (x - loc) )
  endif

end function f_dst_exp_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_exp_cdf(x, lambda, loc, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_exp_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

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

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume lambda = 1, overwrite if specified
  lambda_w = 1.0_wp
  if (present(lambda)) lambda_w = lambda

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if lambda value is valid
  if (lambda_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ---- compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_exp_cdf_core(x, lambda_w, loc_w, tail_w)

end function f_dst_exp_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_exp_cdf_core(x, lambda, loc, tail) result(p)

! ==== Description
!! Cumulative distribution function for exponential distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x        !! sample position
  real(wp)        , intent(in) :: loc      !! location parameter
  real(wp)        , intent(in) :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  character(len=*), intent(in) :: tail     !! tail options
  real(wp)                     :: p        !! returned probability integral

! ==== Instructions

! ---- compute CDF

  ! compute integral (left tailed)
  if (x .lt. loc) then
     p = 0.0_wp
  elseif (lambda .lt. 0.0_wp) then
     p = 0.0_wp
  else
     p = 1.0_wp - exp( - lambda * (x - loc) )
  endif

  ! tail options
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
   end select

end function f_dst_exp_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_exp_ppf(p, lambda, loc) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_exp_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)   , intent(in)           :: p        !! probability between 0.0 - 1.0
  real(wp)   , intent(in), optional :: loc      !! location parameter
  real(wp)   , intent(in), optional :: lambda   !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  real(wp)                          :: loc_w    !! final value for loc
  real(wp)                          :: lambda_w !! final value for lambda
  real(wp)                          :: x        !! sample position

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume lambda = 1, overwrite if specified
  lambda_w = 1.0_wp
  if (present(lambda)) lambda_w = lambda

  ! check if lambda value is valid
  if (lambda_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if value invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_exp_ppf_core(p, lambda_w, loc_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_exp_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_exp_ppf_core(p, lambda, loc) result(x)

! ==== Description
!! Percent point function/quantile function for exponential distribution.

! ==== Declarations
  real(wp)   , intent(in) :: p                  !! probability between 0.0 - 1.0
  real(wp)   , intent(in) :: loc                !! location parameter
  real(wp)   , intent(in) :: lambda             !! lambda parameter, beta(scale) = 1/lambda = mu/mean
  integer(i4), parameter  :: i_max = c_bisect_i !! max. number of iterations
  real(wp)   , parameter  :: tol = c_bisect_tol !! p deviation tolerance
  real(wp)                :: a, b               !! section bounds for bisection algorithm
  real(wp)                :: x_mid, p_mid       !! x and p mid points in bisection algorithm
  integer(i4)             :: i                  !! for iteration
  real(wp)                :: x                  !! sample position

! ==== Instructions

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = loc
  b = loc - log(1.0_wp - p) / lambda * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     ! difference between passed p and new mid point p
     p_mid = f_dst_exp_cdf_core(x_mid, lambda=lambda, loc=loc, tail="left") - p
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

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_exp_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_chi2_pdf(x, df, loc, scale) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_chi2_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: df      !! degrees of freedom
  real(wp), intent(in), optional :: loc     !! location parameter
  real(wp), intent(in), optional :: scale   !! scale parameter
  real(wp)                       :: loc_w   !! final value for loc
  real(wp)                       :: scale_w !! final value for scale
  real(wp)                       :: fx      !! resulting PDF value

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ----compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_chi2_pdf_core(x, df, loc_w, scale_w)

end function f_dst_chi2_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_chi2_pdf_core(x, df, loc, scale) result(fx)

! ==== Description
!! Probability density function for the chi-squared distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: df      !! degrees of freedom
  real(wp), intent(in) :: loc     !! location parameter
  real(wp), intent(in) :: scale   !! scale parameter
  real(wp)             :: z       !! standardised variable
  real(wp)             :: fx      !! resulting PDF value

! ==== Instructions

! ----compute PDF
  z = (x - loc) / scale
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. scale .le. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = (1.0_wp / (2.0_wp ** (df / 2.0_wp) * &
        & gamma(df / 2.0_wp) * scale)) * &
        & z ** (df / 2.0_wp - 1.0_wp) * exp(-z / 2.0_wp)
  endif

end function f_dst_chi2_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_chi2_cdf(x, df, loc, scale, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_chi2_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

  ! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in)           :: df      !! degrees of freedom
  real(wp)        , intent(in), optional :: loc     !! location parameter
  real(wp)        , intent(in), optional :: scale   !! scale parameter
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: loc_w   !! final value for loc
  real(wp)                               :: scale_w !! final value for scale
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: p       !! resulting CDF value

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ----compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_chi2_cdf_core(x, df, loc_w, scale_w, tail_w)

end function f_dst_chi2_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_chi2_cdf_core(x, df, loc, scale, tail) result(p)

! ==== Description
!! Cumulative distribution function for the chi-squared distribution.

  ! ==== Declarations
  real(wp)        , intent(in) :: x       !! sample position
  real(wp)        , intent(in) :: df      !! degrees of freedom
  real(wp)        , intent(in) :: loc     !! location parameter
  real(wp)        , intent(in) :: scale   !! scale parameter
  character(len=*), intent(in) :: tail    !! tail options
  real(wp)                     :: z       !! standardised variable
  real(wp)                     :: p       !! resulting CDF value

! ==== Instructions

! ----compute CDF

  ! compute integral (left tailed)
  z = (x - loc) / scale
  if (z .le. 0.0_wp .or. df .le. 0.0_wp .or. scale .le. 0.0_wp) then
     p = 0.0_wp
  else
     p = f_dst_gammai_core(z / 2.0_wp, df / 2.0_wp)
  endif

  ! tail options
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (x .gt. loc) then
           p = 2.0_wp * (1.0_wp - p)
        elseif (x .le. loc) then
           p = 2.0_wp * p
        endif
     ! confidence interval
     case("confidence")
        if (x .gt. loc) then
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        elseif (x .le. loc) then
           p = 1.0_wp - 2.0_wp * p
        endif
   end select

end function f_dst_chi2_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_chi2_ppf(p, df, loc, scale) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_chi2_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: p                !! probability between 0.0 and 1.0
  real(wp), intent(in)           :: df               !! degrees of freedom
  real(wp), intent(in), optional :: loc              !! location parameter
  real(wp), intent(in), optional :: scale            !! scale parameter
  real(wp)                       :: loc_w            !! final value for loc
  real(wp)                       :: scale_w          !! final value for scale
  real(wp)                       :: x                !! sample position

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if degrees of freedom value is valid
  if (df .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_chi2_ppf_core(p, df, loc_w, scale_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_chi2_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_chi2_ppf_core(p, df, loc, scale) result(x)

! ==== Description
!! Percent point function/quantile functionfor the chi-squared distribution.

! ==== Declarations
  real(wp), intent(in)   :: p                  !! probability between 0.0 and 1.0
  real(wp), intent(in)   :: df                 !! degrees of freedom
  real(wp), intent(in)   :: loc                !! location parameter
  real(wp), intent(in)   :: scale              !! scale parameter
  integer(i4), parameter :: i_max = c_bisect_i !! max. number of iterations
  real(wp), parameter    :: tol = c_bisect_tol !! tolerance for convergence
  real(wp)               :: a, b               !! interval bounds
  real(wp)               :: x_mid, p_mid       !! midpoint and corresponding CDF value
  integer(i4)            :: i                  !! iteration counter
  real(wp)               :: x                  !! sample position

! ==== Instructions

! ---- compute PPF

  ! set initial section (min possible approaches 0)
  a = loc
  b = loc + df * 10.0_wp

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     p_mid = f_dst_chi2_cdf_core(x_mid, df, loc=loc&
          &, scale=scale, tail="left") - p
     if (abs(p_mid) .lt. tol) then
        x = x_mid
        return
     else if (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     end if
  end do

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_chi2_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_f_pdf(x, d1, d2, loc, scale) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_f_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: d1      !! numerator degrees of freedom
  real(wp), intent(in)           :: d2      !! denominator degrees of freedom
  real(wp), intent(in), optional :: loc     !! location parameter
  real(wp), intent(in), optional :: scale   !! scale parameter
  real(wp)                       :: loc_w   !! final location
  real(wp)                       :: scale_w !! final scale
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

  ! check if numerator degrees of freedom value is valid
  if (d1 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

  ! check if denominator degrees of freedom value is valid
  if (d2 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ----compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_f_pdf_core(x, d1, d2, loc_w, scale_w)

end function f_dst_f_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_f_pdf_core(x, d1, d2, loc, scale) result(fx)

! ==== Description
!! Probability density function for the F distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: d1      !! numerator degrees of freedom
  real(wp), intent(in) :: d2      !! denominator degrees of freedom
  real(wp), intent(in) :: loc     !! location parameter
  real(wp), intent(in) :: scale   !! scale parameter
  real(wp)             :: B       !! beta(0.5*d1, 0.5*d2)
  real(wp)             :: z       !! standardised variable
  real(wp)             :: fx

! ==== Instructions

! ----compute PDF

  ! get z score (standardise)
  z = (x - loc) / scale

  ! calculate beta(0.5*d1, 0.5*d2) from gamma functions
  B = gamma(0.5_wp * d1) * gamma(0.5_wp * d2) / &
    & gamma(0.5_wp * d1 + 0.5_wp * d2)

  ! calculate probability/fx
  if (z .le. 0.0_wp) then
     fx = 0.0_wp
  else
     fx = sqrt(( (d1 * z) ** d1 * d2 ** d2 ) / &
       & ( (d1 * z + d2)**(d1 + d2) )) / z / (scale * B)
  endif

end function f_dst_f_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_f_cdf(x, d1, d2, loc, scale, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_f_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in)           :: d1      !! numerator degrees of freedom
  real(wp)        , intent(in)           :: d2      !! denominator degrees of freedom
  real(wp)        , intent(in), optional :: loc     !! location parameter
  real(wp)        , intent(in), optional :: scale   !! scale parameter
  character(len=*), intent(in), optional :: tail    !! tail option
  real(wp)                               :: loc_w   !! final location
  real(wp)                               :: scale_w !! final scale
  character(len=16)                      :: tail_w  !! tail option final
  real(wp)                               :: p       !! output probability

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if numerator degrees of freedom value is valid
  if (d1 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if denominator degrees of freedom value is valid
  if (d2 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ----compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_f_cdf_core(x, d1, d2, loc_w, scale_w, tail_w)


end function f_dst_f_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_f_cdf_core(x, d1, d2, loc, scale, tail) result(p)

! ==== Description
!! Cumulative density function for the F distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x       !! sample position
  real(wp)        , intent(in) :: d1      !! numerator degrees of freedom
  real(wp)        , intent(in) :: d2      !! denominator degrees of freedom
  real(wp)        , intent(in) :: loc     !! location parameter
  real(wp)        , intent(in) :: scale   !! scale parameter
  character(len=*), intent(in) :: tail    !! tail option
  real(wp)                     :: z       !! standardised variable
  real(wp)                     :: xbeta   !! beta variable
  real(wp)                     :: p       !! output probability

! ==== Instructions

! ----compute CDF

  ! get z score (standardise)
  z = (x - loc) / scale

  ! z must be positive non-zero; return 0 if not
  if (z .le. 0.0_wp) then
    p = 0.0_wp
    return
  endif

  ! transform to beta domain
  xbeta = (d1 * z) / (d1 * z + d2)
  p = f_dst_betai_core(xbeta, 0.5_wp * d1, 0.5_wp * d2)

  ! tail options
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
     ! two-tailed
     case("two")
        if (z .lt. 1.0_wp) then
           p = 2.0_wp * p
        else
           p = 2.0_wp * (1.0_wp - p)
        endif
     ! confidence interval
     case("confidence")
        if (z .lt. 1.0_wp) then
           p = 1.0_wp - 2.0_wp * p
        else
           p = 1.0_wp - 2.0_wp * (1.0_wp - p)
        endif
  end select

end function f_dst_f_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_f_ppf(p, d1, d2, loc, scale) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_f_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: p                !! probability (0.0 < p < 1.0)
  real(wp), intent(in)           :: d1               !! numerator degrees of freedom
  real(wp), intent(in)           :: d2               !! denominator degrees of freedom
  real(wp), intent(in), optional :: loc              !! location parameter
  real(wp), intent(in), optional :: scale            !! scale parameter
  real(wp)                       :: loc_w            !! effective location
  real(wp)                       :: scale_w          !! effective scale
  real(wp)                       :: x                !! result: quantile at p

! ==== Instructions

! ---- handle input

  ! assume loc = 0, overwrite if specified
  loc_w = 0.0_wp
  if (present(loc)) loc_w = loc

  ! assume scale = 1, overwrite if specified
  scale_w = 1.0_wp
  if (present(scale)) scale_w = scale

  ! check if scale value is valid
  if (scale_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if numerator degrees of freedom value is valid
  if (d1 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if denominator degrees of freedom value is valid
  if (d2 .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_f_ppf_core(p, d1, d2, loc_w, scale_w)

  ! issue warning in case of suspicious result
  if (x .eq. c_sentinel_r) call s_err_warn(fsml_warning(1))

end function f_dst_f_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_f_ppf_core(p, d1, d2, loc, scale) result(x)

! ==== Description
!! Percent point function / quantile function for the F distribution.

! ==== Declarations
  real(wp), intent(in)   :: p                  !! probability (0.0 < p < 1.0)
  real(wp), intent(in)   :: d1                 !! numerator degrees of freedom
  real(wp), intent(in)   :: d2                 !! denominator degrees of freedom
  real(wp), intent(in)   :: loc                !! location parameter
  real(wp), intent(in)   :: scale              !! scale parameter
  integer(i4), parameter :: i_max = c_bisect_i !! max. number of iterations
  real(wp)   , parameter :: tol = c_bisect_tol !! tolerance for convergence
  real(wp)               :: a, b               !! search bounds
  real(wp)               :: x_mid, p_mid       !! midpoint and its CDF
  integer(i4)            :: i                  !! iteration counter
  real(wp)               :: x                  !! result: quantile at p

! ==== Instructions

! ---- compute PPF

  ! set initial section
  a = loc
  b = loc + scale * 100.0_wp  ! heuristically large upper bound

  ! iteratively refine with bisection method
  do i = 1, i_max
     x_mid = 0.5_wp * (a + b)
     p_mid = f_dst_f_cdf_core(x_mid, d1, d2&
          &, loc=loc, scale=scale, tail="left") - p
     if (abs(p_mid) .lt. tol) then
        x = x_mid
        return
     else if (p_mid .lt. 0.0_wp) then
        a = x_mid
     else
        b = x_mid
     end if
  end do

  ! if x not found in iterations, pass sentinel
  if (i .ge. i_max) x = c_sentinel_r

end function f_dst_f_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gpd_pdf(x, xi, mu, sigma) result(fx)

! ==== Description
!! Impure wrapper function for `f_dst_gpd_pdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp), intent(in)           :: x       !! sample position
  real(wp), intent(in)           :: xi      !! distribution shape parameter
  real(wp), intent(in), optional :: mu      !! distribution location
  real(wp), intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)                       :: mu_w    !! final value of mu
  real(wp)                       :: sigma_w !! final value of sigma
  real(wp)                       :: fx

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     fx = c_sentinel_r
     return
  endif

! ---- compute PDF

  ! call pure function to calculate probability/fx
  fx = f_dst_gpd_pdf_core(x, xi, mu_w, sigma_w)

end function f_dst_gpd_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gpd_pdf_core(x, xi, mu, sigma) result(fx)

! ==== Description
!! Probability density function for generalised pareto distribution.

! ==== Declarations
  real(wp), intent(in) :: x       !! sample position
  real(wp), intent(in) :: xi      !! distribution shape parameter
  real(wp), intent(in) :: mu      !! distribution location
  real(wp), intent(in) :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)             :: z       !! z-score
  real(wp)             :: fx

! ==== Instructions

! ---- compute PDF

  ! compute z-score
  z = (x - mu) / sigma

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

end function f_dst_gpd_pdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gpd_cdf(x, xi, mu, sigma, tail) result(p)

! ==== Description
!! Impure wrapper function for `f_dst_gpd_cdf_core`.
!! Handles optional arguments and invalid values for arguments.

! ==== Declarations
  real(wp)        , intent(in)           :: x       !! sample position
  real(wp)        , intent(in)           :: xi      !! distribution shape parameter
  real(wp)        , intent(in), optional :: mu      !! distribution location
  real(wp)        , intent(in), optional :: sigma   !! distribution dispersion/scale (must be positive)
  character(len=*), intent(in), optional :: tail    !! tail options
  real(wp)                               :: mu_w    !! final value of mu
  real(wp)                               :: sigma_w !! final value of sigma
  character(len=16)                      :: tail_w  !! final tail option
  real(wp)                               :: p       !! returned probability integral

! ==== Instructions

! ---- handle input

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! assume left-tailed, overwrite if specified
  tail_w = "left"
  if (present(tail)) tail_w = tail

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     p = c_sentinel_r
     return
  endif

  ! check if tail options are valid
  if (tail_w .ne. "left" .and. tail_w .ne. "right" .and. &
     &tail_w .ne. "two" .and. tail_w .ne. "confidence") then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(2))
     p = c_sentinel_r
     return
  endif

! ---- compute CDF

  ! call pure function to calculate probability integral
  p = f_dst_gpd_cdf_core(x, xi, mu_w, sigma_w, tail_w)


end function f_dst_gpd_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gpd_cdf_core(x, xi, mu, sigma, tail) result(p)

! ==== Description
!! Cumulative distribution function for generalised pareto distribution.

! ==== Declarations
  real(wp)        , intent(in) :: x       !! sample position
  real(wp)        , intent(in) :: xi      !! distribution shape parameter
  real(wp)        , intent(in) :: mu      !! distribution location
  real(wp)        , intent(in) :: sigma   !! distribution dispersion/scale (must be positive)
  character(len=*), intent(in) :: tail    !! tail options
  real(wp)                     :: z       !! z-score
  real(wp)                     :: p       !! returned probability integral

! ==== Instructions

! ---- compute CDF

  ! compute z-score
  z = (x - mu) / sigma

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
  select case(tail)
    ! left-tailed; P(z<x)
     case("left")
        p = p
     ! right-tailed; P(z>x)
     case("right")
        p = 1.0_wp - p
   end select

end function f_dst_gpd_cdf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
impure function f_dst_gpd_ppf(p, xi, mu, sigma) result(x)

! ==== Description
!! Impure wrapper function for `f_dst_gpd_ppf_core`.
!! Handles optional arguments and invalid values for arguments.

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

  ! assume location/mean = 0, overwrite if specified
  mu_w = 0.0_wp
  if (present(mu)) mu_w = mu

  ! assume sigma = 1, overwrite if specified
  sigma_w = 1.0_wp
  if (present(sigma)) sigma_w = sigma

  ! check if sigma value is valid
  if (sigma_w .le. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

  ! check if p value is valid
  if (p .gt. 1.0_wp .or. p .lt. 0.0_wp) then
     ! write error message and assign sentinel value if invalid
     call s_err_print(fsml_error(1))
     x = c_sentinel_r
     return
  endif

! ---- compute PPF

  ! call pure function to calculate x
  x = f_dst_gpd_ppf_core(p, xi, mu_w, sigma_w)

end function f_dst_gpd_ppf


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gpd_ppf_core(p, xi, mu, sigma) result(x)

! ==== Description
!! Percent point function/quantile function for generalised pareto distribution.

! ==== Declarations
  real(wp)        , intent(in) :: p       !! probability between 0.0 - 1.0
  real(wp)        , intent(in) :: mu      !! distribution location
  real(wp)        , intent(in) :: sigma   !! distribution dispersion/scale (must be positive)
  real(wp)        , intent(in) :: xi      !! distribution shape parameter
  real(wp)                     :: x       !! sample position

! ==== Instructions

! ---- compute PPF

  ! compute inverse cdf based on xi
  if (abs(xi) .lt. 1.0e-12_wp) then
     ! if xi is approximately zero, use exponential distribution
     x = mu - sigma * log(1.0_wp - p)
  else
     ! if xi is not zero, use general formula
     x = mu + (sigma / xi) * ( (1.0_wp - p) ** (-xi) - 1.0_wp )
  endif

end function f_dst_gpd_ppf_core


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_gammai_core(x, a) result(p)

! ==== Description
!! Incomplete gamma function (lower). Needed by gamma and chi-squared cdf.
!! Uses Fortran 2008+ intrinsics.

! ==== Declarations
  real(wp), intent(in)   :: x                     !! integration limit
  real(wp), intent(in)   :: a                     !! shape parameter
  real(wp)               :: p                     !! probability
  real(wp)               :: sum, term, lngamma_a
  real(wp)               :: ap, del, b, c, d, h
  real(wp)   , parameter :: eps = 1.0e-12_wp      !! convergence threshold
  real(wp)   , parameter :: fpmin = 1.0e-30_wp    !! small number to prevent division by zero
  integer(i4), parameter :: i_max = c_bisect_i    !! max. number of iterations
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

end function f_dst_gammai_core


! ==================================================================== !
! -------------------------------------------------------------------- !
elemental function f_dst_betai_core(x, a, b) result(betai)

! ==== Description
!! Computes the regularised incomplete beta function. beta_inc and beta_cf
!! algorithms are based on several public domain Fortran and C code,
!! Lentz's algorithm (1976), and modified to use 2008+ intrinsics.

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
    cf = betai_cf_core(x, a, b)
    betai = bt * cf / a
  else
    cf = betai_cf_core(1.0_wp - x, b, a)
    betai = 1.0_wp - bt * cf / b
  endif

  contains

  ! ------------------------------------------------------------------ !
  elemental function betai_cf_core(x, a, b) result(cf)

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
     integer(i4), parameter :: i_max = c_bisect_i !! max. number of iterations
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

  end function betai_cf_core

end function f_dst_betai_core

end module fsml_dst

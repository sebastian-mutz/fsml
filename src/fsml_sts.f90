module fsml_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for basic statistics.                                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic sample statistics.

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

! ==== Description
!! Computes arithmetic mean.
!! $$ \bar{x} = \frac{1}{n} \cdot \sum_{i=1}^{n} x_i $$
!! where \( n \) is the size of (or number of observations in) vector x,
!! \( x_i \) are individual elements in x, and
!! \( \bar{x} \) is the arithmetic mean of x.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: mean   !! arithmetic mean

! ==== Instructions
  mean = sum(x) / real(size(x), kind=wp)

end function f_sts_mean


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_var(x) result(var)

! ==== Description
!! Computes variance.
!! $$ \operatorname{var}(x) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x})^2 $$
!! where \( n \) is the size of (or number of observations in) vector x,
!! \( x_i \) are individual elements in x, and
!! \( \bar{x} \) is the arithmetic mean of x.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: var    !! variance

! ==== Instructions
  var = sum( ( x - f_sts_mean(x) ) * (x - f_sts_mean(x) ) ) / &
      & real(size(x), kind=wp)

end function f_sts_var


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_std(x) result(std)

! ==== Description
!! Computes standard deviation.
!! $$ \sigma = \sqrt{\operatorname{var}(x)} $$
!! where \( \operatorname{var}(x) \) is the variance of vector x.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp)             :: std    !! standard deviation

! ==== Instructions
  std = sqrt( f_sts_var(x) )

end function f_sts_std


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_cov(x,y) result(cov)

! ==== Description
!! Computes covariance. x and y must be the same size.
!! $$ \operatorname{cov}(x, y) = \frac{1}{n} \cdot \sum_{i=1}^{n} (x_i - \bar{x}) \cdot (y_i - \bar{y}) $$
!! where \( n \) is the size of (or number of observations in) vectors x and y,
!! \( x_i \) and \( y_i \) are individual elements in x and y, and
!! \( \bar{x} \) and \( \bar{y} \) are the arithmetic means of x and y.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: cov    !! covariance

! ==== Instructions
  cov = sum ( ( x - f_sts_mean(x) ) * (y - f_sts_mean(y) ) ) / &
      & real(size(x), kind=wp)

end function f_sts_cov


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_trend(x,y) result(trend)

! ==== Description
!! Computes regression coefficient/trend. x and y must be the same size.
!! $$ m = \frac{\operatorname{cov}(x, y)}{\operatorname{var}(x)} $$
!! where \( m \) is the slope of the regression line (linear trend),
!! \( \operatorname{cov}(x, y) \) is the covariance of x and y, and
!! \( \operatorname{var}(x) \) is the variance of x.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: trend  !! trend/regression slope

! ==== Instructions
  trend = f_sts_cov(x,y) / f_sts_var(x)

end function f_sts_trend


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_sts_corr(x,y) result(corr)

! ==== Description
!! Computes Pearson correlation coefficient. x and y must be the same size.
!! $$ \rho_{x,y} = \frac{\operatorname{cov}(x, y)}{\sigma_x \cdot \sigma_y} $$
!! where \( \rho_{x,y} \) is the Pearson correlation coefficient for vectors x and y,
!! \( \operatorname{cov}(x, y) \) is the covariance of x and y, and
!! \( \sigma_{x} \) and \( \sigma_{y} \) are the standard deviations of x and y.

! ==== Declarations
  real(wp), intent(in) :: x(:)   !! x vector (assumed size array)
  real(wp), intent(in) :: y(:)   !! y vector (assumed size array)
  real(wp)             :: corr   !! Pearson correlation coefficient

! ==== Instructions
  corr = f_sts_cov(x,y) / sqrt( f_sts_var(x) * f_sts_var(y) )

end function f_sts_corr

end module fsml_sts

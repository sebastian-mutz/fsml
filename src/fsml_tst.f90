module fsml_tst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for common statistical tests.                               |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for common statistical tests.

  ! load modules
  use :: fsml_typ

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_tst_t1s

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_tst_t1s(x, mu, t, p, alt)

! ==== Description
!! One sample t-test.
!! $$ t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$
!! Returns t or p value, depending on passed option.

! ==== Declarations
  real(wp)        , intent(in)           :: x(:) !! x vector (samples)
  real(wp)        , intent(in)           :: mu   !! population mean (null hypothesis expected value)
  character(len=*), intent(in), optional :: alt  !! options for alternative hypothesis
  real(wp)        , intent(out)          :: t    !! sample number
  real(wp)        , intent(out)          :: p    !! p-value
  real(wp)                               :: s    !! sample standard deviation

! ==== Instructions
  t = 0.0_wp
  p = 0.0_wp

end subroutine s_tst_t1s


end module fsml_tst

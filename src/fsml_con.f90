module fsml_con

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for constants.                                              |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for computational and mathematical constants.

  ! load modules
  use :: fsml_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: c_pi
  public :: sentinel_r
  public :: bisect_tol, bisect_i

  ! mathematical constants
  real(wp), parameter    :: c_pi = 3.1415926535897932384626433832795028841972_wp !! pi

  ! computational constants
  real(wp), parameter    :: sentinel_r = -999.0_wp  !! real sentinel value
  real(wp), parameter    :: bisect_tol = 1.0e-12_wp !! convergence tolerance for bisection method (dst)
  integer(i4), parameter :: bisect_i   = 200        !! max. number of iterations for bisection method (dst)

end module fsml_con

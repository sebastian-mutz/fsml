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
  public :: c_sentinel_r, c_sentinel_i
  public :: c_conv_tol, c_bisect_i, c_kmeans_i

  ! mathematical constants
  real(wp), parameter    :: c_pi = 3.1415926535897932384626433832795028841972_wp !! pi

  ! computational constants
  real(wp)   , parameter :: c_sentinel_r = -999.0_wp  !! real sentinel value
  integer(i4), parameter :: c_sentinel_i = 999        !! integer sentinel value
  real(wp)   , parameter :: c_conv_tol   = 1.0e-12_wp !! convergence tolerance (e.g., for bisection method (dst) and kmeans(nlp))
  integer(i4), parameter :: c_bisect_i   = 200        !! max. number of iterations for bisection method (dst)
  integer(i4), parameter :: c_kmeans_i   = 1000       !! max. number of iterations for kmeans clustering convergence

end module fsml_con

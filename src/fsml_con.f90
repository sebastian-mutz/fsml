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
!! Module for mathematical constants.

! load modules
  use :: fsml_typ

! basic options
  implicit none
  private

! declare public procedures
  public :: c_pi

! constants
  real(wp), parameter :: c_pi = 3.1415926535897932384626433832795028841972_wp

end module fsml_con

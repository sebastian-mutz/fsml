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
!! Constants module.

! load modules
  use :: fsml_typ

! basic options
  implicit none
  private

! declare public procedures
  public :: c_pi

! constants
  real(wp), parameter :: c_pi = 3.14159265358979323846_wp  !! pi (real64 can handle)

end module fsml_con

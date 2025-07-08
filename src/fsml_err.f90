module fsml_err
!
! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Errors handling.                                                   |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for everything related to error handling.

  ! load modules
  use :: fsml_ini
  use :: fsml_typ, only: fsml_typ_error

  ! basic options
  implicit none
  private

  ! declare public
  public :: fsml_error

! ==== Declarations

  ! error types
  type(fsml_typ_error) :: fsml_error(4)

! ==== Data

  ! errors
  data fsml_error(1)%msg/ "[fsml error] Argument (sigma) value not in valid&
                        & range (> 0.0). Returning sentinel."/
  data fsml_error(1)%sv/-1.0_wp/

  data fsml_error(2)%msg/ "[fsml error] Argument (p) value not in valid&
                        & range (0.0 - 1.0). Returning sentinel."/
  data fsml_error(2)%sv/-999.0_wp/

  data fsml_error(3)%msg/"[fsml error] Argument (tail) value not in list of&
                        & valid options. Returning sentinel."/
  data fsml_error(3)%sv/-1.0_wp/

  data fsml_error(4)%msg/ "[fsml error] Argument (degrees of freedom) value&
                        & not in valid range (> 1.0). Returning sentinel."/
  data fsml_error(4)%sv/-1.0_wp/


end module fsml_err

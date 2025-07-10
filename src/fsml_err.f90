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
  use :: fsml_utl
  use :: fsml_typ, only: fsml_typ_error

  ! basic options
  implicit none
  private

  ! declare public
  public :: fsml_error, s_err_print

! ==== Declarations

  ! error types
  type(fsml_typ_error) :: fsml_error(2)

! ==== Data

  ! errors
  data fsml_error(1)%msg/ "[fsml error] Argument value out of valid&
                        & range. Returning sentinel."/
  data fsml_error(1)%sv/-999.0_wp/

  data fsml_error(2)%msg/"[fsml error] Argument value not in list of&
                        & valid options. Returning sentinel."/
  data fsml_error(2)%sv/-999.0_wp/

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_err_print(error)

! ==== Description
!! Prints error message in specific format.

! ==== Declarations
  type(fsml_typ_error) :: error
  character(len=128)   :: fstring

! ==== Instructions
  fstring = trim(error%msg) // " (" // trim(f_utl_r2c(error%sv)) // ")"
!  fstring = txt_error // trim(error%msg) // txt_info // &
!          & " (" // trim(f_utl_r2c(error%sv)) // ")" // txt_reset
  write(std_e, '(A)') fstring

end subroutine

end module fsml_err

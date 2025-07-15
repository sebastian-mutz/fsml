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
  use :: fsml_con

  ! basic options
  implicit none
  private

  ! declare public
  public :: fsml_error, s_err_print
  public :: fsml_warning, s_err_warn

! ==== Declarations

  ! error messages
  character(len=128), parameter :: fsml_error(2) = [ character(len=128) :: &
                                  & "[fsml error] Argument value out of valid&
                                  & range. Returning sentinel.          ", &
                                  & "[fsml error] Argument value not in list&
                                  & of valid options. Returning sentinel." ]
  ! warning messages
  character(len=128), parameter :: fsml_warning(1) = [ character(len=128) :: &
                                  & "[fsml warning] Suspicious value returned.&
                                  & Convergence may not have been reached in&
                                  & bisection iterations." ]

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_err_print(error)

! ==== Description
!! Prints error message in specific format.

! ==== Declarations
  character(len=*), intent(in) :: error
  character(len=128)           :: fstring

! ==== Instructions
  fstring = trim(error) // " (" // trim(f_utl_r2c(sentinel_r)) // ")"
!  fstring = txt_error // trim(error) // txt_info // &
!          & " (" // trim(f_utl_r2c(sentinel_r)) // ")" // txt_reset
  write(std_e, '(A)') fstring

end subroutine s_err_print

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_err_warn(warning)

! ==== Description
!! Prints warning message in specific format.

! ==== Declarations
  character(len=*), intent(in) :: warning
  character(len=128)           :: fstring

! ==== Instructions
  fstring = warning
!  fstring = txt_warn // trim(warning) // txt_reset
  write(std_e, '(A)') fstring

end subroutine s_err_warn

end module fsml_err

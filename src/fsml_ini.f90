module fsml_ini
!
! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Initialisation module, includes kinds FSML works with.             |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Initialisation module, includes kinds FSML works with.

  ! load modules
  use :: iso_fortran_env, only: int32, int64, real32, real64, real128 &
                            &, input_unit, output_unit, error_unit
  ! load stdlib procedures for ansi terminal colours
!  use :: stdlib_ansi, only: fg_color_blue, fg_color_cyan, fg_color_magenta&
!                            &, ansi_code, style_bold, style_reset, operator(//)
  ! load stblib procedures for linear algebra
  use :: stdlib_linalg, only: eigh, chol

  ! basic options
  implicit none
  private

  ! declare public
  public :: hp, qp, dp, sp, wp, i4, i8
  public :: std_i, std_o, std_e, std_rw
  public :: eigh, chol
!  public :: txt_info, txt_error, txt_warn, txt_bold, txt_reset, ansi_code, operator(//)

! ==== Declarations

  ! define kinds (used consistently and explicitly in derived types and entire project)
  integer, parameter :: hp = selected_real_kind(p=33, r=4931) !! very high precision (for testing)
  integer, parameter :: qp = real128                          !! quadruple precision
  integer, parameter :: dp = real64                           !! double precision
  integer, parameter :: sp = real32                           !! single precision
  integer, parameter :: wp = dp                               !! working precision
  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64

  ! standard i/o
  integer, parameter :: std_i  = input_unit
  integer, parameter :: std_o  = output_unit
  integer, parameter :: std_e  = error_unit
  integer, parameter :: std_rw = 21

  ! text colours
!   type(ansi_code), parameter :: txt_info  = fg_color_blue
!   type(ansi_code), parameter :: txt_error = fg_color_magenta
!   type(ansi_code), parameter :: txt_warn  = fg_color_cyan
!   type(ansi_code), parameter :: txt_bold  = style_bold
!   type(ansi_code), parameter :: txt_reset = style_reset

end module fsml_ini

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

  ! basic options
  implicit none
  private

  ! declare public
  public :: hp, qp, dp, sp, wp, i4, i8
  public :: std_i, std_o, std_e, std_rw

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

end module fsml_ini

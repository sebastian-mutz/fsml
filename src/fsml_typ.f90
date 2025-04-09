module fsml_typ
!
! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Kinds and derived types for fsml.                                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! FSML kinds and derived types.

  ! load modules
  use :: iso_fortran_env, only: int32, int64, real32, real64&
                            &, input_unit, output_unit, error_unit

  ! basic options
  implicit none
  private

  ! declare public
  public :: dp, sp, wp, i4, i8
  public :: std_i, std_o, std_e, std_rw
  public :: fsml_typ_df

! ==== Declarations

  ! define kinds (used consistently and explicitly in derived types and entire project)
  integer, parameter :: sp = real32 !! single precision
  integer, parameter :: dp = real64 !! double precision
  integer, parameter :: wp = dp     !! working precision
  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64

  ! standard i/o
  integer, parameter :: std_i  = input_unit
  integer, parameter :: std_o  = output_unit
  integer, parameter :: std_e  = error_unit
  integer, parameter :: std_rw = 21

! ==== Definitions

  ! simple dataframe
  type :: fsml_typ_df
     !! Derived type for dataframe.
     integer(i4)                    :: id        !! ID/index for data frame
     character(len=64), allocatable :: nm        !! dataframe name
     real(wp)         , allocatable :: data(:,:) !! data
     integer(i4)      , allocatable :: row_id(:) !! ID/key/index for rows
     integer(i4)      , allocatable :: col_id(:) !! ID/key/index for columns
     character(len=64), allocatable :: row_nm(:) !! names/labels for rows
     character(len=64), allocatable :: col_nm(:) !! names/labels for columns
  end type fsml_typ_df

end module fsml_typ

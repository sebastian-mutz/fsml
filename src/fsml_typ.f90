module fsml_typ
!
! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Derived types for fsml.                                            |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for kinds and derived types; sets working precision.

  ! load modules
  use :: fsml_ini

  ! basic options
  implicit none
  private

  ! declare public
  public :: fsml_typ_df, fsml_typ_error
  public :: fsml_error

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

  ! error types and messages
  type :: fsml_typ_error
     !! Derived type for errors.
     integer(i4)       :: id  !! error ID
     character(len=64) :: msg !! error message
  end type fsml_typ_error

! ==== Data

  ! errors
  type(fsml_typ_error) :: fsml_error(2)

  ! errors
  data fsml_error(1)%id/1/
  data fsml_error(1)%msg/"argument not in valid value range"/
  data fsml_error(2)%id/2/
  data fsml_error(2)%msg/"unrecognised argument value"/


end module fsml_typ

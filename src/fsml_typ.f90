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
!! Module for kinds and derived types; sets working precision.

  ! load modules
  use :: fsml_ini

  ! basic options
  implicit none
  private

  ! declare public
  public :: fsml_typ_df

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

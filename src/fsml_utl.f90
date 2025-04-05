module fsml_utl

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Test application for fsml lib.                                     |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic statistics.

! load modules
  use :: fsml_typ
  use :: stdlib_io, only: loadtxt

! basic options
  implicit none
  private

! declare public procedures
  public :: s_utl_readcsv

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_utl_readcsv(infile, labelcol, labelrow)

! ==== Description
!! Uses fortran-gmt interface for creating a geographical map or heat map.

! ==== Declarations
  character(len=*) , intent(in)  :: infile     !! read csv file
!  type(fsml_typ_df), intent(out) :: df       !! dataframe
  logical, optional, intent(in)  :: labelcol   !! true if first column contains row labels
  logical, optional, intent(in)  :: labelrow   !! true if first row contains column lavels
  logical                        :: w_labelcol !! final value of labelcol
  logical                        :: w_labelrow !! inal value of labeleow
  character(len=1000)            :: line
  integer(i4)                    :: nrow, ncol, ios

! ==== Instructions

! ---- handle optional arguments

! is first column a list of labels
  if (present(labelcol)) then
     ! use passed value
     w_labelcol = labelcol
  else
     ! use defaul value
     w_labelcol = .true.
  endif

! is first row a list of labels
  if (present(labelrow)) then
     ! use passed value
     w_labelrow = labelrow
  else
     ! use defaul value
     w_labelrow = .true.
  endif

! --- get dims

! open existing file (check if there, stop in case of error)
  open(unit=std_rw, file=infile, status="old", action="read", iostat=ios)
  if (ios .ne. 0) then
    write(std_o,*) "Error opening file:", infile
    stop
  endif

! get number of rows
  nrow = 0
  do
    read(std_rw, '(A)', iostat=ios) line
    if (ios .ne. 0) exit
    if (len_trim(line) .gt. 0) nrow = nrow + 1
  enddo

! correct data row count if there is a label row
  if (w_labelrow) nrow = nrow -1
  write(std_o,*) "number of rows:", nrow

! get number of columns



end subroutine s_utl_readcsv





end module fsml_utl

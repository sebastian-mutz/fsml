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

! basic options
  implicit none
  private

! declare public procedures
  public :: s_utl_readcsv

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_utl_readcsv(infile, df, labelcol, labelrow, delimiter)

! ==== Description
!! Read CSV file directly into dataframe.

! ==== Declarations
  character(len=*) , intent(in)           :: infile      !! read csv file
  type(fsml_typ_df), intent(inout)        :: df          !! dataframe
  logical          , intent(in), optional :: labelcol    !! true if first column contains row labels
  logical          , intent(in), optional :: labelrow    !! true if first row contains column lavels
  character(len=1) , intent(in), optional :: delimiter   !! single char delimiter
  logical                                 :: w_labelcol  !! final value of labelcol
  logical                                 :: w_labelrow  !! final value of labeleow
  character(len=1)                        :: w_delimiter !! final value of delimiter
  character(len=1000)                     :: line        !! row will be read into line string
  character(len=64), allocatable          :: cells(:)    !! cells in row
  integer(i4)                             :: nrow        !! number of rows
  integer(i4)                             :: ncol        !! number of columns
  integer(i4)                             :: ios         !! io status
  integer(i4)                             :: i, j, k, p

! ==== Instructions

  write(std_o, '(a, a)') "Reading file: ", trim(infile)

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

! is first row a list of labels
  if (present(delimiter)) then
     ! use passed value
     w_delimiter = delimiter
  else
     ! use defaul value
     w_delimiter = ","
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
    read(std_rw, '(a)', iostat=ios) line
    if (ios .ne. 0) exit
    if (len_trim(line) .gt. 0) nrow = nrow + 1
  enddo

! get number of columns
  rewind(std_rw)
  read(std_rw, '(a)') line
  ncol = 1
  do i = 1, len_trim(line)
     if (line(i:i) .eq. w_delimiter) ncol = ncol + 1
  enddo

! summary
  write(std_o, '(a18,i5)') "number of rows:   ", nrow
  write(std_o, '(a18,i5)') "number of columns:", ncol
  write(std_o, *)

! --- allocate

! column and row ids and names
  allocate(df%row_id(nrow))
  allocate(df%row_nm(nrow))
  allocate(df%col_id(ncol))
  allocate(df%col_nm(ncol))
  allocate(cells(ncol))

! update row number for data only
  if (w_labelrow) then
     i = nrow - 1
  else
     i = nrow
  endif

! update column number for data only
  if (w_labelcol) then
     j = ncol - 1
  else
     j = ncol
  endif

! data matrix dimensions
  allocate(df%data(i,j))

! default values
  df%id        = 0
  df%nm        = "no_name"
  df%row_id(:) = 0
  df%row_nm(:) = "no_name"
  df%col_id(:) = 0
  df%col_nm(:) = "no_name"
  df%data(:,:) = 0.0_wp

! ---- read data into dataframe

  rewind(std_rw)
  ! loop through rows
  do i = 1, 1!nrow
     ! read entire line
     read(std_rw, '(a)') line
     ! split line up into cells
     call split_line(trim(line), w_delimiter, ncol, cells)
     print*, cells(:)
     ! pass cell contents to dataframe

  enddo

!   if (labelrow) then
!
!   if (labelcol) then
!      do j = 1, ncol
!         df%col_nm(j) = cell(j)
!      enddo
!      do i = 1, nrow
!         df%row_nm(i) = cell(1)
!         do j = 2, ncol
!            df%data(i,j) = cell(j)
!         enddo
!      enddo
!   else
!      do j = 1, ncol
!         df%col_nm(j) = cell(j)
!      enddo
!      do i = 2, nrow
!         do j = 1, ncol
!            df%data(i,j) = cell(j)
!         enddo
!      enddo
!   endif
!
!
!   else
!
!   if (labelcol) then
!      do i = 1, nrow
!         df%row_nm(i) = cell(1)
!         do j = 2, ncol
!            df%data(i,j) = cell(j)
!         enddo
!      enddo
!   else
!      do i = 1, nrow
!         do j = 1, ncol
!            df%data(i,j) = cell(j)
!         enddo
!      enddo
!   endif
!
!   endif

  close(std_rw)


  deallocate(cells)

  contains
     subroutine split_line(line, delimiter, ncol, cells)
     !! Splits passed line up into cells by delimiter.
     character(len=*) , intent(in)    :: line        !! row read into line string
     character(len=1) , intent(in)    :: delimiter   !! single char delimiter
     integer(i4)      , intent(in)    :: ncol        !! number of columns
     character(len=64), intent(inout) :: cells(ncol) !! cells between delimiter (in line)
     integer(i4)                      :: i, j, k, p
     ! reset initial cell position (p) and column number (j)
     p = 1
     j = 1
     ! go through line, left to right
     do k = p, len_trim(line)
        if (line(k:k) .eq. delimiter) then
           ! pass data
           cells(j) = line(p:k-1)
           ! update initial cell position and column number
           p = k + 1
           j = j + 1
           ! if last column, pass remaining data
           if (j .eq. ncol) cells(j) = line(p:len_trim(line))
        endif
     enddo
     end subroutine split_line

end subroutine s_utl_readcsv



end module fsml_utl

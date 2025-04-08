module fsml_utl

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Utilities module; includes io procedures.                          |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Utilities module; includes io procedures.

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
  logical          , optional, intent(in) :: labelcol    !! true if first column contains row labels
  logical          , optional, intent(in) :: labelrow    !! true if first row contains column lavels
  character(len=1) , optional, intent(in) :: delimiter   !! single char delimiter
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

  write(std_o, *)
  write(std_o, '(a, a)') "Reading file: ", trim(infile)

! ---- handle optional arguments

  ! pass labelcol if present, otherwise set to true
  if (present(labelcol)) then
     w_labelcol = labelcol
  else
     w_labelcol = .false.
  endif

  ! pass labelrow if present, otherwise set to true
  if (present(labelrow)) then
     w_labelrow = labelrow
  else
     w_labelrow = .false.
  endif

  ! pass delimiter if present, otherwise set to comma
  if (present(delimiter)) then
     w_delimiter = delimiter
  else
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

! --- allocate

  ! number of actual cells per row in file
  allocate(cells(ncol))

  ! update row number for data if labelrow
  i = merge(nrow - 1, nrow, w_labelrow)

  ! update column number for data if labelcol
  j = merge(ncol - 1, ncol, w_labelcol)

  ! column and row ids and names
  allocate(df%row_id(i))
  allocate(df%row_nm(i))
  allocate(df%col_id(j))
  allocate(df%col_nm(j))

  ! data matrix dimensions
  allocate(df%data(i,j))

  ! default values
  df%id = 0
  df%nm = "no_name"
  df%row_nm(:) = "no_name"
  df%col_nm(:) = "no_name"
  df%data(:,:) = 0.0_wp
  do k = 1, i
     df%row_id(k) = k
  enddo
  do k = 1, j
     df%col_id(k) = k
  enddo

! ---- read data into dataframe

  rewind(std_rw)

  ! read column names if there is a row of labels
  if (w_labelrow) then
     ! read entire line
     read(std_rw, '(a)') line
     ! split line up into cells
     call split_line(trim(line), w_delimiter, ncol, cells)
     ! pass cell values to dataframe
     if (w_labelcol) then
        df%col_nm(:) = cells(2:ncol)
     else
        df%col_nm(:) = cells(:)
     endif
  endif

  ! determine starting row for data
  k = merge(2, 1, w_labelrow)

  ! determine index correction for data
  p = merge(1, 0, w_labelrow)

  if (w_labelcol) then
     do i = k, nrow
        ! read entire line
        read(std_rw, '(a)') line
        ! split line up into cells
        call split_line(trim(line), w_delimiter, ncol, cells)
        ! pass cell values label column
        df%row_nm(i-p) = cells(1)
        ! pass cell values to data matrix
        do j = 2, ncol
           df%data(i-p,j-1) = s2r(cells(j))
        enddo
     enddo
  else
     do i = k, nrow
        ! read entire line
        read(std_rw, '(a)') line
        ! split line up into cells
        call split_line(trim(line), w_delimiter, ncol, cells)
        ! pass cell values to data matrix
        do j = 1, ncol
           df%data(i-p,j) = s2r(cells(j))
        enddo
     enddo
  endif

! ---- finish

  ! close file
  close(std_rw)

  ! deallocate
  deallocate(cells)

  ! summary
  write(std_o, '(a18,i5)') "number of rows:   ", nrow
  write(std_o, '(a18,i5)') "number of columns:", ncol
  write(std_o, *)

! ---- conatined procedures
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

  pure function s2r(s) result(r)
     !! Converts string to real.
     character(len=*), intent(in) :: s
     real(wp)                     :: r
     read(s, *) r
  end function

end subroutine s_utl_readcsv



end module fsml_utl

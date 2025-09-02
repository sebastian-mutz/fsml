program example_climate_clusters

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Research application example for climate science.                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini ! import wp; alternatively: iso_fortran_env, wp => real64

  implicit none

! ==== Description
!! The programme demonstrates the use of nonlinear procedures.

! ==== Declarations
  type(fsml_typ_df)        :: df        ! dataframe
  character(len=256)       :: infile    ! data file
  integer(i4)              :: nc        ! number of target clusters
  integer(i4)              :: nd        ! number of data points
  integer(i4)              :: nv        ! number of variables ("features")
  real(wp)   , allocatable :: x(:,:)    ! data
  real(wp)   , allocatable :: gm(:)     ! global mean per variable
  real(wp)   , allocatable :: cm(:,:)   ! cluster centroids
  real(wp)   , allocatable :: cov(:,:)  ! covariance matrix
  real(wp)   , allocatable :: sigma(:)  ! standard dev per variable
  integer(i4), allocatable :: cl(:)     ! cluster assignment per data point
  integer(i4), allocatable :: cc(:)     ! cluster sizes
  integer(i4)              :: i, j      ! for loops

! ==== Instructions

  ! set target cluster numbers
  nc = 5

  infile = "./example/research/data/Mutz_and_Ehlers_2019/&
  &Mutz_and_Ehlers_2019_e007_2_PI_t159l31_alterm.csv"

  call fsml_read_csv(infile, df, labelcol=.false., labelrow=.true., delimiter=",")

  nd = 100!size(df%row_id)
  nv = size(df%col_id)-2 ! remove lat and lon

  allocate( x(nd, nv)   )
  allocate( gm(nv)      )
  allocate( cm(nv, nc)  )
  allocate( cov(nv, nv) )
  allocate( sigma(nv)   )
  allocate( cl(nd)      )
  allocate( cc(nc)      )

  ! pass data (except lat and lon)
  x = df%data(101:200,3:nv+2)

  ! file summary
  write(6, *)
  write(6, '(a)') "file summary"
  write(6, '(a,i8)') "number of rows (data points):          ", nd
  write(6, '(a,i8)') "number of columns (variables/features):", nv
  write(6, *)

  ! print column id and names
  write(6, '(a)') "dataframe columns"
  do i = 1, nv
     write(6, '(i8, a2, a40)') df%col_id(i), ": ", df%col_nm(i+2)
  enddo
  write(6, *)

! ---- clustering

  ! conduct clustering
  call fsml_hclust(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)
  write(*,'(A)') "> clustering"
  print*
  write(*,'(A)')  "  cluster centroids:  "
  do j = 1, nc
    write(*,'(A,I2,A,I5,A,7F12.8)'), "  cluster", j, " (size: ", cc(j), ")", cm(1:nv, j)
  enddo
  print*
  write(*,'(A)'), '  global means (gm):'
  write(*,'(7F12.6)'), gm(1:nv)
  print*
  write(*,'(A)'), '  covariance matrix:'
  do i = 1, nv
     write(*,'(7F12.6)'), cov(i, 1:nv)
  enddo
  print*
  write(*,'(A)'), '  standard deviations (sigma):'
  write(*,'(7F12.6)'), sigma(1:nv)
  print*

end program

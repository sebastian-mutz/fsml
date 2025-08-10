program example_nlp

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for nonlinear procedures (nlp module).                    |
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
  integer(i4), parameter :: nd = 100 ! number of data points
  integer(i4), parameter :: nv = 5   ! number of variables ("features")
  integer(i4), parameter :: nc = 4   ! number of target clusters for hcluster
  real(wp) :: x(nd, nv)              ! data
  real(wp) :: gm(nv)                 ! global mean per variable
  real(wp) :: cm_h(nv, nc)           ! cluster centroids (hcluster)
  real(wp) :: cm_k(nv, nc)           ! cluster centroids (k-means)
  real(wp) :: cov(nv, nv)            ! covariance matrix (hcluster)
  real(wp) :: cov_k(nv, nv)
  real(wp) :: sigma(nv)              ! standard dev per variable
  integer(i4) :: cl(nd)              ! cluster assignment per data point
  integer(i4) :: cc(nc)              ! cluster sizes
  integer(i4) :: i, j                ! for loops

! ==== Instructions

  ! create data
  do i = 1, nd
    do j = 1, nv
      x(i,j) = sin(real(i*j, wp))  ! example data, replace with real input
    enddo
  enddo

! ---- hierarchical clustering

  ! conduct hierarchical clustering
  call fsml_hcluster(x, nd, nv, nc, gm, cm_h, cl, cc, cov, sigma)
  write(*,'(A)') "> hierarchical cluster analysis"
  print*
  write(*,'(A)')  "  cluster centroids:  "
  do j = 1, nc
    write(*,'(A,I2,A,I2,A,5F12.6)'), "  cluster", j, " (size: ", cc(j), ")", cm_h(1:nv, j)
  enddo
  print*
  write(*,'(A)'), '  global means (gm):'
  write(*,'(5F12.6)'), gm(1:nv)
  print*
  write(*,'(A)'), '  covariance matrix:'
  do i = 1, nv
     write(*,'(5F12.6)'), cov(i, 1:nv)
  enddo
  print*
  write(*,'(A)'), '  standard deviations (sigma):'
  write(*,'(5F12.6)'), sigma(1:nv)
  print*

! ---- k-means clustering

  ! conduct k-means clustering with clusters and covariance matrix from hcluster
  call fsml_kmeans(x, nd, nv, nc, cm_h, gm, cm_k, cl, cc, cov_k, sigma, cov_in=cov)
  write(*,'(A)') "> k-means cluster analysis"
  print*
  write(*,'(A)')  "  cluster centroids:  "
  do j = 1, nc
    write(*,'(A,I2,A,I2,A,5F12.6)'), "  cluster", j, " (size: ", cc(j), ")", cm_k(1:nv, j)
  enddo
  print*
  write(*,'(A)'), '  global means (gm):'
  write(*,'(5F12.6)'), gm(1:nv)
  print*
  write(*,'(A)'), '  covariance matrix:'
  do i = 1, nv
     write(*,'(5F12.6)'), cov(i, 1:nv)
  enddo
  print*
  write(*,'(A)'), '  standard deviations (sigma):'
  write(*,'(5F12.6)'), sigma(1:nv)
  print*




end program example_nlp

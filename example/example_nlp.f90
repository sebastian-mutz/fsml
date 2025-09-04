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
  call fsml_hclust(x, nd, nv, nc, gm, cm_h, cl, cc, cov, sigma)
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
  ! cluster centroids:
  ! cluster 1 (size: 10)    1.114982    1.280285    0.426739   -0.620822   -0.979638
  ! cluster 2 (size: 74)    0.046810   -0.284541   -0.044039   -0.064682    0.347077
  ! cluster 3 (size:  9)   -0.540211    0.961225   -1.147117    1.125295   -0.846290
  ! cluster 4 (size:  7)   -1.393120   -0.056833    1.330794    0.123859   -1.181532
  !
  ! global means (gm):
  !  -0.001272   -0.002720   -0.004636   -0.007745   -0.014948
  !
  ! covariance matrix:
  !   1.000000    0.000134    0.000297    0.000800    0.005803
  !   0.000134    1.000000    0.000932    0.006108   -0.003596
  !   0.000297    0.000932    1.000000   -0.003435   -0.002111
  !   0.000800    0.006108   -0.003435    1.000000   -0.001955
  !   0.005803   -0.003596   -0.002111   -0.001955    1.000000
  !
  ! standard deviations (sigma):
  !   0.712573    0.712679    0.714744    0.711827    0.711722


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
     write(*,'(5F12.6)'), cov_k(i, 1:nv)
  enddo
  print*
  write(*,'(A)'), '  standard deviations (sigma):'
  write(*,'(5F12.6)'), sigma(1:nv)
  print*
  ! cluster centroids:
  ! cluster 1 (size: 20)    0.897497    1.076016    0.566456    0.025259   -0.064641
  ! cluster 2 (size: 50)    0.268886   -0.806361   -0.226511   -0.142963    0.169572
  ! cluster 3 (size: 16)   -0.772046    1.132262   -0.927534    0.439275    0.033617
  ! cluster 4 (size: 14)   -1.360107    0.048683    1.059784   -0.027532   -0.551692
  !
  ! global means (gm):
  ! -0.001272   -0.002720   -0.004636   -0.007745   -0.014948
  !
  ! covariance matrix:
  !   1.000000    0.000134    0.000297    0.000800    0.005803
  !   0.000134    1.000000    0.000932    0.006108   -0.003596
  !   0.000297    0.000932    1.000000   -0.003435   -0.002111
  !   0.000800    0.006108   -0.003435    1.000000   -0.001955
  !   0.005803   -0.003596   -0.002111   -0.001955    1.000000
  !
  ! standard deviations (sigma):
  !   0.712573    0.712679    0.714744    0.711827    0.711722


! ---- hybrid hierarchical - k-means clustering

  ! conduct clustering (combines 2 steps above)
  call fsml_hkmeans(x, nd, nv, nc, gm, cm_k, cl, cc, cov_k, sigma)
  write(*,'(A)') "> hierarchical clustering with k-means refinements"
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
     write(*,'(5F12.6)'), cov_k(i, 1:nv)
  enddo
  print*
  write(*,'(A)'), '  standard deviations (sigma):'
  write(*,'(5F12.6)'), sigma(1:nv)
  print*
  ! cluster centroids:
  ! cluster 1 (size: 20)    0.897497    1.076016    0.566456    0.025259   -0.064641
  ! cluster 2 (size: 50)    0.268886   -0.806361   -0.226511   -0.142963    0.169572
  ! cluster 3 (size: 16)   -0.772046    1.132262   -0.927534    0.439275    0.033617
  ! cluster 4 (size: 14)   -1.360107    0.048683    1.059784   -0.027532   -0.551692
  !
  ! global means (gm):
  ! -0.001272   -0.002720   -0.004636   -0.007745   -0.014948
  !
  ! covariance matrix:
  !   1.000000    0.000134    0.000297    0.000800    0.005803
  !   0.000134    1.000000    0.000932    0.006108   -0.003596
  !   0.000297    0.000932    1.000000   -0.003435   -0.002111
  !   0.000800    0.006108   -0.003435    1.000000   -0.001955
  !   0.005803   -0.003596   -0.002111   -0.001955    1.000000
  !
  ! standard deviations (sigma):
  !   0.712573    0.712679    0.714744    0.711827    0.711722

end program example_nlp

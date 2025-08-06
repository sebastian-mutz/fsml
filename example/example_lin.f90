program example_lin

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for linear (algebra) statistical procedures (lin module). |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini ! import wp; alternatively: iso_fortran_env, wp => real64

  implicit none

! ==== Description
!! The programme demonstrates the use of procedures relying heavily on linear
!! algebra.

! ==== Declarations

! ---- general
  integer(i4)            :: i

! ---- PCA
  integer(i4), parameter :: m = 5, n = 3
  real(wp)   , parameter :: x1(m,n) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, &
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, &
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  &
                                     & ], shape=[5,3])
  real(wp)               ::  pc(m,n), eof(n,n), ew(n), eof_scaled(n,n), r2(n), wt(n)

! ---- LDA
  integer(i4), parameter :: nd = 5, nv = 3, nc = 2  !! no. of samples per class, 3 variables, 2 classes
  real(wp)   , parameter :: x2(nc,nv, nd) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! class 1, var 1
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! class 1, var 2
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp, & ! class 1, var 3
                                     & 1.1_wp, 1.3_wp, 1.5_wp, 1.4_wp, 1.6_wp, & ! class 2, var 1
                                     & 1.2_wp, 1.5_wp, 1.6_wp, 1.7_wp, 1.8_wp, & ! class 2, var 2
                                     & 1.0_wp, 1.2_wp, 1.3_wp, 1.1_wp, 1.4_wp  & ! class 2, var 3
                                     & ], shape=[nc,nv,nd])
  real(wp) :: score, g, mh
  real(wp) :: sa(n)


! ==== Instructions

! ---- Principal Component Analysis

  call fsml_pca(x1, m, n, pc=pc, ev=eof, ew=ew, r2=r2)
  write(*,'(A)') "> principal component analysis"
  print*
  write(*,'(A)')  "  principal components:  "
  write(*,'(3F10.5)')  (pc(i,:), i=1,m)
  print*
  write(*,'(A)')   "  eigenvectors:         "
  write(*,'(3F10.5)')  (eof(i,:), i=1,n)
  print*
  write(*,'(A)')   "  eigenvalues:          "
  write(*,'(3F10.5)')  ew
  print*
  write(*,'(A)')   "  explained variance:   "
  write(*,'(3F10.5)')  r2
  print*
  ! principal components:
  ! -0.54693  -0.07722   0.13896
  !  0.42315  -0.07015   0.27645
  ! -0.42586  -0.05419  -0.13456
  !  0.13778   0.43289  -0.06344
  !  0.41185  -0.23133  -0.21741
  !
  ! eigenvectors:
  !  0.28392   0.36310   0.88744
  !  0.47190   0.75276  -0.45897
  !  0.83469  -0.54909  -0.04238
  !
  ! eigenvalues:
  !  0.21204   0.06368   0.04128
  !
  ! explained variance:
  !  0.66888   0.20089   0.13023

! ---- Empirical Orthogonal Functions

  ! applying weights of 1 (detault = 1/n)
  wt = 1.0_wp

  ! call eof procedure with weights of 1 and opt=0 (covariance matrix) to makes it consistent with fmsl_pca and other commona pca implementations
  call fsml_eof(x1, m, n, pc=pc, eof=eof, ew=ew, opt=0, &
                 wt=wt, eof_scaled=eof_scaled, r2=r2)
  write(*,'(A)') "> empirical orthogonal function analysis"
  print*
  write(*,'(A)')  "  principal components:  "
  write(*,'(3F10.5)')  (pc(i,:), i=1,m)
  print*
  write(*,'(A)')   "  EOFs (eigenvectors):  "
  write(*,'(3F10.5)')  (eof(i,:), i=1,n)
  print*
  write(*,'(A)')   "  eigenvalues:          "
  write(*,'(3F10.5)')  ew
  print*
  write(*,'(A)')   "  explained variance:   "
  write(*,'(3F10.5)')  r2
  print*
  write(*,'(A)')   "  scaled EOFs:          "
  write(*,'(3F10.5)')  (eof_scaled(i,:), i=1,n)
  print*
  ! principal components:
  ! -0.54693  -0.07722   0.13896
  !  0.42315  -0.07015   0.27645
  ! -0.42586  -0.05419  -0.13456
  !  0.13778   0.43289  -0.06344
  !  0.41185  -0.23133  -0.21741
  !
  ! eigenvectors:
  !  0.28392   0.36310   0.88744
  !  0.47190   0.75276  -0.45897
  !  0.83469  -0.54909  -0.04238
  !
  ! eigenvalues:
  !  0.21204   0.06368   0.04128
  !
  ! explained variance:
  !  0.66888   0.20089   0.13023
  !
  ! scaled EOFs:
  !  0.13074   0.09163   0.18031
  !  0.21730   0.18996  -0.09325
  !  0.38435  -0.13856  -0.00861

! ---- 2-Class Multivariate Linear Discriminant Analysis

  call fsml_lda_2class(x2, nc, n, m, sa, g, score, mh)
  write(*,'(A)') "> linear discriminant analysis (2-class)"
  print*
  write(*,'(A,F10.5)') "  classification score: ", score
  write(*,'(A,F10.5)') "  Mahalanobis distance: ", mh
  write(*,'(A,F10.5)') "  discriminant value g: ", g
  print*
  write(*,'(A)') "  standardised coefficients:"
  write(*,'(3F10.5)') sa
  print*
  ! classification score:    0.80000
  ! Mahalanobis distance:    1.43391
  ! discriminant value g:   -2.58344
  !
  ! standardised coefficients:
  ! 0.38562   1.43199  -2.13615


end program example_lin

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
  integer(i4), parameter :: m = 5, n = 3
  integer(i4)            :: i
  real(wp)   , parameter :: x(m,n) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, &
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, &
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  &
                                     & ], shape=[5,3])
  real(wp)               ::  pc(m,n), eof(n,n), ew(n), eof_scaled(n,n), r2(n), wt(n)

! ==== Instructions

! ---- Principal Component Analysis

  call fsml_pca(x, m, n, pc=pc, ev=eof, ew=ew, r2=r2)
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
  call fsml_eof(x, m, n, pc=pc, eof=eof, ew=ew, opt=0, &
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


end program example_lin

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
  integer(i4)            :: i, j

! ---- data

  ! dims
  integer(i4), parameter :: nd = 5, nv = 3, nc = 2

  ! PCA/EOF
  real(wp)   , parameter :: x1(nd,nv) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, &
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, &
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  &
                                     & ], shape=[nd,nv])

  ! LDA
  real(wp)   , parameter :: x2(nd,nv,nc) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! class 1, var 1
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! class 1, var 2
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp, & ! class 1, var 3
                                     & 1.1_wp, 1.3_wp, 1.5_wp, 1.4_wp, 1.6_wp, & ! class 2, var 1
                                     & 1.2_wp, 1.5_wp, 1.6_wp, 1.7_wp, 1.8_wp, & ! class 2, var 2
                                     & 1.0_wp, 1.2_wp, 1.3_wp, 1.1_wp, 1.4_wp  & ! class 2, var 3
                                     & ], shape=[nd,nv,nc])

  ! OLS/Ridge
  real(wp), parameter :: x3(nd,nv) = reshape([ &
                                   & 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, & ! var 1
                                   & 2.0_wp, 1.0_wp, 4.0_wp, 3.0_wp, 2.0_wp, & ! var 2
                                   & 3.0_wp, 3.5_wp, 2.0_wp, 1.0_wp, 2.5_wp  & ! var 3
                                   & ], shape=[nd,nv])
  real(wp), parameter :: y(nd) = [10.0_wp, 12.0_wp, 13.0_wp, 14.0_wp, 15.0_wp]  ! target values

  ! Mahalanobis distance
  real(wp) :: dist         ! distance

! ---- results

  ! PCA
  real(wp) :: pc(nd,nv), eof(nv,nv), ew(nv), eof_scaled(nv,nv), r2(nv), wt(nv)

  ! LDA
  real(wp) :: score, g, mh, sa(nv)

  ! OLS and Ridge
  real(wp) :: b(nv), b0, yhat(nd), se(nv), covb(nv,nv), rsq

! ==== Instructions

! ---- Principal Component Analysis

  call fsml_pca(x1, nd, nv, pc=pc, ev=eof, ew=ew, r2=r2)
  write(*,'(A)') "> principal component analysis"
  print*
  write(*,'(A)')  "  principal components:  "
  write(*,'(3F10.5)')  (pc(i,:), i=1,nd)
  print*
  write(*,'(A)')   "  eigenvectors:         "
  write(*,'(3F10.5)')  (eof(i,:), i=1,nv)
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
  call fsml_eof(x1, nd, nv, pc=pc, eof=eof, ew=ew, opt=0, &
                 wt=wt, r2=r2, eof_scaled=eof_scaled)
  write(*,'(A)') "> empirical orthogonal function analysis"
  print*
  write(*,'(A)')  "  principal components:  "
  write(*,'(3F10.5)')  (pc(i,:), i=1,nd)
  print*
  write(*,'(A)')   "  EOFs (eigenvectors):  "
  write(*,'(3F10.5)')  (eof(i,:), i=1,nv)
  print*
  write(*,'(A)')   "  eigenvalues:          "
  write(*,'(3F10.5)')  ew
  print*
  write(*,'(A)')   "  explained variance:   "
  write(*,'(3F10.5)')  r2
  print*
  write(*,'(A)')   "  scaled EOFs:          "
  write(*,'(3F10.5)')  (eof_scaled(i,:), i=1,nv)
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

  call fsml_lda_2class(x2, nd, nv, nc, sa, g, score, mh)
  write(*,'(A)') "> linear discriminant analysis (2-class)"
  print*
  write(*,'(A,F10.5)') "  classification score: ", score
  write(*,'(A,F10.5)') "  Mahalanobis distance: ", mh
  write(*,'(A,F10.5)') "  discriminant value g: ", g
  print*
  write(*,'(A)') "  standardised coefficients:"
  write(*,'(3F10.5)') sa
  print*
  ! classification score:    1.00000
  ! Mahalanobis distance:    5.53247
  ! discriminant value g:   53.34966
  !
  ! standardised coefficients:
  ! 2.07805   9.14221   5.42794

! ---- Multiple Linear Regression (Ordinary Least Squares)

  call fsml_ols(x3, y, nd, nv, b0, b, rsq, yhat, se, covb)
  write(*,'(A)') "> multiple linear regression (OLS)"
  print*
  write(*,'(A,F10.5)') "  R²:        ", rsq
  write(*,'(A,F10.5)') "  Intercept: ", b0
  print*
  write(*,'(A)') "  Coefficients:"
  write(*,'(3F10.5)') b
  print*
  write(*,'(A)') "  Predicted values:"
  write(*,'(5F10.5)') yhat
  print*
  write(*,'(A)') "  Standard errors:"
  write(*,'(3F10.5)') se
  print*
  ! R²:           0.97361
  ! Intercept:    8.78516
  !
  ! Coefficients:
  ! 1.22266   0.05078   0.09375
  !
  ! Predicted values:
  ! 10.39062  11.60938  12.84375  13.92187  15.23437
  !
  ! Standard errors:
  ! 0.25240   0.43454   0.60515

! ---- Multiple Linear Regression (Ridge)

  ! call with lambda set to 0.2
  call fsml_ridge(x3, y, nd, nv, 0.2_wp, b0, b, rsq, yhat, se, covb)
  write(*,'(A)') "> multiple linear regression (Ridge)"
  print*
  write(*,'(A,F10.5)') "  R²:        ", rsq
  write(*,'(A,F10.5)') "  Intercept: ", b0
  print*
  write(*,'(A)') "  Coefficients:"
  write(*,'(3F10.5)') b
  print*
  write(*,'(A)') "  Predicted values:"
  write(*,'(5F10.5)') yhat
  print*
  write(*,'(A)') "  Standard errors:"
  write(*,'(3F10.5)') se
  print*
  ! R²:           0.97293
  ! Intercept:    9.11522
  !
  ! Coefficients:
  ! 1.18224   0.02590   0.03162
  !
  ! Predicted values:
  ! 10.44413  11.61628  12.82879  13.95351  15.15729
  !
  ! Standard errors:
  ! 0.23481   0.36909   0.49025

! ---- Mahalanobis distance

  ! call without passing covariance matrix
  dist = fsml_mahalanobis(x3(:,1), x3(:,2))
  write(*,'(A)') "> Mahalanobis distance"
  print*
  write(*,'(A,F10.5)') "  Distance:  ", dist
  print*
  ! Distance:     1.41421





end program example_lin

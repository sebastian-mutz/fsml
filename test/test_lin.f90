program test_lin

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for linear procedures (lin module)                           |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini

! ==== Declarations
  logical             :: status        !! status (test passed = true)
  real(wp), parameter :: tol = 1.0E-12 !! deviation tolerance; sp: 1.0E-6, dp: 1.0E-12, qp: 1.0E-30

! ==== Instructions

  print*
  print*, "LIN tests"

  print*, "> lin: testing PCA"
  status = test_pca(tol)
  call handle_status(status)

  print*, "> lin: testing EOF analysis"
  status = test_eof(tol)
  call handle_status(status)

  print*, "> lin: testing LDA"
  status = test_lda(tol)
  call handle_status(status)

  print*, "> lin: testing OLS regression"
  status = test_ols(tol)
  call handle_status(status)

  print*, "> lin: testing ridge regression"
  status = test_ridge(tol)
  call handle_status(status)

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine handle_status(status)

! ==== Description
!! Handles test status. (Message printing and stopping if needed)

! ==== Declarations
  logical, intent(in) :: status

! Instructions
  if (status) then
     print*, "  [✓] passed"
  else
     print*, "  [x] one or more failed"
     stop
  endif

end subroutine handle_status


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_pca(tol) result(status)

! ==== Description
!! Tests for PCA.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i, j
  integer(i4), parameter :: nd = 5, nv = 3
  real(wp)   , parameter :: x1(nd,nv) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, &
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, &
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  &
                                     & ], shape=[nd,nv])
  real(wp) :: res_pc(nd,nv), res_eof(nv,nv), res_ew(nv), res_r2(nv)
  real(wp), parameter :: ans_pc(nd,nv) = reshape([  &
                     & -0.54692539750109925_wp,     &
                     & 0.42315070723686443_wp,      &
                     & -0.42586030649496504_wp,     &
                     & 0.13778484019820619_wp,      &
                     & 0.41185015656099278_wp,      &
                     & -7.7216628990462305E-002_wp, &
                     & -7.0146451550147676E-002_wp, &
                     & -5.4193391785707812E-002_wp, &
                     & 0.43288620251018844_wp,      &
                     & -0.23132973018387198_wp,     &
                     & 0.13896115203759024_wp,      &
                     & 0.27644701897445756_wp,      &
                     & -0.13455881850988216_wp,     &
                     & -6.3441890639094542E-002_wp, &
                     & -0.21740746186307131_wp      &
                     & ], shape=[nd,nv])
  real(wp), parameter :: ans_eof(nv,nv) = reshape([ &
                     & 0.28391840658445738_wp,      &
                     & 0.47190119706668648_wp,      &
                     & 0.83468532909688298_wp,      &
                     & 0.36310044104714717_wp,      &
                     & 0.75276383524334656_wp,      &
                     & -0.54909441634485356_wp,     &
                     & 0.88743924192809365_wp,      &
                     & -0.45897262288371443_wp,     &
                     & -4.2375975851105398E-002_wp  &
                     & ], shape=[nv,nv])
  real(wp), parameter :: ans_ew(nv) = &
                     & [0.21203653144063464_wp     , 6.3680941140400460E-002_wp, &
                     & 4.1282527418964869E-002_wp]
  real(wp), parameter :: ans_r2(nv) = &
                     & [0.66888495722597685_wp     , 0.20088624965426016_wp, &
                     & 0.13022879311976301_wp]

! ==== Instructions

  status = .true.

  ! pca
  call fsml_pca(x1, nd, nv, pc=res_pc, ev=res_eof, ew=res_ew, r2=res_r2)
  do i = 1,nv
     do j = 1,nd
        if (abs( res_pc(j,i) - ans_pc(j,i) ) .gt. tol) status = .false.
     enddo
  enddo
  do i = 1, nv
     do j = 1, nv
        if (abs( res_eof(j,i) - ans_eof(j,i) ) .gt. tol) status = .false.
     enddo
     if (abs( res_ew(i) - ans_ew(i) ) .gt. tol) status = .false.
     if (abs( res_r2(i) - ans_r2(i) ) .gt. tol) status = .false.
  enddo

end function test_pca


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_eof(tol) result(status)

! ==== Description
!! Tests for EOF analysis.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i, j
  integer(i4), parameter :: nd = 5, nv = 3
  real(wp)   , parameter :: x1(nd,nv) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, &
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, &
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  &
                                     & ], shape=[nd,nv])
  real(wp) :: res_pc(nd,nv), res_eof(nv,nv), res_ew(nv)
  real(wp) :: res_eof_scaled(nv,nv), res_r2(nv), wt(nv)
  real(wp), parameter :: ans_pc(nd,nv) = reshape([  &
                     & -0.54692539750109925_wp,     &
                     & 0.42315070723686443_wp,      &
                     & -0.42586030649496504_wp,     &
                     & 0.13778484019820619_wp,      &
                     & 0.41185015656099278_wp,      &
                     & -7.7216628990462305E-002_wp, &
                     & -7.0146451550147676E-002_wp, &
                     & -5.4193391785707812E-002_wp, &
                     & 0.43288620251018844_wp,      &
                     & -0.23132973018387198_wp,     &
                     & 0.13896115203759024_wp,      &
                     & 0.27644701897445756_wp,      &
                     & -0.13455881850988216_wp,     &
                     & -6.3441890639094542E-002_wp, &
                     & -0.21740746186307131_wp      &
                     & ], shape=[nd,nv])
  real(wp), parameter :: ans_eof(nv,nv) = reshape([ &
                     & 0.28391840658445738_wp,      &
                     & 0.47190119706668648_wp,      &
                     & 0.83468532909688298_wp,      &
                     & 0.36310044104714717_wp,      &
                     & 0.75276383524334656_wp,      &
                     & -0.54909441634485356_wp,     &
                     & 0.88743924192809365_wp,      &
                     & -0.45897262288371443_wp,     &
                     & -4.2375975851105398E-002_wp  &
                     & ], shape=[nv,nv])
  real(wp), parameter :: ans_ew(nv) = &
                     & [0.21203653144063464_wp     , 6.3680941140400460E-002_wp, &
                     & 4.1282527418964869E-002_wp]
  real(wp), parameter :: ans_r2(nv) = &
                     & [0.66888495722597685_wp     , 0.20088624965426016_wp, &
                     & 0.13022879311976301_wp]

! ==== Instructions

  status = .true.

  ! applying weights of 1 (detault = 1/n)
  wt = 1.0_wp

  ! eof
  call fsml_eof(x1, nd, nv, pc=res_pc, eof=res_eof, ew=res_ew, opt=0, &
                 wt=wt, r2=res_r2, eof_scaled=res_eof_scaled)
  do i = 1,nv
     do j = 1,nd
        if (abs( res_pc(j,i) - ans_pc(j,i) ) .gt. tol) status = .false.
     enddo
  enddo
  do i = 1, nv
     do j = 1, nv
        if (abs( res_eof(j,i) - ans_eof(j,i) ) .gt. tol) status = .false.
     enddo
     if (abs( res_ew(i) - ans_ew(i) ) .gt. tol) status = .false.
     if (abs( res_r2(i) - ans_r2(i) ) .gt. tol) status = .false.
  enddo

end function test_eof




! ==================================================================== !
! -------------------------------------------------------------------- !
function test_lda(tol) result(status)

! ==== Description
!! Tests for LDA.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i
  integer(i4), parameter :: nd = 5, nv = 3, nc = 2
  real(wp)   , parameter :: x1(nd,nv,nc) = reshape([ &
                         & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! class 1, var 1
                         & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! class 1, var 2
                         & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp, & ! class 1, var 3
                         & 1.1_wp, 1.3_wp, 1.5_wp, 1.4_wp, 1.6_wp, & ! class 2, var 1
                         & 1.2_wp, 1.5_wp, 1.6_wp, 1.7_wp, 1.8_wp, & ! class 2, var 2
                         & 1.0_wp, 1.2_wp, 1.3_wp, 1.1_wp, 1.4_wp  & ! class 2, var 3
                         & ], shape=[nd,nv,nc])
  real(wp) :: res_sa(nv), res_g, res_score, res_mh
  real(wp), parameter :: ans_score  = 1.0000000000000000_wp
  real(wp), parameter :: ans_g      = 42.679725844362636_wp
  real(wp), parameter :: ans_sa(nv) = &
          & [1.7523677799718411_wp, 7.7093900223233245_wp, 4.5772446141297660_wp]

! ==== Instructions

  status = .true.

  ! lda
  call fsml_lda_2class(x1, nd, nv, nc, res_sa, res_g, res_score, res_mh)
  if (abs( res_score - ans_score ) .gt. tol) status = .false.
  if (abs( res_g - ans_g ) .gt. tol) status = .false.
  do i = 1, nv
     if (abs( res_sa(i) - ans_sa(i) ) .gt. tol) status = .false.
  enddo

end function test_lda




! ==================================================================== !
! -------------------------------------------------------------------- !
function test_ols(tol) result(status)

! ==== Description
!! Tests for OLS regression.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i
  integer(i4), parameter :: nd = 5, nv = 3
  real(wp), parameter :: x1(nd,nv) = reshape([ &
                                   & 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, & ! var 1
                                   & 2.0_wp, 1.0_wp, 4.0_wp, 3.0_wp, 2.0_wp, & ! var 2
                                   & 3.0_wp, 3.5_wp, 2.0_wp, 1.0_wp, 2.5_wp  & ! var 3
                                   & ], shape=[nd,nv])
  real(wp), parameter :: y(nd) = [10.0_wp, 12.0_wp, 13.0_wp, 14.0_wp, 15.0_wp]  ! target values

  ! OLS and Ridge
  real(wp) :: res_b(nv), res_b0, res_yhat(nd)
  real(wp) :: res_se(nv), res_covb(nv,nv), res_rsq
  real(wp), parameter :: ans_b(nv) = &
          & [1.2226562500000071_wp, 5.0781249999943157E-002_wp, &
          &  9.3750000000028422E-002_wp]
  real(wp), parameter :: ans_b0 = 8.7851562500000000_wp
  real(wp), parameter :: ans_yhat(nd) = &
          & [10.390624999999979_wp, 11.609375000000057_wp, &
          &  12.843749999999851_wp, 13.921874999999886_wp, &
          &  15.234374999999993_wp]
  real(wp), parameter :: ans_se(nv) = &
          & [0.25239937467421891_wp, 0.43454288010325609_wp, &
          &  0.60515364784489889_wp]
  real(wp), parameter :: ans_rsq = 0.97360641891891897_wp

! ==== Instructions

  status = .true.

  ! reg
  call fsml_ols(x1, y, nd, nv, res_b0, res_b, res_rsq, res_yhat, res_se, res_covb)
  if (abs( res_b0 - ans_b0 ) .gt. tol) status = .false.
  if (abs( res_rsq - ans_rsq ) .gt. tol) status = .false.
  do i = 1, nv
     if (abs( res_b(i) - ans_b(i) ) .gt. tol) status = .false.
     if (abs( res_se(i) - ans_se(i) ) .gt. tol) status = .false.
  enddo
  do i = 1, nd
     if (abs( res_yhat(i) - ans_yhat(i) ) .gt. tol) status = .false.
  enddo

end function test_ols




! ==================================================================== !
! -------------------------------------------------------------------- !
function test_ridge(tol) result(status)

! ==== Description
!! Tests for ridge regression.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i
  integer(i4), parameter :: nd = 5, nv = 3
  real(wp), parameter :: x1(nd,nv) = reshape([ &
                                   & 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, & ! var 1
                                   & 2.0_wp, 1.0_wp, 4.0_wp, 3.0_wp, 2.0_wp, & ! var 2
                                   & 3.0_wp, 3.5_wp, 2.0_wp, 1.0_wp, 2.5_wp  & ! var 3
                                   & ], shape=[nd,nv])
  real(wp), parameter :: y(nd) = [10.0_wp, 12.0_wp, 13.0_wp, 14.0_wp, 15.0_wp]  ! target values

  ! OLS and Ridge
  real(wp) :: res_b(nv), res_b0, res_yhat(nd)
  real(wp) :: res_se(nv), res_covb(nv,nv), res_rsq
  real(wp), parameter :: ans_b(nv) = &
          & [1.2226562500000071_wp, 5.0781249999943157E-002_wp, &
          &  9.3750000000028422E-002_wp]
  real(wp), parameter :: ans_b0 = 8.7851562500000000_wp
  real(wp), parameter :: ans_yhat(nd) = &
          & [10.390624999999979_wp, 11.609375000000057_wp, &
          &  12.843749999999851_wp, 13.921874999999886_wp, &
          &  15.234374999999993_wp]
  real(wp), parameter :: ans_se(nv) = &
          & [0.25239937467421891_wp, 0.43454288010325609_wp, &
          &  0.60515364784489889_wp]
  real(wp), parameter :: ans_rsq = 0.97360641891891897_wp

! ==== Instructions

  status = .true.

  ! reg
  call fsml_ridge(x1, y, nd, nv, 0.0_wp, res_b0, res_b, res_rsq, res_yhat, res_se, res_covb)
  if (abs( res_b0 - ans_b0 ) .gt. tol) status = .false.
  if (abs( res_rsq - ans_rsq ) .gt. tol) status = .false.
  do i = 1, nv
     if (abs( res_b(i) - ans_b(i) ) .gt. tol) status = .false.
     if (abs( res_se(i) - ans_se(i) ) .gt. tol) status = .false.
  enddo
  do i = 1, nd
     if (abs( res_yhat(i) - ans_yhat(i) ) .gt. tol) status = .false.
  enddo

end function test_ridge


end program

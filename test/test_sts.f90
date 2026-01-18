program test_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for basic (descriptive) statistics functions (sts module)    |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini
  use :: fsml_utl

! ==== Declarations
  logical             :: status        !! status (test passed = true)
  real(wp), parameter :: tol = 1.0E-12 !! deviation tolerance; sp: 1.0E-6, dp: 1.0E-12, qp: 1.0E-30

! ==== Instructions

  print*
  print*, "STS tests"

  print*, "> sts: testing centre metrics (mean, median)"
  status = test_centre(tol)
  call handle_status(status)

  print*, "> sts: testing dispersion metrics (variance, std. dev.)"
  status = test_dispersion(tol)
  call handle_status(status)

  print*, "> sts: testing basic statistical relationships (cov., trend, Pearson- and Spearman corr.)"
  status = test_relationship(tol)
  call handle_status(status)

  print*, "> sts: testing results of invalid arguments (and print messages)"
  status = test_invalid()
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
function test_centre(tol) result(status)

! ==== Description
!! Tests normal centre estimates (mean and median).

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 3      !! number of tests
  real(wp)   , parameter :: x(5,n) = reshape([ &
                          9781.5_wp, 7892.2_wp, 9880.6_wp, 9891.5_wp, 10872.96_wp, &
                          0.02_wp, -0.25_wp, 0.5_wp, 0.75_wp, 0.99_wp, &
                          -1.02_wp, -0.25_wp, -0.5_wp, -0.75_wp, -0.99_wp], shape=[5,n])
  real(wp), parameter :: ans(n*2) =                            &
                                  [9663.7520000000004_wp,      &
                                   0.40199999999999997_wp,     &
                                   -0.70199999999999996_wp,    &
                                   9880.6000000000004_wp,      &
                                   0.50000000000000000_wp,     &
                                   -0.75000000000000000_wp  ]

! ==== Instructions
  j=0
  status = .true.
  do i = 1, n
     j=j+1
     res = fsml_mean( x(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_median( x(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo

end function test_centre




! ==================================================================== !
! -------------------------------------------------------------------- !
function test_dispersion(tol) result(status)

! ==== Description
!! Tests dispersion metrics (variance, standard deviation).

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 3      !! number of tests
  real(wp)   , parameter :: x(5,n) = reshape([ &
                          9781.5_wp, 7892.2_wp, 9880.6_wp, 9891.5_wp, 10872.96_wp, &
                          0.02_wp, -0.25_wp, 0.5_wp, 0.75_wp, 0.99_wp, &
                          -1.02_wp, -0.25_wp, -0.5_wp, -0.75_wp, -0.99_wp], shape=[5,n])
  real(wp), parameter :: ans(n*2) =                            &
                                  [942667.45481599984_wp,      &
                                   0.20949599999999999_wp,     &
                                   8.6296000000000012E-002_wp, &
                                   970.91063173497071_wp,      &
                                   0.45770733007020981_wp,     &
                                   0.29376180827330162_wp]

! ==== Instructions
  j=0
  status = .true.
  do i = 1, n
     j=j+1
     res = fsml_var( x(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_std( x(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo

end function test_dispersion



! ==================================================================== !
! -------------------------------------------------------------------- !
function test_relationship(tol) result(status)

! ==== Description
!! Tests relationship metrics (covariance, ols linear trend,
!! Pearson correlation, Spearman correlation).

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 3      !! number of tests
  real(wp)   , parameter :: x(5,n) = reshape([ &
                          9781.5_wp, 7892.2_wp, 9880.6_wp, 9891.5_wp, 10872.96_wp, &
                          0.02_wp, -0.25_wp, 0.5_wp, 0.75_wp, 0.99_wp, &
                          -1.02_wp, -0.25_wp, -0.5_wp, -0.75_wp, -0.99_wp], shape=[5,n])
  real(wp)   , parameter :: y(5,n) = reshape([ &
                          9471.5_wp, 8219.2_wp, 8305.6_wp, 9753.5_wp, 10123.31_wp, &
                          0.12_wp, -0.12_wp, 0.3_wp, 0.95_wp, 0.69_wp, &
                          -2.02_wp, 0.25_wp, -0.1_wp, 0.5_wp, 0.3_wp], shape=[5,n])

  real(wp), parameter :: ans(n*4) =                             &
                                  [704522.57221999927_wp,       &
                                   0.19952999999999999_wp,      &
                                   0.15619000000000002_wp,      &
                                   0.59789701542842866_wp,      &
                                   0.76194294879138502_wp,      &
                                   1.4479466023917680_wp,       &
                                   0.75052027185654679_wp,      &
                                   0.90367127810428349_wp,      &
                                   0.46060469163085388_wp,      &
                                   0.90000000000000002_wp,      &
                                   0.90000000000000002_wp,      &
                                   0.20000000000000001_wp]

! ==== Instructions
  j=0
  status = .true.
  do i = 1, n
     j=j+1
     res = fsml_cov( x(:,i), y(:,i), ddf=1.0_wp)
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_trend( x(:,i), y(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_pcc( x(:,i), y(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_scc( x(:,i), y(:,i) )
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo

end function test_relationship



! ==================================================================== !
! -------------------------------------------------------------------- !
function test_invalid() result(status)

! ==== Description
!! Tests on a subset of functions if invalid arguments return NaN values.

! ==== Declarations
  logical             :: status     !! status (test passed = true)
  real(wp)            :: res        !! result
  real(wp), parameter :: x1(5) = [0.02_wp, -0.25_wp, 0.5_wp, 0.75_wp, 0.99_wp]
  real(wp), parameter :: x2(1) = [0.02_wp]

! ==== Instructions

  status = .true.

  ! mean of array wth invalid dimensions
  res = fsml_mean(x2)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! median of array wth invalid dimensions
  res = fsml_median(x2)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! variance with infalid ddf
  res = fsml_var(x1, 0.5_wp)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! covariance of unequal length vectors
  res = fsml_cov(x1, x1(2:5))
  if (.not. f_utl_is_nan(res)) status = .false.

end function test_invalid


end program

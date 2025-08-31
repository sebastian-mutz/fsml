program test_dst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for statistical hypothesis test module (tst module)          |
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

  print*, "> testing all t-tests"
  status = test_ttests(tol)
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
     print*, fg_color_cyan // "  passed" // style_reset
  else
     print*,  fg_color_magenta // "  > error: one or more failed" // style_reset
     stop
  endif

end subroutine handle_status


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_ttests(tol) result(status)

! ==== Description
!! Tests for all types of t-tests.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i, j
  integer(i4), parameter :: n1 = 10, n2 = 13
  real(wp)   , parameter :: mu = 1.25_wp
  real(wp)   , parameter :: x1(n1) =        &
                                  [1.10_wp, &
                                   1.02_wp, &
                                   1.12_wp, &
                                   2.01_wp, &
                                   1.92_wp, &
                                   1.01_wp, &
                                   1.10_wp, &
                                   1.26_wp, &
                                   1.51_wp, &
                                   1.01_wp]
  real(wp)   , parameter :: x2(n2) =        &
                                  [2.10_wp, &
                                   1.02_wp, &
                                   1.22_wp, &
                                   1.05_wp, &
                                   0.95_wp, &
                                   1.02_wp, &
                                   2.00_wp, &
                                   1.05_wp, &
                                   1.12_wp, &
                                   1.20_wp, &
                                   1.12_wp, &
                                   1.01_wp, &
                                   1.12_wp]
  real(wp)            :: res_t, res_p, res_df
  real(wp), parameter :: ans_t(4) =                        &
                                  [0.46723663923080166_wp, &
                                   0.15844543599004465_wp, &
                                   0.48615905360287642_wp, &
                                   0.48504551341216551_wp  ]

  real(wp), parameter :: ans_p(4) =                        &
                                  [0.65143607004113613_wp, &
                                   0.87760409339719425_wp, &
                                   0.63188896449497833_wp, &
                                   0.63309210403053418_wp  ]

  real(wp), parameter :: ans_df(4) =                        &
                                  [9.0000000000000000_wp,   &
                                   9.0000000000000000_wp,   &
                                   21.000000000000000_wp,   &
                                   19.342285575583787_wp    ]

! ==== Instructions

  j=0

  ! 1-sample t-test
  j=j+1
  call fsml_ttest_1sample(x1, mu, res_t, res_df, res_p, h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.
  if (abs( res_df - ans_df(j) ) .gt. tol) status = .false.

  ! paired t-test
  j=j+1
  call fsml_ttest_paired(x1, x2(:n1), res_t, res_df, res_p, h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.
  if (abs( res_df - ans_df(j) ) .gt. tol) status = .false.

  ! 2-sample pooled t-test (assume equal variances)
  j=j+1
  call fsml_ttest_2sample(x1, x2, res_t, res_df, res_p, eq_var=.true., h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.
  if (abs( res_df - ans_df(j) ) .gt. tol) status = .false.

  ! 2-sample t-test for unequal variances (Welch's t-test)
  j=j+1
  call fsml_ttest_2sample(x1, x2, res_t, res_df, res_p, eq_var=.false., h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.
  if (abs( res_df - ans_df(j) ) .gt. tol) status = .false.

end function test_ttests







end program

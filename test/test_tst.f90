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

  print*
  print*, "TST tests"

  print*, "> tst: testing all t-tests"
  status = test_ttests(tol)
  call handle_status(status)

  print*, "> tst: testing all Mann–Whitney U / Wilcoxon tests"
  status = test_ranktests(tol)
  call handle_status(status)

  print*, "> tst: testing ANOVA"
  status = test_anova(tol)
  call handle_status(status)

  print*, "> tst: testing Kruskal Wallis H test"
  status = test_kw(tol)
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


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_ranktests(tol) result(status)

! ==== Description
!! Tests for rank based tests (non-parametric equivalents of the t-tests).

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
  real(wp)            :: res_t, res_p
  real(wp), parameter :: ans_t(3) =                       &
                                  [27.000000000000000_wp, &
                                   21.000000000000000_wp, &
                                   59.500000000000000_wp  ]

  real(wp), parameter :: ans_p(3) =                        &
                                  [0.95935363404016349_wp, &
                                   0.85895492273748242_wp, &
                                   0.73303167364837196_wp  ]

! ==== Instructions

  j=0

  ! 1-sample Wilcoxon signed rank test
  j=j+1
  call fsml_signedrank_1sample(x1, mu, res_t, res_p, h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.

  ! 2-sample (paired sample) Wilcoxon signed rank test
  j=j+1
  call fsml_signedrank_paired(x1, x2(:n1), res_t, res_p, h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.

  ! 2-sample Wilcoxon Mann-Whitney U rank sum test
  j=j+1
  call fsml_ranksum(x1, x2, res_t, res_p, h1="two")
  if (abs( res_t - ans_t(j) ) .gt. tol) status = .false.
  if (abs( res_p - ans_p(j) ) .gt. tol) status = .false.

end function test_ranktests


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_anova(tol) result(status)

! ==== Description
!! Tests for ANOVA.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i
  real(wp)   , parameter :: x2d(5,3) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! Group 1
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! Group 2
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  & ! Group 3
                                     & ], shape=[5,3])
  real(wp)            :: res_t, res_p, res_df1, res_df2
  real(wp), parameter :: ans_t = 4.5867507886435366_wp
  real(wp), parameter :: ans_p = 3.3138386787152863E-002_wp
  real(wp), parameter :: ans_df1 = 2.0000000000000000_wp
  real(wp), parameter :: ans_df2 = 12.000000000000000_wp

! ==== Instructions

  ! ANOVA
  call fsml_anova_1way(x2d, res_t, res_df1, res_df2, res_p)
  if (abs( res_t - ans_t ) .gt. tol) status = .false.
  if (abs( res_p - ans_p ) .gt. tol) status = .false.
  if (abs( res_df1 - ans_df1 ) .gt. tol) status = .false.
  if (abs( res_df2 - ans_df2 ) .gt. tol) status = .false.

end function test_anova


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_kw(tol) result(status)

! ==== Description
!! Tests for Kruskal Wallis H test.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i
  real(wp)   , parameter :: x2d(5,3) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! Group 1
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! Group 2
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  & ! Group 3
                                     & ], shape=[5,3])
  real(wp)            :: res_t, res_p, res_df
  real(wp), parameter :: ans_t = 5.9850000000000065_wp
  real(wp), parameter :: ans_p = 5.0161875149147495E-002_wp
  real(wp), parameter :: ans_df = 2.0000000000000000_wp

! ==== Instructions

  ! Kruskal Wallis H test
  call fsml_kruskalwallis(x2d, res_t, res_df, res_p)
  if (abs( res_t - ans_t ) .gt. tol) status = .false.
  if (abs( res_p - ans_p ) .gt. tol) status = .false.
  if (abs( res_df - ans_df ) .gt. tol) status = .false.

end function test_kw




end program

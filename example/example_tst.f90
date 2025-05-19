program example_tst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for statistical tests (tst module).                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_typ

  implicit none

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
  real(wp)               :: t, p, df

  print*

  ! ---- 1 Sample T-Test

  call fsml_ttest_1sample(x1, mu, t, df, p, h1="two")
  write(*,'(A)') "> 1 sample student's t-test"
  write(*,'(A20,10F10.4)') " samples:", x1
  write(*,'(A20,F10.4)') " test statistic:", t
  write(*,'(A20,F10.4)') " degrees of freedom:", df
  write(*,'(A20,F10.4)') " p-value:", p
  print*

  ! ---- 2 Sample T-Test

  call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.false., h1="two")
  write(*,'(A)') "> 2 sample Welch's t-test"
  write(*,'(A20,10F10.4)') " samples 1:", x1
  write(*,'(A20,13F10.4)') " samples 2:", x2
  write(*,'(A20,F10.4)') " test statistic:", t
  write(*,'(A20,F10.4)') " degrees of freedom:", df
  write(*,'(A20,F10.4)') " p-value:", p
  print*


end program example_tst

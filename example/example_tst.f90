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

  integer(i4), parameter :: n    = 10
  real(wp)   , parameter :: mu   = 1.25_wp
  real(wp)   , parameter :: x(n) =          &
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
  real(wp)               :: t, p
  integer(i4)            :: df

  print*

  ! ---- Student's T-Test (1 Sample)

  call fsml_ttest_1sample(x, mu, t, df, p, "two")
  write(*,'(A)') "> 1 sample student's t-test"
  write(*,'(A20,10F10.4)') "  samples:", x
  write(*,'(A20,F10.4)') "  test statistic:", t
  write(*,'(A20,I10)') "  degrees of freedom:", df
  write(*,'(A20,F10.4)') "  p-value:", p

end program example_tst

program example_sts

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for basic statistics (sts module).                        |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini ! import wp; alternatively: iso_fortran_env, wp => real64

  implicit none

! ==== Description
!! The programme demonstrates the use of several statistics procedures based on
!! the same few data vectors.

! ==== Declarations
  integer(i4), parameter :: n1 = 10
  real(wp)   , parameter :: x1(n1) =        &
                                  [1.10_wp, &
                                   1.22_wp, &
                                   0.11_wp, &
                                   2.31_wp, &
                                   1.72_wp, &
                                   2.21_wp, &
                                   3.10_wp, &
                                   2.27_wp, &
                                   1.21_wp, &
                                   1.01_wp]
  real(wp)   , parameter :: x2(n1) =        &
                                  [2.10_wp, &
                                   1.10_wp, &
                                   0.23_wp, &
                                   1.25_wp, &
                                   1.95_wp, &
                                   1.52_wp, &
                                   1.00_wp, &
                                   1.05_wp, &
                                   2.32_wp, &
                                   1.40_wp]

! ==== Instructions

  ! mean of vector x1
  print*, "mean (x1): ", fsml_mean(x1)
  ! 1.6260000000000001

  ! median of vector x1
  print*, "median (x1): ", fsml_median(x1)
  ! 1.4700000000000000

  ! variance of vector x1
  print*, "variance (x1): ", fsml_var(x1)
  ! 0.66974399999999989

  ! sample variance of vector x1
  print*, "sample variance (x1): ", fsml_var(x1, ddf=1.0_wp)
  ! 0.74415999999999993

  ! standard deviation of vector x1
  print*, "standard deviation (x1): ", fsml_std(x1)
  ! 0.81837888535812064

  ! sample standard deviation of vector x1
  print*, "sample standard deviation (x1): ", fsml_std(x1, ddf=1.0_wp)
  ! 0.86264708890716135

  ! covariance of x1 and x2
  print*, "covariance (x1, x2): ", fsml_cov(x1, x2)
  ! 3.4878000000000006E-002

  ! sample covariance of x1 and x2
  print*, "sample covariance (x1, x2): ", fsml_cov(x1, x2, ddf=1.0_wp)
  ! 3.8753333333333341E-002

  ! linear regression slope for x1 and x2
  print*, "trend (x1, x2): ", fsml_trend(x1, x2)
  ! 5.2076614348168869E-002

  ! Pearson correlation coefficient for x1 and x2
  print*, "Pearson R (x1, x2): ", fsml_pcc(x1, x2)
  ! 7.2912607800353288E-002

  ! Spearman rank correlation coefficient for x1 and x2
  print*, "Spearman R (x1, x2): ", fsml_scc(x1, x2)
  ! -0.19999999999999998

end program example_sts

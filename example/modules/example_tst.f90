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
  use :: fsml_ini ! import wp; alternatively: iso_fortran_env, wp => real64

  implicit none

! ==== Description
!! The programme demonstrates the use of all statistical tests based on
!! the same few data vectors.

! ==== Declarations
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
  real(wp)   , parameter :: x2d(5,3) = reshape([ &
                                     & 2.1_wp, 2.5_wp, 1.9_wp, 2.3_wp, 2.0_wp, & ! Group 1
                                     & 2.4_wp, 2.8_wp, 2.6_wp, 3.2_wp, 2.9_wp, & ! Group 2
                                     & 2.0_wp, 2.8_wp, 2.1_wp, 2.3_wp, 2.9_wp  & ! Group 3
                                     & ], shape=[5,3])
  real(wp)               :: t, p, df, df1, df2

! ==== Instructions

! ---- Parametric Tests

  print*
  write(*,'(A)') " Parametric Tests"
  write(*,'(A)') " ================"
  print*

  ! 1-sample t-test
  call fsml_ttest_1sample(x1, mu, t, df, p, h1="two")
  write(*,'(A)') "> 1-sample student's t-test"
  write(*,'(A,10F10.4)') "  samples:            ", x1
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (t):     0.4672
  ! degrees of freedom:     9.0000
  ! p-value:                0.6514


  ! paired t-test
  call fsml_ttest_paired(x1, x2(:n1), t, df, p, h1="two")
  write(*,'(A)') "> 2-sample paired t-test"
  write(*,'(A,10F10.4)') "  samples 1:          ", x1
  write(*,'(A,10F10.4)') "  samples 2:          ", x2(:n1)
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (t):     0.1584
  ! degrees of freedom:     9.0000
  ! p-value:                0.8776


  ! 2-sample pooled t-test (assume equal variances)
  call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.true., h1="two")
  write(*,'(A)') "> 2-sample pooled t-test"
  write(*,'(A,10F10.4)') "  samples 1:          ", x1
  write(*,'(A,13F10.4)') "  samples 2:          ", x2
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (t):     0.4862
  ! degrees of freedom:    21.0000
  ! p-value:                0.6319


  ! 2-sample t-test for unequal variances (Welch's t-test)
  call fsml_ttest_2sample(x1, x2, t, df, p, eq_var=.false., h1="two")
  write(*,'(A)') "> 2-sample Welch's t-test"
  write(*,'(A,10F10.4)') "  samples 1:          ", x1
  write(*,'(A,13F10.4)') "  samples 2:          ", x2
  write(*,'(A,F10.4)')   "  test statistic (t): ", t
  write(*,'(A,F10.4)')   "  degrees of freedom: ", df
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (t):     0.4850
  ! degrees of freedom:    19.3423
  ! p-value:                0.6331

  ! 1-way ANOVA
  call fsml_anova_1way(x2d, t, df1, df2, p)
  write(*,'(A)') "> 1-way ANOVA"
  write(*,'(A,5F10.4)') "  group 1:            ", x2d(:, 1)
  write(*,'(A,5F10.4)') "  group 2:            ", x2d(:, 2)
  write(*,'(A,5F10.4)') "  group 3:            ", x2d(:, 3)
  write(*,'(A,F10.4)')  "  test statistic (F): ", t
  write(*,'(A,F10.4)')  "  df between:         ", df1
  write(*,'(A,F10.4)')  "  df within:          ", df2
  write(*,'(A,F10.4)')  "  p-value:            ", p
  ! test statistic (F):     4.5868
  ! df between:             2.0000
  ! df within:             12.0000
  ! p-value:                0.0331


! ---- Non-Parametric Tests

  print*
  write(*,'(A)') " Non-Parametric Tests"
  write(*,'(A)') " ===================="
  print*

  ! 1-sample Wilcoxon signed rank test
  call fsml_signedrank_1sample(x1, mu, t, p, h1="two")
  write(*,'(A)') "> 1-sample Wilcoxon signed rank test"
  write(*,'(A,10F10.4)') "  samples:            ", x1
  write(*,'(A,F10.4)')   "  test statistic (w): ", t
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (w):    27.0000
  ! p-value:                0.9594

  ! 2-sample (paired sample) Wilcoxon signed rank test
  call fsml_signedrank_paired(x1, x2(:n1), t, p, h1="two")
  write(*,'(A)') "> paired sample Wilcoxon signed rank test"
  write(*,'(A,10F10.4)') "  samples 1:          ", x1
  write(*,'(A,10F10.4)') "  samples 2:          ", x2(:n1)
  write(*,'(A,F10.4)')   "  test statistic (w): ", t
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (w):    21.0000
  ! p-value:                0.8590


  ! 2-sample Wilcoxon Mann-Whitney U rank sum test
  call fsml_ranksum(x1, x2, t, p, h1="two")
  write(*,'(A)') "> Mann-Whitney U rank sum test"
  write(*,'(A,10F10.4)') "  samples 1:          ", x1
  write(*,'(A,13F10.4)') "  samples 2:          ", x2
  write(*,'(A,F10.4)')   "  test statistic (U): ", t
  write(*,'(A,F10.4)')   "  p-value:            ", p
  print*
  ! test statistic (U):    59.5000
  ! p-value:                0.7330


  ! Kruskal Wallis H test
  call fsml_kruskalwallis(x2d, t, df, p)
  write(*,'(A)') "> Kruskal Wallis H test"
  write(*,'(A,5F10.4)') "  group 1:            ", x2d(:, 1)
  write(*,'(A,5F10.4)') "  group 2:            ", x2d(:, 2)
  write(*,'(A,5F10.4)') "  group 3:            ", x2d(:, 3)
  write(*,'(A,F10.4)')  "  test statistic (H): ", t
  write(*,'(A,F10.4)')  "  degrees of freedom  ", df
  write(*,'(A,F10.4)')  "  p-value:            ", p
  ! test statistic (H):     5.9850
  ! degrees of freedom      2.0000
  ! p-value:                0.0502

end program example_tst

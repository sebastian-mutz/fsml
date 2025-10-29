program main

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Test application for fsml lib.                                     |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml_ini !iso_fortran_env, wp => real64
  use :: fsml
  use :: fsml_sts, only: corr => f_sts_pcc_core
  use :: fsml_dst, only: norm_pdf_elemental => f_dst_norm_pdf_core

  implicit none

  type(fsml_typ_df)  :: df
  character(len=128) :: infile
  integer            :: i
  real(wp)           :: r
  real(wp)           :: x1(10,5), x2(10,5), pcc(5), x(10)

!   infile = "./example/research/data/Mutz_et_al_2021/DMC_Mutz2021_Antofagasta.csv"
!
!   call fsml_read_csv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")
!
!   ! print column id and names
!   print*, "dataframe columns"
!   do i = 1, size(df%col_id)
!      print*, df%col_id(i), df%col_nm(i)
!   enddo
!   print*
!
!   ! mean of first variable (msl - mean sea level pressure)
!   print*, "mean: ", fsml_mean(df%data(:,1))
!
!   ! variance of second variable (t2m - 2m air temperature)
!   print*, "variance: ", fsml_var(df%data(:,2))
!
!   ! correlation of msl and t2m
!   print*, "pearson correlation coefficent: ", fsml_pcc(df%data(:,1), df%data(:,2))

  ! exponential pdf (x=0.8, lambda=0.5)
  print*, fsml_exp_pdf(0.8_wp, lambda=0.5_wp)

  ! left-tailed p-value for normal distribution with specified mean and standard deviation
  print*, fsml_norm_cdf(2.0_wp, mu=0.3_wp, sigma=1.3_wp, tail="left")

  ! left-tailed p-value for t distribution with specified degrees of freedom
  print*, fsml_t_ppf(0.9_wp, df=20.0_wp, mu=0.2_wp, sigma=1.2_wp)

  ! genrealised pareto distribution cdf
  print*, fsml_gpd_cdf(1.9_wp, xi=1.2_wp, mu=0.6_wp, sigma=2.2_wp, tail="left")

  ! gamma distribution pdf
  print*, fsml_gamma_pdf(0.2_wp, alpha=1.2_wp, beta=0.6_wp, loc=0.0_wp)

  ! chi square distribution ppf
  print*, fsml_chi2_ppf(0.2_wp, df=10.0_wp, loc=2.0_wp, scale=1.2_wp)

  ! invalid sigma; returns sentinel and prints error message
  print*, fsml_norm_pdf(2.0_wp, mu=0.0_wp, sigma=0.0_wp)

  ! invalid tail option; returns sentinel and prints error message
  print*, fsml_norm_cdf(2.0_wp, mu=0.0_wp, sigma=1.0_wp, tail="two-sided")

  ! invalid probability option; returns sentinel and prints error message
  print*, fsml_norm_ppf(0.95_wp, mu=0.0_wp, sigma=1.0_wp)

  ! ---- pure procedure in do concurrent loop

  ! generate data
  call random_seed()
  call random_number(x1)
  call random_number(x2)

  ! pearson correlation coefficients for 10 vector pairs
  do concurrent (i=1:5)
     pcc(i) = corr(x1(:,i), x2(:,i))
  enddo
  print*, pcc

  ! ---- use of elemental function for nomal pdf
  call random_number(x)
  print*, norm_pdf_elemental(x, 0.0_wp, 1.0_wp)

end program main

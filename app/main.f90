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

  use :: fsml
  use :: fsml_typ !iso_fortran_env, wp => real64

  implicit none

  type(fsml_typ_df)  :: df
  character(len=128) :: infile
  integer            :: i
  real(wp)           :: r

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  call fsml_read_csv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

  ! print column id and names
  print*, "dataframe columns"
  do i = 1, size(df%col_id)
     print*, df%col_id(i), df%col_nm(i)
  enddo
  print*

  ! mean of first variable (msl - mean sea level pressure)
  print*, "mean: ", fsml_mean(df%data(:,1))

  ! variance of second variable (t2m - 2m air temperature)
  print*, "variance: ", fsml_var(df%data(:,2))

  ! correlation of msl and t2m
  print*, "correlation coefficent: ", fsml_corr(df%data(:,1), df%data(:,2))

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



end program main

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
  use :: iso_fortran_env, dp => real64

  implicit none

  type(fsml_typ_df)  :: df
  character(len=128) :: infile
  integer            :: i
  real(dp)           :: r

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

  ! normal pdf (x=0.8)
  print*, fsml_norm_pdf(0.8_dp)

  ! left-tailed p-value for normal distribution with specified mean and standard deviation
  print*, fsml_norm_cdf(2.0_dp, mu=0.3_dp, sigma=1.3_dp, tail="left")

  ! left-tailed p-value for t distribution with specified degrees of freedom
  print*, fsml_t_ppf(0.9_dp, df=20, mu=0.2_dp, sigma=1.2_dp)

  ! left-tailed p-value for t distribution with specified degrees of freedom
  print*, fsml_exp_ppf(0.9_dp, mu=0.5_dp, lambda=0.8_dp)


end program main

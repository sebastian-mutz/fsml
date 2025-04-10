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

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  call fsml_readcsv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

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
  print*, fsml_pdf_norm(0.8_dp)

  ! left-tailed p-value for normal distribution with specified mean and standard deviation
  print*, fsml_cdf_norm(2.0_dp, mu=0.3_dp, sigma=1.3_dp, tail="left")

end program main

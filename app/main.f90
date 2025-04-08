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

  use fsml
  implicit none

  type(fsml_typ_df)  :: df
  character(len=128) :: infile
  integer            :: i
  real               :: r

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  call fsml_readcsv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

  ! print column id and names
  print*, "dataframe columns"
  do i = 1, size(df%col_id)
     print*, df%col_id(i), df%col_nm(i)
  enddo
  print*

  ! mean
  print*, trim(df%col_nm(1)), " mean: ", fsml_mean(df%data(:,1))

  ! trend
  print*, trim(df%col_nm(2)), " variance: ", fsml_var(df%data(:,2))

  ! correlation
  print*, trim(df%col_nm(1))," - ",trim(df%col_nm(2)), " correlation: "&
       &, fsml_corr(df%data(:,1), df%data(:,2))

end program main

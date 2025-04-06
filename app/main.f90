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

  type(fsml_typ_df) :: df
  character(len=128) :: infile

  infile = "./example/data/DMC_Mutz2021_Antofagasta.csv"

  call fsml_readcsv(infile, df, labelcol=.true., labelrow=.true., delimiter=",")

end program main

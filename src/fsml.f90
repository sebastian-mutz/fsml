module fsml

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Main module.                                                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! FSML module.

! load modules
  use :: fsml_typ
  use :: fsml_sts

! basic options
  implicit none
  private

! declare public procedures
  public :: f_sts_mean, f_sts_var, f_sts_std, f_sts_cov, f_sts_trend, f_sts_corr, say_hello

contains

  subroutine say_hello
    print *, "Hello, fsml!"
  end subroutine say_hello

end module fsml

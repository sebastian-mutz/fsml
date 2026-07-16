program example_dat

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for data manipulation procedures (dat module)             |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini ! import wp; alternatively: iso_fortran_env, wp => real64

  implicit none

  real(wp) :: x(10)
  logical  :: mask(10)


  print*

  ! ---- Sampling

  ! generate array elements
  call random_number(x)

  print*, "> fixed-n sampling"

  ! sample fixed n=5
  call fsml_sample_n(size(x), 5, mask)
  print*, mask

  ! apply mask to array
  print*, pack(x, mask)


  print*, "> Poisson sampling"

  ! poisson sample with fixed probability of inclusion of element (p=0.2)
  call fsml_sample_p(size(x), 0.2_wp, mask)
  print*, mask

  ! apply mask to array
  print*, pack(x, mask)


end program example_dat

program test_dst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for distributions functions (dst module)                     |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_typ
  use stdlib_ansi, only : fg_color_cyan, fg_color_red&
                       &, style_reset, operator(//)
  implicit none

! ==== Declarations
  logical             :: status        !! status (test passed = true)
  real(wp), parameter :: tol = 1.0E-12 !! deviation tolerance; sp: 1.0E-6, dp: 1.0E-12, qp: 1.0E-30

! ==== Instructions

! ---- normal distribution

  print*, "> testing normal distribution pdf"
  status = test_norm_pdf(tol)
  call handle_status(status)

  print*, "> testing normal distribution cdf"
  status = test_norm_cdf(tol)
  call handle_status(status)

  print*, "> testing normal distribution ppf"
  status = test_norm_ppf(tol)
  call handle_status(status)

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine handle_status(status)

! ==== Description
!! Handles test status. (Message printing and stopping if needed)

! ==== Declarations
  logical, intent(in) :: status

! Instructions
  if (status) then
     print*, fg_color_cyan // "  passed" // style_reset
  else
     print*,  fg_color_red // "  > error: one or more failed" // style_reset
     stop
  endif

end subroutine handle_status


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function test_norm_pdf(tol) result(status)

! ==== Description
!! Tests normal distribution pdf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i
  integer(i4), parameter :: n = 5      !! number of tests
  real(wp)   , parameter :: x(n)  =          &
                                  [0.0_wp,   &
                                   -1.12_wp, &
                                   -2.26_wp, &
                                   2.51_wp,  &
                                   3.96_wp]
  real(wp), parameter :: mu(n)    =          &
                                  [0.0_wp,   &
                                   0.0_wp,   &
                                   0.0_wp,   &
                                   1.5_wp,   &
                                   0.0_wp]
  real(wp), parameter :: sigma(n) =          &
                                  [1.0_wp,   &
                                   1.5_wp,   &
                                   1.0_wp,   &
                                   1.0_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n)   =                                                &
                                  [0.398942280401432677940113897949965292_wp,      &
                                   0.201259693603731452039202898503809931_wp,      &
                                   3.10319322150082531601930144037777788E-0002_wp, &
                                   0.239551097728013351942510710686558937_wp,      &
                                   1.56925634065532219466981635072044412E-0004_wp]

! ==== Instructions
  status = .true.
  do i = 1, n
     res = fsml_norm_pdf(x(i), mu(i), sigma(i))
     if (abs( res - ans(i) ) .gt. tol) status = .false.
  enddo

end function test_norm_pdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function test_norm_cdf(tol) result(status)

! ==== Description
!! Tests normal distribution cdf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i
  integer(i4), parameter :: n = 5      !! number of tests
  real(wp)   , parameter :: x(n)  =          &
                                  [0.5_wp,   &
                                   -2.2_wp,  &
                                   0.6_wp,   &
                                   3.50_wp,  &
                                   6.96_wp]
  real(wp), parameter :: mu(n)    =          &
                                  [0.0_wp,   &
                                   0.0_wp,   &
                                   0.0_wp,   &
                                   1.5_wp,   &
                                   0.0_wp]
  real(wp), parameter :: sigma(n) =          &
                                  [1.0_wp,   &
                                   1.5_wp,   &
                                   1.0_wp,   &
                                   1.0_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n)   =                                                &
                                  [0.691462461274013103637704610608337744_wp,      &
                                   7.12333774139861152950218237152826767E-0002_wp, &
                                   0.725746882249926419705637214930230859_wp,      &
                                   0.977249868051820792799717362833466596_wp,      &
                                   0.999999999998298636307804321368231543_wp]

! ==== Instructions
  status = .true.
  do i = 1, n
     res = fsml_norm_cdf(x(i), mu(i), sigma(i))
     if (abs( res - ans(i) ) .gt. tol) status = .false.
  enddo

end function test_norm_cdf


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function test_norm_ppf(tol) result(status)

! ==== Description
!! Tests normal distribution ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i
  integer(i4), parameter :: n = 5      !! number of tests
  real(wp)   , parameter :: p(n)  =          &
                                  [0.01_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: mu(n)    =          &
                                  [0.0_wp,   &
                                   0.0_wp,   &
                                   0.0_wp,   &
                                   1.5_wp,   &
                                   0.0_wp]
  real(wp), parameter :: sigma(n) =          &
                                  [1.0_wp,   &
                                   1.5_wp,   &
                                   1.0_wp,   &
                                   1.0_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n)   =                                           &
                                  [-2.32634787404094822704792022705078125_wp, &
                                   -1.01173462529231983353383839130401611_wp, &
                                   0.00000000000000000000000000000000000_wp,  &
                                   2.17448975019487988902255892753601074_wp,  &
                                   2.32634787404094822704792022705078125_wp]

! ==== Instructions
  status = .true.
  do i = 1, n
     res = fsml_norm_ppf(p(i), mu(i), sigma(i))
     if (abs( res - ans(i) ) .gt. tol) status = .false.
  enddo

end function test_norm_ppf




end program test_dst

program test_nlp

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for non-linear procedures (nlp module)                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini

! ==== Declarations
  logical             :: status        !! status (test passed = true)
  real(wp), parameter :: tol = 1.0E-12 !! deviation tolerance; sp: 1.0E-6, dp: 1.0E-12, qp: 1.0E-30

! ==== Instructions

  print*
  print*, "NLP tests"

  print*, "> nlp: testing kmeans and hierarchical clustering"
  status = test_hkmeans(tol)
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
     print*, "  [✓] passed"
  else
     print*, "  [x] one or more failed"
     stop
  endif

end subroutine handle_status


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_hkmeans(tol) result(status)

! ==== Description
!! Tests for hierarchical clustering.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  integer(i4)            :: i, j
  integer(i4), parameter :: nd = 100   ! number of data points
  integer(i4), parameter :: nv = 5     ! number of variables ("features")
  integer(i4), parameter :: nc = 4     ! number of target clusters for hcluster
  real(wp)    :: x(nd, nv)
  real(wp)    :: res_gm(nv), res_cm(nv, nc)
  real(wp)    :: res_cov(nv, nv), res_sigma(nv)
  integer(i4) :: res_cl(nd), res_cc(nc)
  real(wp), parameter :: ans_gm(nv) = &
             & [-1.2717101366042016E-003_wp,  &
             &  -2.7201214498793836E-003_wp,  &
             &  -4.6363700864136747E-003_wp,  &
             &  -7.7449137582850849E-003_wp,  &
             &  -1.4947915546043295E-002_wp]
  real(wp), parameter :: ans_cm(nv, nc) = reshape([ &
             & 0.89749664912209792_wp, &
             & 1.0760158693824371_wp, &
             & 0.56645559409685020_wp, &
             & 2.5259432856793641E-002_wp, &
             & -6.4640707255686783E-002_wp, &
             & 0.26888589423987058_wp, &
             & -0.80636140696978320_wp, &
             & -0.22651104550077825_wp, &
             & -0.14296277117460343_wp, &
             & 0.16957245259062947_wp, &
             & -0.77204593059465965_wp, &
             & 1.1322621831870201_wp, &
             & -0.92753381993871753_wp, &
             & 0.43927461649247779_wp, &
             & 3.3617479606219809E-002_wp, &
             & -1.3601066289229238_wp, &
             & 4.8682716417721032E-002_wp, &
             & 1.0597843937229565_wp, &
             & -2.7531711591810303E-002_wp, &
             & -0.55169201129409007_wp &
             & ], shape=[nv,nc])
  integer(i4), parameter :: ans_cc(nc) = [20, 50, 16, 14]

! ==== Instructions

  status = .true.

  ! create data
  do i = 1, nd
    do j = 1, nv
      x(i,j) = sin(real(i*j, wp))  ! example data, replace with real input
    enddo
  enddo

  call fsml_hkmeans(x, nd, nv, nc, res_gm, res_cm, res_cl, res_cc, res_cov, res_sigma)

  do i = 1, nc
     do j = 1, nv
        if (abs( res_cm(j,i) - ans_cm(j,i) ) .gt. tol) status = .false.
     enddo
     if (abs( res_cc(i) - ans_cc(i) ) .gt. 0) status = .false.
  enddo
  do j = 1, nv
     if (abs( res_gm(j) - ans_gm(j) ) .gt. tol) status = .false.
  enddo

end function test_hkmeans



end program

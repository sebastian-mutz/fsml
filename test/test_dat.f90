program test_dat

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for statistical hypothesis test module (tst module)          |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini

! ==== Declarations
  logical             :: status        !! status (test passed = true)

! ==== Instructions

  print*
  print*, "DAT tests"

  print*, "> dat: testing fixed-n sampling"
  status = sample_n()
  call handle_status(status)

  print*, "> dat: testing Poisson sampling"
  status = sample_p()
  call handle_status(status)

  print*, "> dat: testing k-fold sampling"
  status = sample_k()
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
function sample_n() result(status)

! ==== Description
!! Tests fixed-n sampling.

! ==== Declarations
  logical                  :: status          !! status (test passed = true)
  integer(i4)              :: i, r
  integer(i4), allocatable :: seed(:)
  integer(i4), parameter   :: m = 10         !! population = 10
  integer(i4), parameter   :: n(3) = [2,0,5] !! sample sizes
  logical                  :: res_mask(m)
  logical    , parameter   :: ans_mask(m,3) = reshape([ &
                            & .false., .true. , .false., .false., .true. , & ! test 1
                            & .false., .false., .false., .false., .false., &
                            & .false., .false., .false., .false., .false., & ! test 2
                            & .false., .false., .false., .false., .false., &
                            & .true. , .false., .false., .true. , .true. , & ! test 3
                            & .false., .false., .true. , .true. , .false.  &
                            & ], shape=[m,3])

! ==== Instructions

  status = .true.

  ! set random seed for reproducibility
  call random_seed(size=r)
  allocate(seed(r))
  seed = 12345
  call random_seed(put=seed)

  ! call sampling routine for 3 tests
  do i = 1, 3
     call fsml_sample_n(m, n(i), res_mask)
     if ( .not. all(res_mask .eqv. ans_mask(:,i)) ) status = .false.
  enddo

end function sample_n


! ==================================================================== !
! -------------------------------------------------------------------- !
function sample_p() result(status)

! ==== Description
!! Tests Poisson sampling.

! ==== Declarations
  logical                  :: status          !! status (test passed = true)
  integer(i4)              :: i, r
  integer(i4), allocatable :: seed(:)
  integer(i4), parameter   :: m = 10         !! population = 10
  real(wp)   , parameter   :: p(3) = [0.1_wp,-0.2_wp,0.8_wp]
  logical                  :: res_mask(m)
  logical    , parameter   :: ans_mask(m,3) = reshape([ &
                            & .false., .false., .false., .false., .false., & ! test 1
                            & .false., .true. , .true. , .false., .false., &
                            & .false., .false., .false., .false., .false., & ! test 2
                            & .false., .false., .false., .false., .false., &
                            & .false., .false., .true. , .true. , .false., & ! test 3
                            & .true. , .false., .true. , .false., .true.   &
                            & ], shape=[m,3])

! ==== Instructions

  status = .true.

  ! set random seed for reproducibility
  call random_seed(size=r)
  allocate(seed(r))
  seed = 12345
  call random_seed(put=seed)

  ! call sampling routine for 3 tests
  do i = 1, 3
     call fsml_sample_p(m, p(i), res_mask)
     if ( .not. all(res_mask .eqv. ans_mask(:,i)) ) status = .false.
  enddo

end function sample_p


! ==================================================================== !
! -------------------------------------------------------------------- !
function sample_k() result(status)

! ==== Description
!! Tests Poisson sampling.

! ==== Declarations
  logical                  :: status          !! status (test passed = true)
  integer(i4)              :: i, r
  integer(i4), allocatable :: seed(:)
  integer(i4), parameter   :: m = 10         !! population = 10
  integer(i4), parameter   :: k = 3          !! k = 3 folds
  logical                  :: res_mask(m,k)
  logical    , parameter   :: ans_mask(m,k) = reshape([ &
                            & .true. , .true. , .false., .false., .true. , & ! mask 1
                            & .true. , .true. , .false., .true. , .true. , &
                            & .true. , .false., .true. , .true. , .true. , & ! mask 2
                            & .true. , .true. , .true. , .false., .false., &
                            & .false., .true. , .true. , .true. , .false., & ! mask 3
                            & .false., .false., .true. , .true. , .true.   &
                            & ], shape=[m,3])

! ==== Instructions

  status = .true.

  ! set random seed for reproducibility
  call random_seed(size=r)
  allocate(seed(r))
  seed = 12345
  call random_seed(put=seed)

  ! call sampling routine
  call fsml_sample_k(m, k, res_mask)
  do i = 1, k
     if ( .not. all(res_mask(:,i) .eqv. ans_mask(:,i)) ) status = .false.
  enddo

end function sample_k


end program

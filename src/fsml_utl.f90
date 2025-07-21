module fsml_utl

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Utilities module.                                                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Utilities module.

  ! load modules
  use :: fsml_ini

  ! basic options
  implicit none
  private

  ! public array operations
  public :: s_utl_rank
  ! public procedures for data type conversion
  public :: f_utl_r2c, f_utl_i2c, f_utl_c2r

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_utl_rank(x, ranks)

! ==== Description
!! Ranks all samples such that the smallest value obtains rank 1
!! and the largest rank n. Handles tied ranks and assigns average
!! rank to tied elements within one group of tied elements.

! ==== Declarations
  real(wp)                , intent(in)  :: x(:)     !! x array
  real(wp)   , allocatable, intent(out) :: ranks(:) !! ranks of x
  integer(i4), allocatable              :: idx(:)   !! index vector to sort x
  real(wp)                              :: rank_sum !! sum of tied ranks
  integer(i4)                           :: cnt      !! counter
  integer(i4)                           :: n        !! size of x
  integer(i4)                           :: i, j, k  !! loop control & flexible

! ==== Instructions

  ! allocate
  n = size(x)
  allocate(idx(n))
  allocate(ranks(n))

! ---- create index vector

  ! create index vector
  do i = 1, n
     idx(i) = i
  enddo

  ! sort index based on x
  do i = 2, n
     do j = i, 2, -1
        if (x(idx(j)) .lt. x(idx(j-1))) then
           k = idx(j)
           idx(j) = idx(j-1)
           idx(j-1) = k
        else
           exit
        endif
     enddo
  enddo

! ---- get rank sums

  ! assign ranks (with tie averaging)
  i = 1
  do while (i .le. n)

     ! initialise rank sum and reset counter for tie group
     rank_sum = real(i, kind=wp)
     cnt = 1

     ! check for ties
     do j = i + 1, n
        if (x(idx(j)) .eq. x(idx(i))) then
           rank_sum = rank_sum + real(j, kind=wp)
           cnt = cnt + 1
        else
           exit
        endif
     enddo

     ! average rank for tie group
     rank_sum = rank_sum / real(cnt, kind=wp)

     ! assign average rank to all tied elements
     do k = i, i + cnt - 1
        ranks(idx(k)) = rank_sum
     enddo

     ! advance to next group
     i = i + cnt
  enddo

  ! deallocate
  deallocate(idx)

end subroutine s_utl_rank


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_utl_r2c(r) result(c)

! ==== Description
!! Convert real to char.

! ==== Declarations
  real(wp), intent(in) :: r
  character(len=256)   :: c

! ==== Instructions
  write(c, '(F7.2)') r
  c = adjustl(c)

end function f_utl_r2c


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_utl_i2c(i) result(c)

! ==== Description
!! Convert integer to char.

! ==== Declarations
  integer(i4), intent(in) :: i
  character(len=256)      :: c

! ==== Instructions
  write(c, '(I3)') i
  c = adjustl(c)

end function f_utl_i2c


! ==================================================================== !
! -------------------------------------------------------------------- !
function f_utl_c2r(c) result(r)

! ==== Description
!! Converts char to real.

! ==== Declarations
  character(len=*), intent(in) :: c
  real(wp)                     :: r

! ==== Instructions
  read(c, *) r

end function f_utl_c2r



end module fsml_utl

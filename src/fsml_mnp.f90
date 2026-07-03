module fsml_mnp

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for data manipulation (e.g., sub-sampling, filters).        |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for basic sample statistics.

  ! load modules
  use :: fsml_ini, only: wp, i4, ieee_quiet_nan, ieee_is_nan, ieee_value
  use :: fsml_utl, only: f_utl_assign_nan

  ! basic options
  implicit none
  private

  ! declare public procedures
  ! public array operations
  public :: s_mnp_rank, s_mnp_sort

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_mnp_rank(x, ranks)

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

end subroutine s_mnp_rank


! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_mnp_sort(a_in, n, mode, idx_in, a_out, idx_out)

! ==== Description
!! Sort real array in ascending (mode=1) or descending (mode=2) order.
!! Preserves the input array. Outputs sorted array and index mapping.

! ==== Declarations
  integer(i4), intent(in)    :: n          !! number of elements
  integer(i4), intent(in)    :: mode       !! 1=ascending, 2=descending
  real(wp)   , intent(in)    :: a_in(n)    !! input array (unchanged)
  integer(i4), intent(in)    :: idx_in(n)  !! initial index mapping
  real(wp)   , intent(out)   :: a_out(n)   !! sorted output array
  integer(i4), intent(out)   :: idx_out(n) !! updated index mapping
  real(wp)                   :: tmp_a      !! swap buffer for a
  integer(i4)                :: tmp_idx    !! swap buffer for idx
  integer(i4)                :: i, j

! ==== Instructions

  ! make working copies
  a_out   = a_in
  idx_out = idx_in

  select case (mode)
  ! ascending
  case (1)
     do i = 1, n - 1
        do j = i + 1, n
           if (a_out(j) .lt. a_out(i)) then
              tmp_a      = a_out(i)
              a_out(i)   = a_out(j)
              a_out(j)   = tmp_a
              tmp_idx    = idx_out(i)
              idx_out(i) = idx_out(j)
              idx_out(j) = tmp_idx
           endif
        enddo
     enddo
  ! descending
  case (2)
     do i = 1, n - 1
        do j = i + 1, n
           if (a_out(j) .gt. a_out(i)) then
              tmp_a      = a_out(i)
              a_out(i)   = a_out(j)
              a_out(j)   = tmp_a
              tmp_idx    = idx_out(i)
              idx_out(i) = idx_out(j)
              idx_out(j) = tmp_idx
           endif
        enddo
     enddo
  ! invalid option returns sentinel
  case default
     a_out = f_utl_assign_nan()
  end select

end subroutine s_mnp_sort



end module fsml_mnp

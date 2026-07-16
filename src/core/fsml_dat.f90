module fsml_dat

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
  use :: fsml_err, only: s_err_print, fsml_error
  use :: fsml_utl, only: f_utl_assign_nan

  ! basic options
  implicit none
  private

  ! declare public procedures
  ! public array operations
  public :: s_dat_rank, s_dat_sort
  public :: s_dat_sample_n, s_dat_sample_p, s_dat_sample_k

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_dat_rank(x, ranks)

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

end subroutine s_dat_rank


! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_dat_sort(a_in, n, mode, idx_in, a_out, idx_out)

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

end subroutine s_dat_sort


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_dat_sample_n(m, n, mask)

! ==== Description
!! Subroutine for sampling a rank 1 array. It shuffles indeces
!! using the forward Fisher-Yates algorithm (as needed given n),
!! then generates an index mask from it.
!! The mask can simply be applied using the pack intrinsic function:
!! new_array = pack (old_array, mask)

! ==== Declarations
  integer(i4), intent(in)  :: m       !! size of population (array)
  integer(i4), intent(in)  :: n       !! size of sample
  logical    , intent(out) :: mask(m) !! index mask for sampled data
  integer(i4)              :: idx(m)  !! index mask for sampled data
  integer(i4)              :: i, j, k
  real(wp)                 :: r

! ==== Instructions

! ---- handle input and initialise

  ! check if size is valid
  if (n .lt. 0 .or. n .gt. m) then
     ! write error message and assign false if invalid
     call s_err_print(fsml_error(4))
     mask = .false.
     return
  endif

  ! build indeces
  do i = 1, m
     idx(i) = i
  enddo

! ---- create subsample mask

  ! partial Fisher-Yates shuffle
  do i = 1, n
     call random_number(r)
     j = i + int(r * (m - i + 1))
     k = idx(i)
     idx(i) = idx(j)
     idx(j) = k
  enddo

  ! create mask
  mask = .false.
  do i = 1, n
     mask(idx(i)) = .true.
  enddo

end subroutine s_dat_sample_n


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_dat_sample_p(m, p, mask)

! ==== Description
!! Subroutine for sampling a rank 1 array using Poisson sampling
!! (subjecting individual elements independently to Bernoulli experiments),
!! then generates an index mask for sampling.
!! The mask can simply be applied using the pack intrinsic function:
!! new_array = pack (old_array, mask)

! ==== Declarations
  integer(i4), intent(in)  :: m       !! size of population (array)
  real(wp)   , intent(in)  :: p       !! inclusion probability
  logical    , intent(out) :: mask(m) !! index mask for sampled data
  integer(i4)              :: i
  real(wp)                 :: r

! ==== Instructions

! ---- handle input and initialise

  ! check if population size is valid
  if (m .le. 0) then
     ! write error message and assign false if invalid
     call s_err_print(fsml_error(4))
     mask = .false.
     return
  endif

  ! check if value is valid
  if (p .lt. 0.0_wp .or. p .gt. 1.0_wp) then
     ! write error message and assign false if invalid
     call s_err_print(fsml_error(1))
     mask = .false.
     return
  endif

! ---- create subsample mask

  ! Poisson subsampling
  do i = 1, m
     call random_number(r)
     mask(i) = (r .lt. p)
  enddo

end subroutine s_dat_sample_p


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_dat_sample_k(m, k, mask)

! ==== Description
!! Subroutine for creating k ~equal-sized samples of a rank-1 array.
!! The array indices are shuffled using the Fisher–Yates algorithm.
!! Then, k logical masks are constructed. In each mask, the indices
!! belonging to one of the k folds (not part of the sample) are set
!! to .false. and the remaining indices are set to .true.,
!! making the masks directly suitable for k-fold cross-validation.
!!
!! In cases where m cannot be divided into exactly equal sized k folds,
!! the remainder is added to the last fold. Therefore, the last sample
!! may be smaller than the others.
!!
!! The masks can be applied using the pack intrinsic function, e.g.
!! new_array = pack(old_array, mask(:, i)) for the ith of the k samples.

! ==== Declarations
  integer(i4), intent(in)  :: m         !! size of population (array)
  integer(i4), intent(in)  :: k         !! number of subsample sets
  logical    , intent(out) :: mask(m,k) !! index mask for sampled data
  integer(i4)              :: idx(m)    !! index mask for sampled data
  integer(i4)              :: i, j, n
  real(wp)                 :: r

! ==== Instructions

! ---- handle input and initialise

  ! check if size is valid
  if (k .le. 0 .or. k .gt. m) then
     ! write error message and assign false if invalid
     call s_err_print(fsml_error(4))
     mask = .false.
     return
  endif

  ! build indeces
  do i = 1, m
     idx(i) = i
  enddo

! ---- create subsample mask

  ! full Fisher-Yates shuffle
  do i = m, 2, -1
     call random_number(r)
     j = int(r * i) + 1
     n = idx(i)
     idx(i) = idx(j)
     idx(j) = n
  enddo

  ! create masks
  n = m/k
  mask = .true.
  do j = 1, k
     if (j .eq. k) then
        do i = 1 + (j-1)*n, m
           mask(idx(i),j) = .false.
        enddo
     else
        do i = 1 + (j-1)*n, j*n
           mask(idx(i),j) = .false.
        enddo
     endif
  enddo

end subroutine s_dat_sample_k


end module fsml_dat

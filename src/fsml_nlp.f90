module fsml_nlp

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for non-linear procedures.                                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for nonlinbear procedures.

  ! load fsml modules
  use :: fsml_ini
  use :: fsml_utl
  use :: fsml_sts
  use :: fsml_err
  use :: fsml_con

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_nlp_cluster_h, s_nlp_cluster_kmeans

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_nlp_cluster_h(x, nd, nv, nc, gm, cm, corr, sigma)

! ==== Description
!! Perform agglomerative hierarchical clustering using centroid linkage
!! and the Mahalanobis distance.
!! The resulting cluster centroids can be passed to a separate k-means
!! procedure for refinement.
!! TODO: substitute distance calculation with procedure for mahalanobis distance
!! TODO: comput correlation on standardised data, if needed at all. CHECK
!! TODO: solve inefficient cluster membership storage
!! TODO: write standardiser

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)             !! input data matrix (samples × variables)
  integer(i4), intent(in)  :: nd                    !! number of data points
  integer(i4), intent(in)  :: nv                    !! number of variables
  integer(i4), intent(in)  :: nc                    !! number of clusters (target)
  real(wp)   , intent(out) :: gm(nv)                !! global means for each variable
  real(wp)   , intent(out) :: cm(nv, nc)            !! cluster centroids
  real(wp)   , intent(out) :: corr(nv, nv)          !! correlation matrix
  real(wp)   , intent(out) :: sigma(nv)             !! standard deviation per variable
  real(wp)                 :: x1(nd, nv)            !! working matrix (standardised data, mutable)
  real(wp)                 :: x2(nd, nv)            !! working matrix (standardised data, preserved)
  real(wp)                 :: vec(nd)               !! temp vectors
  real(wp)                 :: tmp_r1, tmp_r2
  real(wp)                 :: mini, total           !! minima and totals
  integer(i4)              :: cluster(nd, nd)       !! membership table
  integer(i4)              :: nde(nd)               !! number of members per cluster
  integer(i4)              :: new(nd)               !! updated member counts
  logical                  :: remain(nd)            !! flag for active clusters
  integer(i4)              :: cl1, cl2, gone
  integer(i4)              :: idx, i, j, k, l, s
  integer(i4)              :: re_row(nc)
  real(wp)                 :: pre_vector(nc), post_vector(nc)
  integer(i4)              :: cidx(nc)

! ==== Instructions

  ! ----standardise variables and compute correlation matrix
  x1 = x
  x2 = x

  ! standardise
  do i = 1, nv
     vec = x1(:, i)
     ! mean
     tmp_r1 = f_sts_mean_core(vec)
     gm(i) = tmp_r1
     ! standard deviation
     tmp_r2 = sqrt(f_sts_var_core(vec, ddf=1.0_wp))
     sigma(i) = tmp_r2
     ! standardise
     x1(:, i) = (x1(:, i) - tmp_r1) / tmp_r2
     x2(:, i) = (x2(:, i) - tmp_r1) / tmp_r2
  enddo

  ! correlation matrix
  do i = 1, nv
     do j = 1, nv
        corr(i, j) = f_sts_pcc_core(x(:, i), x(:, j))
     enddo
  enddo

  ! ---- initialise cluster membership
  do j = 1, nd
     remain(j)     = .true.
     cluster(j, 1) = j
     nde(j)        = 1
     cluster(j, 2:nd) = 0
  enddo

  ! ---- agglomerative merging
  do s = 1, nd - nc
     mini = huge(1.0_wp)
     do j = 2, nd
        do k = 1, j - 1
           if (remain(j) .and. remain(k)) then
              total = 0.0_wp
              do i = 1, nv
                 do l = 1, nv
                    total = total + (x1(j, i) - x1(k, i)) * corr(i, l) * (x1(j, l) - x1(k, l))
                 enddo
              enddo
              if (sqrt(total) .le. mini) then
                 mini = sqrt(total)
                 cl1  = j
                 cl2  = k
              endif
           endif
        enddo
     enddo

     idx = min(cl1, cl2)
     gone = max(cl1, cl2)
     remain(gone) = .false.

     ! update centroid
     do i = 1, nv
        total = 0.0_wp
        do j = 1, nde(idx)
           total = total + x2(cluster(idx, j), i)
        enddo
        do j = 1, nde(gone)
           total = total + x2(cluster(gone, j), i)
        enddo
        x1(idx, i) = total / real(nde(idx) + nde(gone), wp)
     enddo

     ! merge membership lists
     new(idx) = nde(idx) + nde(gone)
     cluster(idx, nde(idx) + 1:new(idx)) = cluster(gone, 1:nde(gone))
     nde(idx) = new(idx)
  enddo

  ! ---- sort clusters by first variable
  k = 0
  do j = 1, nd
     if (remain(j)) then
        k = k + 1
        pre_vector(k) = x1(j,1)
        re_row(k) = j
     endif
  enddo

  ! sort in descending order (2 = descending)
  call s_utl_sort(pre_vector, nc, 2, re_row, post_vector, cidx)

  do i = 1, nc
     cm(:,i) = x1(re_row(i),:)
  enddo

end subroutine s_nlp_cluster_h




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_nlp_cluster_kmeans(x, nd, nv, nc, cm_in, cm_out, cl, cc)

! ==== Description
!! Perform k-means clustering using Mahalanobis distance.
!! Starts from initial centroids (e.g. from hierarchical clustering) and iteratively
!! reassigns samples until convergence or maximum iterations reached.
!! TODO: replace distance calculation with external procedure; see above
!! TODO: think about further changes needed for successive refinement if hcluster not used

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)             !! standardised data (samples × variables)
  integer(i4), intent(in)  :: nd                    !! number of data points
  integer(i4), intent(in)  :: nv                    !! number of variables
  integer(i4), intent(in)  :: nc                    !! number of clusters
  real(wp)   , intent(in)  :: cm_in(nv, nc)         !! initial cluster centroids (variables × clusters)
  real(wp)   , intent(out) :: cm_out(nv, nc)        !! refined cluster centroids (variables × clusters)
  integer(i4), intent(out) :: cl(nd)                !! cluster assignments per sample
  integer(i4), intent(out) :: cc(nc)                !! cluster sizes
  real(wp)   , parameter   :: tol = c_conv_tol      !! convergence tolerance
  integer(i4), parameter   :: i_max = c_kmeans_i    !! max iterations
  real(wp)                 :: total, mini
  integer(i4)              :: i, j, k, l, iter, cl1
  real(wp)                 :: tmp_cm(nv, nc)        !! temp centroids
  integer(i4)              :: counter(nc)           !! count samples per cluster
  real(wp)                 :: pre_vector(nc), post_vector(nc)
  integer(i4)              :: idx_in(nc), idx_out(nc)
  real(wp)                 :: corr(nv, nv)          !! correlation matrix

! ==== Instructions

  ! ---- compute correlation matrix
  ! NOTE: make corr optional input argument. In case kmeans is called many times
  ! it prevents recalculation. Preferred way if original centroids are randomised.
  do i = 1, nv
     do j = 1, nv
        corr(i, j) = f_sts_pcc_core(x(:, i), x(:, j))
     enddo
  enddo

  ! initialise centroids output
  cm_out = cm_in

  do iter = 1, i_max
     ! assign each point to nearest centroid using Mahalanobis distance
     do j = 1, nd
        mini = huge(1.0_wp)
        cl1  = 1
        do i = 1, nc
           total = 0.0_wp
           do k = 1, nv
              do l = 1, nv
                 total = total + (x(j, l) - cm_out(l, i)) * corr(k, l) * &
                     & (x(j, k) - cm_out(k, i))
              enddo
           enddo
           if (sqrt(total) .lt. mini) then
              mini = sqrt(total)
              cl1 = i
           endif
        enddo
        cl(j) = cl1
     enddo

     ! reset tmp centroids and counts
     tmp_cm  = 0.0_wp
     counter = 0

     ! sum points in each cluster
     do j = 1, nd
        i = cl(j)
        counter(i) = counter(i) + 1
        do l = 1, nv
           tmp_cm(l, i) = tmp_cm(l, i) + x(j, l)
        enddo
     enddo

     ! average to get new centroids, avoid division by zero
     do i = 1, nc
        if (counter(i) .gt. 0) then
           do l = 1, nv
              tmp_cm(l, i) = tmp_cm(l, i) / real(counter(i), wp)
           enddo
        else
           tmp_cm(:, i) = cm_out(:, i)  ! Keep old centroid if cluster empty
        endif
     enddo

     ! check convergence
     if (all(abs(tmp_cm - cm_out) .lt. tol)) exit

     cm_out = tmp_cm
  enddo

  ! prepare for sorting clusters by first variable (descending)
  do i = 1, nc
     pre_vector(i) = cm_out(1, i)
     idx_in(i) = i
  enddo

  ! sort in decensing order (2 = descending)
  call s_utl_sort(pre_vector, nc, 2, idx_in, post_vector, idx_out)

  ! reorder centroids and counts according to sorted order
  tmp_cm = cm_out
  do i = 1, nc
     cm_out(:, i) = tmp_cm(:, idx_out(i))
     cc(i) = counter(idx_out(i))
  enddo

  ! remap cluster assignments to new sorted cluster indices
  do j = 1, nd
     do i = 1, nc
        if (cl(j) .eq. idx_out(i)) then
           cl(j) = i
           exit
        endif
     enddo
  enddo

end subroutine s_nlp_cluster_kmeans





! ! ==================================================================== !
! ! -------------------------------------------------------------------- !
! pure function f_sts_mahalanobis_core(x, y, cov) result(dist)
!
! ! ==== Description
! !! Computes the Mahalanobis distance between vectors x and y.
!
! ! ==== Declarations
!   real(wp), intent(in)  :: x(:), y(:)    ! data vectors
!   real(wp), intent(in)  :: cov(:,:)      ! covariance matrix
!   real(wp)              :: dist          ! Mahalanobis distance
!   real(wp), allocatable :: diff(:), m(:) ! temporary vectors
!   integer(i4)           :: n             ! common vector length
!
! ! ==== Instructions
!
!   ! get sizes and allocate
!   n = size(x)
!   allocate(diff(n), temp(n))
!
!   ! calculate distance
!   diff = x - y
!   temp = matmul(cov, diff)
!   dist = sqrt(dot_product(diff, temp))
!
!   ! deallocate
!   deallocate(diff, temp)
!
! end function f_sts_mahalanobis_core


end module fsml_nlp

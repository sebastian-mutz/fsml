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
  use :: fsml_err
  use :: fsml_con
  use :: fsml_sts
  use :: fsml_lin

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_nlp_hclust, s_nlp_hclust_core
  public :: s_nlp_kmeans, s_nlp_kmeans_core
  public :: s_nlp_hkmeans, s_nlp_hkmeans_core

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_nlp_hclust(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

! ==== Description
!! Impure wrapper procedure for `s_nlp_hclust_core`.

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)  !! input data matrix (samples, variables)
  integer(i4), intent(in)  :: nd         !! number of data points
  integer(i4), intent(in)  :: nv         !! number of variables
  integer(i4), intent(in)  :: nc         !! number of clusters (target)
  real(wp)   , intent(out) :: gm(nv)     !! global means for each variable
  real(wp)   , intent(out) :: cm(nv,nc)  !! cluster centroids
  integer(i4), intent(out) :: cl(nd)     !! cluster assignments for each data point
  integer(i4), intent(out) :: cc(nc)     !! cluster sizes
  real(wp)   , intent(out) :: cov(nv,nv) !! covariance matrix
  real(wp)   , intent(out) :: sigma(nv)  !! standard deviation per variable

! ==== Instructions

! ---- handle input

  ! check if argument values are valid - data points
  if (nd .le. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! issue warning for small datasets
  if (nd .le. 15) then
     ! write warning message
     call s_err_warn("[fsml warning] hcluster: small datasets can create&
                    & problems with Cholesky fractionisation.")
  endif

  ! check if argument values are valid - variable/feature number
  if (nv .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number
  if (nc .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number must be smaller than data points
  if (nc .gt. nd) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print("[fsml error] hcluster: cluster number must be&
                    & equal or less than number of data points.")
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

! ---- compute clusters

  ! call pure procedure
  call s_nlp_hclust_core(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

end subroutine s_nlp_hclust




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_nlp_hclust_core(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

! ==== Description
!! Perform agglomerative hierarchical clustering using centroid linkage
!! and the Mahalanobis distance.
!! NOTE: The procedure is exact, but slow for large nd. For most pracitcal
!! purposes, using Lance–Williams algorithm and other distances is advised.
!! TODO: Implement distance switch and L-W algorithm + approx. distance updates.

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)       !! input data matrix (samples, variables)
  integer(i4), intent(in)  :: nd              !! number of data points
  integer(i4), intent(in)  :: nv              !! number of variables
  integer(i4), intent(in)  :: nc              !! number of clusters (target)
  real(wp)   , intent(out) :: gm(nv)          !! global means for each variable
  real(wp)   , intent(out) :: cm(nv,nc)       !! cluster centroids
  integer(i4), intent(out) :: cl(nd)          !! cluster assignments for each data point
  integer(i4), intent(out) :: cc(nc)          !! cluster sizes
  real(wp)   , intent(out) :: cov(nv,nv)      !! covariance matrix
  real(wp)   , intent(out) :: sigma(nv)       !! standard deviation per variable
  real(wp)                 :: x1(nd,nv)       !! working matrix (standardised data, mutable)
  real(wp)                 :: x2(nd,nv)       !! working matrix (standardised data, preserved)
  real(wp)                 :: vec(nd)
  real(wp)                 :: tmp_r1, tmp_r2
  real(wp)                 :: mini, total
  integer(i4)              :: cluster(nd, nd), nde(nd), new(nd)
  logical                  :: remain(nd)
  integer(i4)              :: cl1, cl2, gone
  integer(i4)              :: idx, i, j, k, s
  integer(i4)              :: re_row(nc)
  real(wp)                 :: pre_vec(nc), post_vec(nc)
  integer(i4)              :: cidx(nc)
  real(wp)                 :: cluster_sum(nv, nd) !! running sum for centroids

! ==== Instructions

  ! ---- standardise variables and compute covariance matrix
  x1 = x
  x2 = x
  do i = 1, nv
     vec = x1(:, i)
     tmp_r1 = f_sts_mean_core(vec)
     gm(i)  = tmp_r1
     tmp_r2 = sqrt(f_sts_var_core(vec, ddf=1.0_wp))
     sigma(i) = tmp_r2
     x1(:, i) = (x1(:, i) - tmp_r1) / tmp_r2
     x2(:, i) = (x2(:, i) - tmp_r1) / tmp_r2
  enddo

  do i = 1, nv
     do j = 1, nv
        cov(i,j) = f_sts_cov_core(x1(:,i), x1(:,j), ddf=1.0_wp)
     enddo
  enddo

  ! ---- initialise cluster membership and running sums
  do j = 1, nd
     remain(j) = .true.
     cluster(j,1) = j
     nde(j) = 1
     cluster(j,2:nd) = 0
     cluster_sum(:, j) = x1(j, :)
  enddo

  ! ---- agglomerative merging
  do s = 1, nd - nc
     mini = huge(1.0_wp)
     do j = 2, nd
        if (.not. remain(j)) cycle       ! skip inactive clusters early
        do k = 1, j-1
           if (remain(k)) then
              ! compute Mahalanobis distance between centroids only
              total = f_lin_mahalanobis_core(cluster_sum(:,j)/real(nde(j),wp), &
                                             cluster_sum(:,k)/real(nde(k),wp), cov)
              if (total .le. mini) then
                 mini = total
                 cl1 = j
                 cl2 = k
              endif
           endif
        enddo
     enddo

     idx = min(cl1, cl2)
     gone = max(cl1, cl2)
     remain(gone) = .false.

     ! update centroid using running sum
     cluster_sum(:, idx) = cluster_sum(:, idx) + cluster_sum(:, gone)
     nde(idx) = nde(idx) + nde(gone)
     x1(idx, :) = cluster_sum(:, idx) / real(nde(idx), wp)

     ! merge membership lists
     new(idx) = nde(idx)
     cluster(idx, nde(idx)-nde(gone)+1:new(idx)) = cluster(gone, 1:nde(gone))
  enddo

  ! ---- sort clusters by first variable
  k = 0
  do j = 1, nd
     if (remain(j)) then
        k = k + 1
        pre_vec(k) = x1(j,1)
        re_row(k) = j
     endif
  enddo
  call s_utl_sort(pre_vec, nc, 2, re_row, post_vec, cidx)

  ! ---- assign centroids and memberships
  do i = 1, nc
     cm(:,i) = x1(re_row(i),:)
  enddo

  cc = 0
  do i = 1, nc
     idx = re_row(i)
     cc(i) = nde(idx)
     do j = 1, nde(idx)
        cl(cluster(idx,j)) = i
     enddo
  enddo

end subroutine s_nlp_hclust_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_nlp_kmeans(x, nd, nv, nc, cm_in, gm, cm, cl, cc, &
                             & cov, sigma, cov_in)

! ==== Description
!! Impure wrapper procedure for `s_nlp_kmeans_core`.

! ==== Declarations
  real(wp)   , intent(in)            :: x(nd, nv)      !! raw data (samples, variables)
  integer(i4), intent(in)            :: nd             !! number of data points
  integer(i4), intent(in)            :: nv             !! number of variables
  integer(i4), intent(in)            :: nc             !! number of clusters
  real(wp)   , intent(in)            :: cm_in(nv,nc)   !! initial centroids (raw, not standardised)
  real(wp)   , intent(out)           :: cm(nv, nc)     !! centroids (refined, standardised)
  real(wp)   , intent(out)           :: gm(nv)         !! global means
  integer(i4), intent(out)           :: cl(nd)         !! cluster assignments
  integer(i4), intent(out)           :: cc(nc)         !! cluster sizes
  real(wp)   , intent(out)           :: cov(nv,nv)     !! covariance matrix
  real(wp)   , intent(out)           :: sigma(nv)      !! standard deviations per variable
  real(wp)   , intent(in) , optional :: cov_in(nv,nv)  !! optional covariance matrix

! ==== Instructions

! ---- handle input

  ! check if argument values are valid - data points
  if (nd .le. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! issue warning for small datasets
  if (nd .le. 15) then
     ! write warning message
     call s_err_warn("[fsml warning] k-means: small datasets can create&
                    & problems with Cholesky fractionisation.")
  endif

  ! check if argument values are valid - variable/feature number
  if (nv .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number
  if (nc .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number must be smaller than data points
  if (nc .gt. nd) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print("[fsml error] k-means: cluster number must be&
                    & equal or less than number of data points.")
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

! ---- compute clusters

  ! call pure procedure
  call s_nlp_kmeans_core(x, nd, nv, nc, cm_in, gm, cm, cl, cc, &
                       & cov, sigma, cov_in=cov_in)

end subroutine s_nlp_kmeans




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_nlp_kmeans_core(x, nd, nv, nc, cm_in, gm, &
                                & cm, cl, cc, cov, sigma, cov_in)

! ==== Description
!! K-means clustering using Mahalanobis distance.
!! NOTE: think about only accepting standardised data to avoid redundant computation
!! in successive calls of procedure. This and repeated Cholesky fractionisation are
!! potential performance bottlenecks.

! ==== Declarations
  real(wp)   , intent(in)            :: x(nd, nv)               !! raw data (samples, variables)
  integer(i4), intent(in)            :: nd                      !! number of data points
  integer(i4), intent(in)            :: nv                      !! number of variables
  integer(i4), intent(in)            :: nc                      !! number of clusters
  real(wp)   , intent(in)            :: cm_in(nv, nc)           !! initial centroids (raw, not standardised)
  real(wp)   , intent(out)           :: cm(nv, nc)              !! centroids (refined, standardised)
  real(wp)   , intent(out)           :: gm(nv)                  !! global means
  integer(i4), intent(out)           :: cl(nd)                  !! cluster assignments
  integer(i4), intent(out)           :: cc(nc)                  !! cluster sizes
  real(wp)   , intent(out)           :: cov(nv,nv)              !! covariance matrix
  real(wp)   , intent(out)           :: sigma(nv)               !! standard deviations per variable
  real(wp)   , intent(in) , optional :: cov_in(nv,nv)           !! optional covariance matrix
  real(wp)                           :: x1(nd,nv)               !! standardised data (mutable)
  real(wp)                           :: x2(nd,nv)               !! standardised data (preserved)
  real(wp)                           :: vec(nd)                 !! temp vector for standardisation
  real(wp)                           :: tmp_cm(nv,nc)           !! working centroids
  real(wp)                           :: total, mini             !! total and minima
  real(wp)                           :: tmp_r1, tmp_r2
  real(wp)                           :: pre_vec(nc), post_vec(nc)
  integer(i4)                        :: idx_in(nc), idx_out(nc)
  integer(i4)                        :: counter(nc)             !! cluster counts
  integer(i4)                        :: i, j, k, l, cl1, iter
  integer(i4), parameter             :: i_max = c_kmeans_i      !! max iterations
  real(wp)   , parameter             :: tol = c_conv_tol        !! convergence tolerance

! ==== Instructions

  ! ---- standardise variables and compute covariance (same as cluster_h)
  x1 = x
  x2 = x

  ! standardise
  do i = 1, nv
     vec = x1(:, i)
     ! mean
     tmp_r1 = f_sts_mean_core(vec)
     gm(i)  = tmp_r1
     ! standard deviation
     tmp_r2   = sqrt(f_sts_var_core(vec, ddf=1.0_wp))
     sigma(i) = tmp_r2
     ! standardise
     x1(:, i) = (x1(:, i) - tmp_r1) / tmp_r2
     x2(:, i) = (x2(:, i) - tmp_r1) / tmp_r2
  enddo

  ! covariance matrix on standardised data
  if (present(cov_in)) then
     cov = cov_in
  else
     do i = 1, nv
        do j = 1, nv
           cov(i, j) = f_sts_cov_core(x1(:, i), x1(:, j), ddf=1.0_wp)
        enddo
     enddo
  endif

  ! ---- initialise centroids (standardise cm_in)
  do i = 1, nc
     do j = 1, nv
        cm(j, i) = (cm_in(j, i) - gm(j)) / sigma(j)
     enddo
  enddo

  ! ---- iterative refinement
  do iter = 1, i_max
     ! assign each point to nearest centroid
     do j = 1, nd
        mini = huge(1.0_wp)
        cl1  = 1
        do i = 1, nc
           total = f_lin_mahalanobis_core(x1(j, :), cm(:, i), cov)
           if (total .lt. mini) then
              mini = total
              cl1  = i
           endif
        enddo
        cl(j) = cl1
     enddo

     ! reset accumulators
     tmp_cm  = 0.0_wp
     counter = 0

     ! sum up members
     do j = 1, nd
        i = cl(j)
        counter(i) = counter(i) + 1
        do l = 1, nv
           tmp_cm(l, i) = tmp_cm(l, i) + x1(j, l)
        enddo
     enddo

     ! compute new centroids
     do i = 1, nc
        if (counter(i) .gt. 0) then
           tmp_cm(:, i) = tmp_cm(:, i) / real(counter(i), wp)
        else
           tmp_cm(:, i) = cm(:, i)
        endif
     enddo

     ! check convergence
     if (all(abs(tmp_cm - cm) .lt. tol)) exit
     cm = tmp_cm
  enddo

  ! ---- sort clusters by first variable

  ! sorting
  do i = 1, nc
     pre_vec(i) = cm(1, i)
     idx_in(i)     = i
  enddo
  call s_utl_sort(pre_vec, nc, 2, idx_in, post_vec, idx_out)

  tmp_cm = cm
  do i = 1, nc
     cm(:, i) = tmp_cm(:, idx_out(i))
     cc(i)        = counter(idx_out(i))
  enddo

  ! remap assignments to clusters
  do j = 1, nd
     do i = 1, nc
        if (cl(j) .eq. idx_out(i)) then
           cl(j) = i
           exit
        endif
     enddo
  enddo

end subroutine s_nlp_kmeans_core




! ==================================================================== !
! -------------------------------------------------------------------- !
impure subroutine s_nlp_hkmeans(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

! ==== Description
!! Impure wrapper procedure for `s_nlp_hkmeans_core`.

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)     !! input data matrix (samples, variables)
  integer(i4), intent(in)  :: nd            !! number of data points
  integer(i4), intent(in)  :: nv            !! number of variables
  integer(i4), intent(in)  :: nc            !! number of clusters (target)
  real(wp)   , intent(out) :: gm(nv)        !! global means for each variable
  real(wp)   , intent(out) :: cm(nv,nc)     !! cluster centroids
  integer(i4), intent(out) :: cl(nd)        !! cluster assignments for each data point
  integer(i4), intent(out) :: cc(nc)        !! cluster sizes
  real(wp)   , intent(out) :: cov(nv,nv)    !! covariance matrix
  real(wp)   , intent(out) :: sigma(nv)     !! standard deviation per variable

! ==== Instructions

! ---- handle input

  ! check if argument values are valid - data points
  if (nd .le. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! issue warning for small datasets
  if (nd .le. 15) then
     ! write warning message
     call s_err_warn("[fsml warning] hkmeans: small datasets can create&
                    & problems with Cholesky fractionisation.")
  endif

  ! check if argument values are valid - variable/feature number
  if (nv .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number
  if (nc .lt. 1) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print(fsml_error(1))
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

  ! check if argument values are valid - cluster number must be smaller than data points
  if (nc .gt. nd) then
     ! write error message and assign NaN/sentinel value if invalid
     call s_err_print("[fsml error] hkmeans: cluster number must be&
                    & equal or less than number of data points.")
     gm    = f_utl_assign_nan()
     cm    = f_utl_assign_nan()
     cov   = f_utl_assign_nan()
     sigma = f_utl_assign_nan()
     cl    = c_sentinel_i
     cc    = c_sentinel_i
     return
  endif

! ---- compute clusters

  ! call pure procedure
  call s_nlp_hkmeans_core(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

end subroutine s_nlp_hkmeans




! ==================================================================== !
! -------------------------------------------------------------------- !
pure subroutine s_nlp_hkmeans_core(x, nd, nv, nc, gm, cm, cl, cc, cov, sigma)

! ==== Description
!! Perform agglomerative hierarchical clustering using centroid linkage
!! and the Mahalanobis distance, then passes cluster centroids and
!! covariance matrix to kmeans cluster procedure for refinement.

! ==== Declarations
  real(wp)   , intent(in)  :: x(nd, nv)     !! input data matrix (samples, variables)
  integer(i4), intent(in)  :: nd            !! number of data points
  integer(i4), intent(in)  :: nv            !! number of variables
  integer(i4), intent(in)  :: nc            !! number of clusters (target)
  real(wp)   , intent(out) :: gm(nv)        !! global means for each variable
  real(wp)   , intent(out) :: cm(nv,nc)     !! cluster centroids
  integer(i4), intent(out) :: cl(nd)        !! cluster assignments for each data point
  integer(i4), intent(out) :: cc(nc)        !! cluster sizes
  real(wp)   , intent(out) :: cov(nv,nv)    !! covariance matrix
  real(wp)   , intent(out) :: sigma(nv)     !! standard deviation per variable
  real(wp)                 :: cm_h(nv,nc)   !! cluster centroids passed from hclust to kmeans
  real(wp)                 :: cov_h(nv,nv)  !! covariance matrix passed from hclust to kmeans

! ==== Instructions

  ! hierarchical clustering
  call s_nlp_hclust_core(x, nd, nv, nc, gm, cm_h, cl, cc, cov_h, sigma)

  ! kmeans refinement
  call s_nlp_kmeans_core(x, nd, nv, nc, cm_h, gm, cm, cl, cc, &
                       & cov, sigma, cov_in=cov_h)

end subroutine s_nlp_hkmeans_core

end module fsml_nlp

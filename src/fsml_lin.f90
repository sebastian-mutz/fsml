module fsml_lin

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Module for common statistical tests.                               |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Module for linear algebra procedures. Uses LAPACK routines (through stdlib).

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
  public :: s_lin_eof, s_lin_pca!, s_lin_lda_2class

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_eof(x, m, n, pc, eof, ew, opt, wt, eof_scaled, r2)

! ==== Description
!! Empirical Orthogonal Function (EOF) analysis

! ==== Declarations
  integer(i4), intent(in)            :: m               !! number of rows
  integer(i4), intent(in)            :: n               !! number of columns
  real(wp)   , intent(in)            :: x(m,n)          !! input data
  integer(i4), intent(in) , optional :: opt             !! 0 = covariance, 1 = correlation
  real(wp)   , intent(in) , optional :: wt(n)           !! optional weights (default = 1.0/n)
  real(wp)   , intent(out)           :: pc(m,n)         !! principal components
  real(wp)   , intent(out)           :: eof(n,n)        !! EOFs/eigenvectors (unweighted)
  real(wp)   , intent(out)           :: ew(n)           !! eigenvalues
  real(wp)   , intent(out), optional :: r2(n)           !! explained variance (fraction)
  real(wp)   , intent(out), optional :: eof_scaled(n,n) !! EOFs/eigenvectors scaled for plotting
  real(wp)                           :: eof_tmp(n,n)    !! temp storage for EOFs
  real(wp)                           :: x_w(m,n)        !! working copy of data
  integer(i4)                        :: opt_w           !! final value for opt
  real(wp)                           :: wt_w(n)         !! final values for wt
  real(wp)                           :: c(n,n)          !! covariance/correlation matrix
  real(wp)                           :: w(n)            !! matrix for eigenvalues (diagonal)
  integer(i4)                        :: nn              !! number of nonzero eigenvalues
  integer(i4)                        :: i, j            !! loop indices
  real(wp)                           :: tmp             !! temporary real

! ==== Instructions

  ! ---- handle input

  ! check matrix dimensions
  if (m .lt. 2 .or. n .lt. 1) then
     ! write error message and stop if invalid
     call s_err_print(fsml_error(3))
     error stop
  endif

  ! weight options; will default to 1/n if not specified
  tmp = real(n, kind=wp)
  wt_w = 1.0_wp / tmp
  if (present(wt)) wt_w = wt

  ! matrix options: covariance (0, default) or correlation
  opt_w = 0
  if (present(opt)) opt_w = opt

  ! make working copy of data
  x_w = x

  ! ---- construct covariance or correlation matrix

  ! prepare data
  do j = 1, n
     ! centre (get anomalies)
     tmp = f_sts_mean_core(x_w(:,j))
     x_w(:,j) = x_w(:,j) - tmp
     ! standardise if specified (if correlation matrix is to be used)
     if (opt_w .eq. 1) then
        tmp = f_sts_std_core(x_w(:,j), 0.0_wp)
        if (tmp .gt. 0.0_wp) then
           x_w(:,j) = x_w(:,j) * (1.0_wp / tmp)
        else
           call s_err_print("[fsml error] Standard deviation&
                           & is zero for a column.")
           error stop
        endif
     endif
     ! apply weights
     x_w(:,j) = x_w(:,j) * sqrt(wt_w(j))
  enddo

  ! construct matrix
  tmp = real(m - 1, kind=wp)
  c = matmul(transpose(x_w), x_w) / tmp

  ! ---- calculate outputs

  ! eigen-decomposition using stdlib eigh
  call eigh(c, w, vectors=eof_tmp)

  ! extract and reorder eigenvalues/eigenvectors in descending order
  eof = 0.0_wp
  ew  = 0.0_wp
  nn  = 0
  do i = n, 1, -1
     if (w(i) .gt. 0.0_wp) then
        nn = nn + 1
        ew(nn) = w(i)
        eof(:,nn) = eof_tmp(:,i)
     endif
  enddo

  ! compute principal components
  pc = 0.0_wp
  pc(:,1:nn) = matmul(x_w, eof(:,1:nn))

  ! undo weight scaling for EOFs
  if (nn .gt. 0) then
     do i = 1, n
        tmp = 1.0_wp / sqrt(wt_w(i))
        eof(i,1:nn) = eof(i,1:nn) * tmp
     enddo
  endif

  ! scale EOFs for plotting
  if (present(eof_scaled)) then
     if (nn .gt. 0) then
        eof_scaled = 0.0_wp
        do j = 1, nn
           tmp = sqrt(ew(j))
           eof_scaled(:,j) = eof(:,j) * tmp
        enddo
     endif
  endif

  ! explained variance (fraction)
  if (present(r2)) then
     r2 = 0.0_wp
     r2(1:nn) = ew(1:nn) / sum(ew(1:nn))
  endif

end subroutine s_lin_eof



! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_pca(x, m, n, pc, ev, ew, r2)

! ==== Description
!! Principal Component Analysis (PCA).
!! It is a special (simplified) case of EOF analysis offered as a separate
!! procedure for clarity/familiarity. It calls `s_lin_eof` with equal weights.

! ==== Declarations
  integer(i4), intent(in)            :: m       !! number of rows
  integer(i4), intent(in)            :: n       !! number of columns
  real(wp)   , intent(in)            :: x(m,n)  !! input data
  real(wp)   , intent(out)           :: pc(m,n) !! principal components
  real(wp)   , intent(out)           :: ev(n,n) !! eigenvectors (unweighted)
  real(wp)   , intent(out)           :: ew(n)   !! eigenvalues
  real(wp)   , intent(out), optional :: r2(n)   !! explained variance (fraction)
  real(wp)                           :: wt(n)   !! simple column weights

! ==== Instructions

  ! ---- handle input

  ! check matrix dimensions
  if (m .lt. 2 .or. n .lt. 1) then
     ! write error message and stop if invalid
     call s_err_print(fsml_error(3))
     error stop
  endif

  ! ---- conduct analysis

  ! set weights to 1
  wt = 1.0_wp

  ! call EOF procedure with simple weights and specify use of covariance matrix (opt=0)
  call s_lin_eof(x, m, n, pc=pc, eof=ev, ew=ew, opt=0, wt=wt, r2=r2)

end subroutine s_lin_pca


! ! ==================================================================== !
! ! -------------------------------------------------------------------- !
! subroutine s_lin_lda_2class(x, nc, nv, nd, score, sa, g, mh)
!
! ! ==== Description
! !! 2-class multivariate Linear Discriminant Analysis (LDA)
! !!
! !! Performs classification and returns:
! !! - Standardised discriminant coefficients
! !! - Reclassification accuracy
! !! - Mahalanobis distance
! !! - Discriminant criterion
!
! ! ==== Declarations
!
!   integer(i4), intent(in)           :: nc              !! number of classes (must be 2)
!   integer(i4), intent(in)           :: nv              !! number of variables
!   integer(i4), intent(in)           :: nd              !! number of datapoints per class
!   real(wp)   , intent(in)           :: x(nc,nv,nd)     !! input data (nc classes × nv variables × nd samples)
!   real(wp)   , intent(out)          :: score           !! classification accuracy
!   real(wp)   , intent(out)          :: sa(nv)          !! standardised discriminant coefficients
!   real(wp)   , intent(out)          :: g               !! discriminant criterion
!   real(wp)   , intent(out)          :: mh              !! Mahalanobis distance
!   real(wp)                          :: xmv(nc,nv)      !! group mean vectors
!   real(wp)                          :: s_g(nc,nv,nv)   !! group covariance matrices
!   real(wp)                          :: s_pool(nv,nv)   !! pooled covariance matrix
!   real(wp)                          :: s_pool_i(nv,nv) !! inverse of pooled covariance matrix
!   real(wp)                          :: a(nv)           !! discriminant vector
!   real(wp)                          :: tmpv(nv)        !! temporary vector
!   real(wp)                          :: d_pool(nc*nd)   !! pooled data for std calc
!   real(wp)                          :: tmp1, tmp2      !! temporary scalars
!   integer(i4)                       :: i, j, k, l      !! loop counters
!   ! Matrix inversion helpers
!   real(wp)                          :: m((nv*(nv+1))/2)
!   real(wp)                          :: mi((nv*(nv+1))/2)
!   real(wp)                          :: workspace(nv)
!   integer(i4)                       :: nullty, err
!
! ! ==== Instructions
!
!   ! ---- validate inputs
!   if (nc .ne. 2) then
!      call s_err_print("[fsml error] 2-class LDA: Number of classes must be 2.")
!      error stop
!   endif
!   if (nv .lt. 2) then
!      call s_err_print("[fsml error] 2-class LDA: 2+ variables required.")
!      error stop
!   endif
!
!   ! ---- compute group means and covariance matrices
!   xmv = 0.0_wp
!   s_g = 0.0_wp
!   do i = 1, nc
!      do j = 1, nv
!         tmpv(:) = x(i,j,:)  ! get variable values
!         xmv(i,j) = f_sts_mean(tmpv, nd)
!      enddo
!      do j = 1, nv
!         do k = 1, nv
!            s_g(i,j,k) = f_sts_cov(x(i,j,:), x(i,k,:), nd)
!         enddo
!      enddo
!   enddo
!
!   ! ---- compute pooled covariance matrix (unweighted, equal nd)
!   s_pool = 0.5_wp * (s_g(1,:,:) + s_g(2,:,:))
!
!   ! ---- invert pooled covariance matrix
!   k = 0
!   do i = 1, nv
!      do j = 1, i
!         k = k + 1
!         m(k) = s_pool(i,j)
!      enddo
!   enddo
!
!   call syminv(m, nv, mi, workspace, nullty, err)
!   if (err .ne. 0) then
!      call s_err_print("[fsml error] 2-class LDA: Matrix inversion failed.")
!      mh    = c_sentinel_r
!      g     = c_sentinel_r
!      sa(:) = c_sentinel_r
!      score = c_sentinel_r
!      return
!   endif
!
!   s_pool_i = 0.0_wp
!   k = 0
!   do i = 1, nv
!      do j = 1, i
!         k = k + 1
!         s_pool_i(i,j) = mi(k)
!         s_pool_i(j,i) = mi(k)
!      enddo
!   enddo
!
!   ! ---- compute discriminant vector a = S_pool⁻¹ * (μ1 - μ2)
!   a = 0.0_wp
!   do i = 1, nv
!      do j = 1, nv
!         a(i) = a(i) + s_pool_i(i,j) * (xmv(1,j) - xmv(2,j))
!      enddo
!   enddo
!
!   ! ---- standardise coefficients
!   do i = 1, nv
!      do j = 1, nd
!         d_pool((i - 1)*nd + j)         = x(1,i,j)
!         d_pool(nv*nd + (i - 1)*nd + j) = x(2,i,j)
!      enddo
!      tmp1 = f_variance(d_pool(1:2*nd), 2*nd)
!      sa(i) = a(i) * sqrt(tmp1)
!   enddo
!
!   ! ---- compute Mahalanobis distance:
!   tmp1 = 0.0_wp
!   do i = 1, nv
!      do j = 1, nv
!         tmp1 = tmp1 + (xmv(1,i) - xmv(2,i)) * s_pool_i(i,j) * (xmv(1,j) - xmv(2,j))
!      enddo
!   enddo
!   mh = sqrt(tmp1)
!
!   ! ---- compute discriminant criterion g
!   tmp1 = dot_product(a, xmv(1,:))
!   tmp2 = dot_product(a, xmv(2,:))
!   g = 0.5_wp * (tmp1 + tmp2)
!
!   ! ---- re-classification and scoring
!   score = 0.0_wp
!   do i = 1, nd
!      do j = 1, nc
!         tmp1 = dot_product(a, x(j,:,i))
!         if ((tmp1 .ge. g .and. j .eq. 1) .or. (tmp1 .lt. g .and. j .eq. 2)) then
!            score = score + 1.0_wp
!         endif
!      enddo
!   enddo
!   score = score / real(nc * nd, kind=wp)
!
! end subroutine s_lin_lda_2class




end module fsml_lin

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
  public :: s_lin_eof, s_lin_pca, s_lin_lda_2c, s_lin_ols, s_lin_ridge

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_eof(x, nd, nv, pc, eof, ew, opt, wt, eof_scaled, r2)

! ==== Description
!! Empirical Orthogonal Function (EOF) analysis

! ==== Declarations
  integer(i4), intent(in)            :: nd                !! number of rows
  integer(i4), intent(in)            :: nv                !! number of columns
  real(wp)   , intent(in)            :: x(nd,nv)          !! input data
  integer(i4), intent(in) , optional :: opt               !! 0 = covariance, 1 = correlation
  real(wp)   , intent(in) , optional :: wt(nv)            !! optional weights (default = 1.0/n)
  real(wp)   , intent(out)           :: pc(nd,nv)         !! principal components
  real(wp)   , intent(out)           :: eof(nv,nv)        !! EOFs/eigenvectors (unweighted)
  real(wp)   , intent(out)           :: ew(nv)            !! eigenvalues
  real(wp)   , intent(out), optional :: r2(nv)            !! explained variance (fraction)
  real(wp)   , intent(out), optional :: eof_scaled(nv,nv) !! EOFs/eigenvectors scaled for plotting
  real(wp)                           :: eof_tmp(nv,nv)    !! temp storage for EOFs
  real(wp)                           :: x_w(nd,nv)        !! working copy of data
  integer(i4)                        :: opt_w             !! final value for opt
  real(wp)                           :: wt_w(nv)          !! final values for wt
  real(wp)                           :: c(nv,nv)          !! covariance/correlation matrix
  real(wp)                           :: w(nv)             !! matrix for eigenvalues (diagonal)
  integer(i4)                        :: nn                !! number of nonzero eigenvalues
  integer(i4)                        :: i, j              !! loop indices
  real(wp)                           :: tmp               !! temporary real

! ==== Instructions

  ! ---- handle input

  ! check matrix dimensions
  if (nd .lt. 2 .or. nv .lt. 1) then
     ! write error message and stop if invalid
     call s_err_print(fsml_error(3))
     error stop
  endif

  ! weight options; will default to 1/nv if not specified
  tmp = real(nv, kind=wp)
  wt_w = 1.0_wp / tmp
  if (present(wt)) wt_w = wt

  ! matrix options: covariance (0, default) or correlation
  opt_w = 0
  if (present(opt)) opt_w = opt

  ! make working copy of data
  x_w = x

  ! ---- construct covariance or correlation matrix

  ! prepare data
  do j = 1, nv
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
  tmp = real(nd - 1, kind=wp)
  c = matmul(transpose(x_w), x_w) / tmp

  ! ---- calculate outputs

  ! eigen-decomposition using stdlib eigh
  call eigh(c, w, vectors=eof_tmp)

  ! extract and reorder eigenvalues/eigenvectors in descending order
  eof = 0.0_wp
  ew  = 0.0_wp
  nn  = 0
  do i = nv, 1, -1
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
     do i = 1, nv
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
subroutine s_lin_pca(x, nd, nv, pc, ev, ew, r2)

! ==== Description
!! Principal Component Analysis (PCA).
!! It is a special (simplified) case of EOF analysis offered as a separate
!! procedure for clarity/familiarity. It calls `s_lin_eof` with equal weights.

! ==== Declarations
  integer(i4), intent(in)            :: nd        !! number of rows
  integer(i4), intent(in)            :: nv        !! number of columns
  real(wp)   , intent(in)            :: x(nd,nv)  !! input data
  real(wp)   , intent(out)           :: pc(nd,nv) !! principal components
  real(wp)   , intent(out)           :: ev(nv,nv) !! eigenvectors (unweighted)
  real(wp)   , intent(out)           :: ew(nv)    !! eigenvalues
  real(wp)   , intent(out), optional :: r2(nv)    !! explained variance (fraction)
  real(wp)                           :: wt(nv)    !! simple column weights

! ==== Instructions

  ! ---- handle input

  ! check matrix dimensions
  if (nd .lt. 2 .or. nv .lt. 1) then
     ! write error message and stop if invalid
     call s_err_print(fsml_error(3))
     error stop
  endif

  ! ---- conduct analysis

  ! set weights to 1
  wt = 1.0_wp

  ! call EOF procedure with simple weights and specify use of covariance matrix (opt=0)
  call s_lin_eof(x, nd, nv, pc=pc, eof=ev, ew=ew, opt=0, wt=wt, r2=r2)

end subroutine s_lin_pca


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_lda_2c(x, nd, nv, nc, sa, g, score, mh)

! ==== Description
!! 2-class multivariate Linear Discriminant Analysis (LDA)
!!
!! Performs classification and returns:
!! - Standardised discriminant coefficients
!! - Reclassification accuracy
!! - Mahalanobis distance
!! - Discriminant criterion

! ==== Declarations
  integer(i4), intent(in)            :: nc              !! number of classes (must be 2)
  integer(i4), intent(in)            :: nv              !! number of variables
  integer(i4), intent(in)            :: nd              !! number of datapoints per class
  real(wp)   , intent(in)            :: x(nd,nv,nc)     !! input data (nd samples × nv variables × nc classes)
  real(wp)   , intent(out)           :: score           !! classification score
  real(wp)   , intent(out)           :: sa(nv)          !! standardised discriminant coefficients
  real(wp)   , intent(out)           :: g               !! discriminant criterion
  real(wp)   , intent(out), optional :: mh              !! Mahalanobis distance
  real(wp)                           :: xmv(nv,nc)      !! group mean vectors
  real(wp)                           :: s_g(nv,nv,nc)   !! group covariance matrices
  real(wp)                           :: s_pool(nv,nv)   !! pooled covariance matrix
  real(wp)                           :: s_pool_i(nv,nv) !! inverse of pooled covariance matrix
  real(wp)                           :: a(nv)           !! discriminant vector
  real(wp)                           :: d_pool(nc*nd)   !! pooled data for std calc
  real(wp)                           :: tmpv(nd)        !! temporary vector
  real(wp)                           :: tmp             !! temporary scalars
  integer(i4)                        :: i, j, k         !! loop counters
  real(wp)                           :: ew(nv)          !! eigenvalues of pooled covariance
  real(wp)                           :: ev(nv,nv)       !! eigenvectors of pooled covariance
  real(wp)                           :: ew_diag(nv,nv)  !! diagonal matrix of inverted eigenvalues

! ==== Instructions

  ! ---- validate inputs
  if (nc .ne. 2) then
     call s_err_print("[fsml error] 2-class LDA: Number of classes must be 2.")
     error stop
  endif
  if (nv .lt. 2) then
     call s_err_print("[fsml error] 2-class LDA: 2+ variables required.")
     error stop
  endif

  ! ---- compute group means and covariance matrices
  xmv = 0.0_wp
  s_g = 0.0_wp
  do i = 1, nc
     do j = 1, nv
        tmpv(:) = x(:,j,i)
        xmv(j,i) = f_sts_mean_core(tmpv)
     enddo
     do j = 1, nv
        do k = 1, nv
           ! get sample covariance (ddf set to 1)
           s_g(k,j,i) = f_sts_cov_core(x(:,j,i), x(:,k,i), ddf=1.0_wp)
        enddo
     enddo
  enddo

  ! ---- invert pooled covariance matrix using eigen-decomposition

  ! compute pooled covariance matrix (equal nd → unweighted average)
  s_pool = 0.5_wp * (s_g(:,:,1) + s_g(:,:,2))

  ! eigendecomposition of s_pool
  call eigh(s_pool, ew, vectors=ev)

  ! set returns to sentinel
  do i = 1, nv
     if (ew(i) .le. 0.0_wp) then
        if (present(mh)) mh = c_sentinel_r
        g     = c_sentinel_r
        sa(:) = c_sentinel_r
        score = c_sentinel_r
        return
     endif
  enddo

  ! construct diagonal matrix from inverted eigenvalues
  ew_diag = 0.0_wp
  do i = 1, nv
     ew_diag(i,i) = 1.0_wp / ew(i)
  enddo

  ! compute inverse of pooled covariance matrix: V * D⁻¹ * Vᵗ
  s_pool_i = matmul(ev, matmul(ew_diag, transpose(ev)))

  ! ---- discriminant coefficients

  ! compute discriminant vector a = S_pool⁻¹ * (μ1 - μ2)
  a = matmul(s_pool_i, xmv(:,1) - xmv(:,2))

  ! standardise coefficients
  do i = 1, nv
     d_pool(1:nd)       = x(1:nd,i,1)
     d_pool(nd+1:nc*nd) = x(1:nd,i,2)
     ! get sample variance (ddf set to 1)
     tmp   = f_sts_var_core(d_pool, ddf=1.0_wp)
     sa(i) = a(i) * sqrt(tmp)
  enddo

  ! ---- compute Mahalanobis distance
  if (present(mh)) then
     mh = sqrt( dot_product( xmv(:,1) - xmv(:,2), &
          & matmul( s_pool_i, xmv(:,1) - xmv(:,2) ) ) )
  endif

  ! ---- compute discriminant criterion g
  g = 0.5_wp * (dot_product(a, xmv(:,1)) + dot_product(a, xmv(:,2)))

  ! ---- re-classification and scoring
  score = 0.0_wp
  do i = 1, nd
     do j = 1, nc
        tmp = dot_product(a, x(i,:,j))
        if ((tmp .ge. g .and. j .eq. 1) .or. &
         & (tmp .lt. g .and. j .eq. 2)) then
           score = score + 1.0_wp
        endif
     enddo
  enddo
  score = score / real(nc * nd, kind=wp)

end subroutine s_lin_lda_2c




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_ols(x, y, nd, nv, b0, b, r2, y_hat, se, cov_b)

! ==== Description
!! Multiple Linear Ordinary Least Squares (OLS) Regression with intercept.
!! NOTE: OLS could be wrapper for ridge (with lambda = 0 or presence checks
!! if mande an optional argument). However, it would increase computation
!! slightly and make code less readable. OLS is often used in teaching and
!! therefore, an easily readable standalone is kept.
!!
!! Computes:
!!  - Intercept b0 (scalar)
!!  - Predictor coefficients b(nv)
!!  - Coefficient of determination R²
!!  - Standard errors se(nv)
!!  - Covariance matrix of predictors covb(nv,nv)

! ==== Declarations
  integer(i4), intent(in)            :: nd                !! number of datapoints
  integer(i4), intent(in)            :: nv                !! number of predictors (excluding intercept)
  real(wp)   , intent(in)            :: x(nd,nv)          !! predictor data matrix (no intercept column)
  real(wp)   , intent(in)            :: y(nd)             !! response vector
  real(wp)   , intent(out)           :: b0                !! intercept coefficient
  real(wp)   , intent(out)           :: b(nv)             !! predictor coefficients
  real(wp)   , intent(out)           :: r2                !! coefficient of determination R²
  real(wp)   , intent(out), optional :: y_hat(nd)         !! predicted y values
  real(wp)   , intent(out), optional :: se(nv)            !! standard errors of predictor coefficients
  real(wp)   , intent(out), optional :: cov_b(nv,nv)      !! covariance matrix of predictor coefficients
  real(wp)                           :: x1(nd,nv+1)       !! matrix with intercept column
  real(wp)                           :: xt(nv+1,nd)       !! transpose of x1
  real(wp)                           :: xtx(nv+1,nv+1)    !! XᵗX matrix
  real(wp)                           :: xtx_i(nv+1,nv+1)  !! inverse of XᵗX
  real(wp)                           :: xty(nv+1)         !! Xᵗy
  real(wp)                           :: res(nd)           !! residuals
  real(wp)                           :: sse               !! sum of squared errors
  real(wp)                           :: sst               !! total sum of squares
  real(wp)                           :: y_bar             !! mean of y
  real(wp)                           :: ew(nv+1)          !! eigenvalues
  real(wp)                           :: ev(nv+1,nv+1)     !! eigenvectors
  real(wp)                           :: ew_diag(nv+1,nv+1)!! diag matrix of 1/λ
  integer(i4)                        :: i, j              !! loop counters

! ==== Instructions

  ! ---- validate input
  if (nd .le. nv + 1) then
     call s_err_print("[fsml error] OLS: Number of observations must&
                    & exceed number of predictors + intercept.")
     error stop
  endif

  ! ---- construct matrix with intercept column (first column = 1.0)
  do i = 1, nd
     x1(i,1) = 1.0_wp
     do j = 1, nv
        x1(i,j+1) = x(i,j)
     enddo
  enddo

  ! ---- eigen-decomposition for inversion of XᵗX

  ! compute transposed matrix and XᵗX
  xt   = transpose(x1)
  xtx  = matmul(xt, x1)

  ! eigendecomposition of s_pool
  call eigh(xtx, ew, vectors=ev)

  ! check result and return sentinel if needed
  do i = 1, nv+1
     if (ew(i) .le. 0.0_wp) then
        b0        = c_sentinel_r
        b(:)      = c_sentinel_r
        r2        = c_sentinel_r
        if (present(y_hat)) y_hat = c_sentinel_r
        if (present(se))       se = c_sentinel_r
        if (present(cov_b)) cov_b = c_sentinel_r
        return
     endif
  enddo

  ! construct diagonal matrix from inverted eigenvalues
  ew_diag = 0.0_wp
  do i = 1, nv+1
     ew_diag(i,i) = 1.0_wp / ew(i)
  enddo

  ! invert XᵗX
  xtx_i = matmul(ev, matmul(ew_diag, transpose(ev)))

  ! ---- compute coefficients: full vector including intercept
  xty = matmul(xt, y)
  b0  = 0.0_wp
  b   = 0.0_wp

  ! get intercept coefficient
  b0 = dot_product(xtx_i(1,:), xty)

  ! get predictor coefficients
  do j = 1, nv
     b(j) = dot_product(xtx_i(j+1, :), xty)
  enddo

  ! ---- predicted values and residuals
  y_hat = matmul(x1, [b0, b])  ! [b0, b] concatenates b0 and b
  res   = y - y_hat

  ! ---- R² calculation

  ! sum of squared errors
  sse = sum(res**2)

  ! mean of y
  y_bar = f_sts_mean_core(y)

  ! total sum of squares
  sst = sum((y - y_bar)**2)

  ! R²
  r2 = 1.0_wp - sse / sst

  ! ---- covariance matrix of full coefficients
  if (present(cov_b)) then
     cov_b = xtx_i(2:nv+1, 2:nv+1) * (sse / real(nd - nv - 1, kind=wp))
  endif

  ! ---- standard errors for predictors (exclude intercept)
  if (present(se)) then
     do i = 1, nv
        se(i) = sqrt(cov_b(i,i))
     enddo
  endif

end subroutine s_lin_ols




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_ridge(x, y, nd, nv, lambda, b0, b, r2, y_hat, se, cov_b)

! ==== Description
!! Multiple Linear Ridge Regression (λ >= 0) with intercept.
!!
!! Computes:
!!  - Intercept b0 (scalar)
!!  - Predictor coefficients b(nv)
!!  - Coefficient of determination R²
!!  - Standard errors se(nv)   (ridge-adjusted)
!!  - Covariance matrix of predictors covb(nv,nv) (ridge-adjusted)
!!
!! Notes:
!!  - When lambda (λ) = 0, this reduces to ordinary least squares (OLS).
!!  - Ridge modifies the variance-covariance formula:
!!       cov(β) = σ² (XᵀX + λI)⁻¹ XᵀX (XᵀX + λI)⁻¹
!!    This shrinks coefficients and affects SEs.

! ==== Declarations
  integer(i4), intent(in)            :: nd                     !! number of datapoints
  integer(i4), intent(in)            :: nv                     !! number of predictors (excluding intercept)
  real(wp)   , intent(in)            :: x(nd,nv)               !! predictor data matrix (no intercept column)
  real(wp)   , intent(in)            :: y(nd)                  !! response vector
  real(wp)   , intent(in)            :: lambda                 !! ridge penalty parameter (≥ 0, non-optional)
  real(wp)   , intent(out)           :: b0                     !! intercept coefficient
  real(wp)   , intent(out)           :: b(nv)                  !! predictor coefficients
  real(wp)   , intent(out)           :: r2                     !! coefficient of determination R²
  real(wp)   , intent(out), optional :: y_hat(nd)              !! predicted y values
  real(wp)   , intent(out), optional :: se(nv)                 !! standard errors of predictor coefficients
  real(wp)   , intent(out), optional :: cov_b(nv,nv)           !! covariance matrix of predictor coefficients
  real(wp)                           :: x1(nd,nv+1)            !! matrix with intercept column
  real(wp)                           :: xt(nv+1,nd)            !! transpose of x1
  real(wp)                           :: xtx(nv+1,nv+1)         !! XᵀX matrix
  real(wp)                           :: xtx_ridge(nv+1,nv+1)   !! ridge-adjusted XᵀX matrix
  real(wp)                           :: xtx_ridge_i(nv+1,nv+1) !! inverse of ridge matrix
  real(wp)                           :: xty(nv+1)              !! Xᵀy
  real(wp)                           :: res(nd)                !! residuals
  real(wp)                           :: sse                    !! sum of squared errors
  real(wp)                           :: sst                    !! total sum of squares
  real(wp)                           :: y_bar                  !! mean of y
  real(wp)                           :: s2                     !! residual variance estimate
  real(wp)                           :: ew_r(nv+1)             !! eigenvalues (ridge)
  real(wp)                           :: ev_r(nv+1,nv+1)        !! eigenvectors (ridge)
  real(wp)                           :: ew_diag_r(nv+1,nv+1)   !! diag matrix of 1/λ (ridge)
  real(wp)                           :: cov(nv+1,nv+1)         !! covariance matrix including intercept
  integer(i4)                        :: i, j                   !! loop counters

! ==== Instructions

  ! ---- validate input
  if (nd .le. nv + 1) then
     call s_err_print("[fsml error] Ridge: Number of observations must&
                    & exceed number of predictors + intercept.")
     error stop
  endif
  if (lambda .lt. 0.0_wp) then
     call s_err_print("[fsml error] Ridge: lambda must be non-zero positive.")
     error stop
  endif

  ! ---- construct matrix with intercept column (first column = 1.0)
  do i = 1, nd
     x1(i,1) = 1.0_wp
     do j = 1, nv
        x1(i,j+1) = x(i,j)
     enddo
  enddo

  ! ---- eigen-decomposition for inversion of XᵗX

  ! compute transposed matrix and XᵀX
  xt   = transpose(x1)
  xtx  = matmul(xt, x1)

  ! build ridge-penalised XᵀX (lambda added to diagonal, except intercept term)
  ! (remains unadjusted OLS matrix if lamnbda = 0)
  xtx_ridge = xtx
  do i = 2, nv+1
     xtx_ridge(i,i) = xtx_ridge(i,i) + lambda
  enddo

  ! eigen-decomposition for inversion of ridge-adjusted XᵀX
  call eigh(xtx_ridge, ew_r, vectors=ev_r)

  ! check result and return sentinel if needed
  do i = 1, nv+1
     if (ew_r(i) .le. 0.0_wp) then
        b0        = c_sentinel_r
        b(:)      = c_sentinel_r
        r2        = c_sentinel_r
        if (present(y_hat)) y_hat = c_sentinel_r
        if (present(se))       se = c_sentinel_r
        if (present(cov_b)) cov_b = c_sentinel_r
        return
     endif
  enddo

  ! ---- construct diagonal matrix from inverted eigenvalues
  ew_diag_r = 0.0_wp
  do i = 1, nv+1
     ew_diag_r(i,i) = 1.0_wp / ew_r(i)
  enddo

  ! invert ridge-adjusted XᵀX
  xtx_ridge_i = matmul(ev_r, matmul(ew_diag_r, transpose(ev_r)))

  ! ---- compute coefficients: full vector including intercept
  xty = matmul(xt, y)
  b0  = 0.0_wp
  b   = 0.0_wp

  ! get intercept coefficient
  b0  = dot_product(xtx_ridge_i(1,:), xty)

  ! get predictor coefficients
  do j = 1, nv
     b(j) = dot_product(xtx_ridge_i(j+1,:), xty)
  enddo

  ! ---- predicted values and residuals
  if (present(y_hat)) then
     y_hat = matmul(x1, [b0, b])
     res   = y - y_hat
  else
     res   = y - matmul(x1, [b0, b])
  endif

  ! ---- R² calculation

  ! sum of squared errors
  sse   = sum(res**2)

  ! mean of y
  y_bar = f_sts_mean_core(y)

  ! total sum of squares
  sst   = sum((y - y_bar)**2)

  ! R²
  r2    = 1.0_wp - sse / sst

  ! ---- covariance matrix (ridge-adjusted)

  ! calculate only when needed; reduces to OLS covariance when lambda = 0.
  if (present(cov_b) .or. present(se)) then
     s2  = sse / real(nd - nv - 1, kind=wp)
     cov = matmul(xtx_ridge_i, matmul(xtx, xtx_ridge_i)) * s2
  endif

  ! ---- covariance matrix of full coefficients
  if (present(cov_b)) cov_b = cov(2:nv+1, 2:nv+1)

  ! ---- standard errors for predictors (exclude intercept)
  if (present(se)) then
     do i = 1, nv
        se(i) = sqrt(cov(i+1, i+1))
     enddo
  endif

end subroutine s_lin_ridge


end module fsml_lin

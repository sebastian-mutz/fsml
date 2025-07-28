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

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_lin_pca

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_lin_pca(x, m, n, opt, wt, pc, eof, ev, eof_scaled, r2)

! ==== Description
!! Empirical Orthogonal Function (EOF) analysis / Principal Component Analysis (PCA)

! ==== Declarations
  integer(i4), intent(in)            :: m               !! number of rows
  integer(i4), intent(in)            :: n               !! number of columns
  real(wp)   , intent(in)            :: x(m,n)          !! input data
  integer(i4), intent(in) , optional :: opt             !! 0 = covariance, 1 = correlation
  real(wp)   , intent(in) , optional :: wt(n)           !! optional weights (default = 1.0/n)
  real(wp)   , intent(out)           :: pc(m,n)         !! principal components
  real(wp)   , intent(out)           :: eof(n,n)        !! EOFs/eigenvectors (unweighted)
  real(wp)   , intent(out)           :: ev(n)           !! eigenvalues
  real(wp)   , intent(out), optional :: r2(n)           !! explained variance (%)
  real(wp)   , intent(out), optional :: eof_scaled(n,n) !! EOFs/eigenvectors scaled for plotting
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
     call s_err_print("[fsml error] Passed array has invalid dimensions.")
     error stop
  endif

  ! weight options
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
        tmp = f_sts_std(x_w(:,j))
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
  call eigh(c, w, vectors=eof)

  ! ---- normalise eigenvectors
  ! NOTE: not needed; eigh returns orthonormal vectors. CHECK
  !do i = 1, n
  !   tmp = sqrt(sum(eof(:,i)**2))
  !   eof(:,i) = eof(:,i) / tmp
  !enddo

  ! extract and reorder eigenvalues/eigenvectors in descending order
  eof = 0.0_wp
  ev  = 0.0_wp
  nn  = 0
  do i = n, 1, -1
     if (w(i) .le. 0.0_wp) exit
     nn = nn + 1
     ev(nn) = w(i)
     eof(:,nn) = eof(:,i)
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
           tmp = sqrt(ev(j))
           eof_scaled(:,j) = eof(:,j) * tmp
        enddo
     endif
  endif

  ! explained variance (%)
  if (present(r2)) then
     r2 = 0.0_wp
     r2(1:nn) = ev(1:nn) / sum(ev(1:nn)) * 100.0_wp
  endif

end subroutine s_lin_pca

end module fsml_lin

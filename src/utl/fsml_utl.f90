module fsml_utl

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Utilities/helper module.                                           |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! Utilities/helper module.

  ! load modules
  use :: fsml_ini, only: wp, i4, ieee_quiet_nan, ieee_is_nan, ieee_value, chol

  ! basic options
  implicit none
  private

  ! solvers
  public :: s_utl_cholesky_solve
  ! public procedures for data type conversion
  public :: f_utl_r2c, f_utl_i2c, f_utl_c2r
  ! NaN handling
  public :: f_utl_assign_nan, f_utl_is_nan

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function s_utl_cholesky_solve(a, b, n) result(x)

! ==== Description
!! Solve a * x = b for x using Cholesky factor returned by stdlib's chol().
!! a : (n,n) symmetric positive-definite
!! b : (n)

! ==== Declarations
  integer(i4), intent(in)  :: n      !! length of vectors
  real(wp)   , intent(in)  :: a(n,n) !! input square matrix, symmetric positive-definite
  real(wp)   , intent(in)  :: b(n)   !! input right-hand side vector (n)
  real(wp)                 :: x(n)   !! output solution vector (n)
  real(wp)                 :: c(n,n) !! lower-triangular Cholesky factor of a
  real(wp)                 :: s      !! temporary accumulator for partial sums in substitutions
  integer(i4)              :: i, j   !! loop indices for rows and columns

! ==== Instructions

  ! compute lower-triangular Cholesky factor c such that a = c * transpose(c)
  c = chol(a, lower = .true., other_zeroed = .true.)

  ! forward substitution: c * y = b  ->  y
  i = 1
  s = b(i)
  x(i) = s / c(i, i)
  do i = 2, n
     s = b(i)
     do j = 1, i - 1
        s = s - c(i, j) * x(j)
     enddo
  x(i) = s / c(i, i)
  enddo

  ! backward substitution: transpose(c) * x = y  ->  x
  i = n
  s = x(i)
  x(i) = s / c(i, i)
  do i = n-1, 1, -1
     s = x(i)
     do j = i + 1, n
        s = s - c(j, i) * x(j)
     enddo
     x(i) = s / c(i, i)
  enddo

end function s_utl_cholesky_solve


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


! ==================================================================== !
! -------------------------------------------------------------------- !
pure elemental function f_utl_assign_nan() result(x)

! ==== Description
!! Return quiet NaN.

! ==== Declarations
  real(wp) :: x

! ==== Instructions
  x = ieee_value(0.0_wp, ieee_quiet_nan)

end function f_utl_assign_nan


! ==================================================================== !
! -------------------------------------------------------------------- !
pure function f_utl_is_nan(x) result(res)

! ==== Description
!! Returns true if is NaN.

! ==== Declarations
  real(wp), intent(in) :: x
  logical              :: res

! ==== Instructions
  res = ieee_is_nan(x)

end function f_utl_is_nan



end module fsml_utl

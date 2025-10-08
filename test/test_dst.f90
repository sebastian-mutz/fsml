program test_dst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Tests for distributions functions (dst module)                     |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: fsml_ini
  use :: fsml_utl

! ==== Declarations
  logical             :: status        !! status (test passed = true)
  real(wp), parameter :: tol = 1.0E-12 !! deviation tolerance; sp: 1.0E-6, dp: 1.0E-12, qp: 1.0E-30

! ==== Instructions

  print*
  print*, "DST tests"

  print*, "> dst: testing all normal distribution functions (pdf, cdf, ppf)"
  status = test_norm(tol)
  call handle_status(status)

  print*, "> dst: testing all t distribution functions (pdf, cdf, ppf)"
  status = test_t(tol)
  call handle_status(status)

  print*, "> dst: testing all gamma distribution functions (pdf, cdf, ppf)"
  status = test_gamma(tol)
  call handle_status(status)

  print*, "> dst: testing all exponential distribution functions (pdf, cdf, ppf)"
  status = test_exp(tol)
  call handle_status(status)

  print*, "> dst: testing all chi-squared distribution functions (pdf, cdf, ppf)"
  status = test_chi2(tol)
  call handle_status(status)

  print*, "> dst: testing all f distribution functions (pdf, cdf, ppf)"
  status = test_f(tol)
  call handle_status(status)

  print*, "> dst: testing all generalised pareto distribution functions (pdf, cdf, ppf)"
  status = test_gpd(tol)
  call handle_status(status)

  print*, "> dst: testing results of invalid arguments (and print messages)"
  status = test_invalid()
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
     print*, "  passed"
  else
     print*, "  [error] one or more failed"
     stop
  endif

end subroutine handle_status


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_norm(tol) result(status)

! ==== Description
!! Tests normal distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp)   , parameter :: x(n)  =          &
                                  [0.5_wp,   &
                                   -2.2_wp,  &
                                   0.6_wp,   &
                                   3.50_wp,  &
                                   6.96_wp]
  real(wp)   , parameter :: p(n)  =          &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: mu(n)    =          &
                                  [0.0_wp,   &
                                   0.0_wp,   &
                                   0.0_wp,   &
                                   1.5_wp,   &
                                   0.0_wp]
  real(wp), parameter :: sigma(n) =          &
                                  [1.0_wp,   &
                                   1.5_wp,   &
                                   1.0_wp,   &
                                   1.0_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n*3) =                            &
                                  [0.35206532676429952_wp,     &
                                   9.0721654941518695E-002_wp, &
                                   0.33322460289179967_wp,     &
                                   5.3990966513188063E-002_wp, &
                                   1.2076755629578172E-011_wp, &
                                   0.69146246127401312_wp,     &
                                   7.1233377413986165E-002_wp, &
                                   0.72574688224992645_wp,     &
                                   0.97724986805182079_wp,     &
                                   0.99999999999829869_wp,     &
                                   -2.0537489106209250_wp,     &
                                   -1.0117346252923198_wp,     &
                                   0.0000000000000000_wp,      &
                                   2.1744897501948799_wp,      &
                                   2.3263478740409482_wp  ]

! ==== Instructions
  j=0
  status = .true.
  do i = 1, n
     j=j+1
     res = fsml_norm_pdf(x(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_norm_cdf(x(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo
  do i = 1, n
     j=j+1
     res = fsml_norm_ppf(p(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  enddo

end function test_norm


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_t(tol) result(status)

! ==== Description
!! Tests Student's t-distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   -2.2_wp,  &
                                   0.6_wp,   &
                                   3.50_wp,  &
                                   6.96_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: df(n) = &
                                  [1.0_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   10.0_wp,  &
                                   30.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.25464790894703254_wp,    &
                                   5.5900519948967296E-002_wp,&
                                   0.30814100972341979_wp,    &
                                   4.7836071267013227E-003_wp,&
                                   1.3402974390335678E-007_wp,&
                                   0.64758361765043326_wp,    &
                                   7.9404487903970061E-002_wp,&
                                   0.71266985630959701_wp,    &
                                   0.99713674728505741_wp,    &
                                   0.99999995057329660_wp,    &
                                   -15.894544844049960_wp,    &
                                   -0.81649658092828759_wp,   &
                                   0.0000000000000000_wp,     &
                                   0.69981206131160434_wp,    &
                                    2.4572615424403921_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_t_pdf(x(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_t_cdf(x(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_t_ppf(p(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_t


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_gamma(tol) result(status)

! ==== Description
!! Tests Gamma distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   7.5_wp,   &
                                   10.0_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: alpha(n) = &
                                  [1.0_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   7.0_wp,   &
                                   9.0_wp]
  real(wp), parameter :: beta(n) = &
                                  [1.0_wp,   &
                                   0.5_wp,   &
                                   2.0_wp,   &
                                   1.5_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.60653065971263342_wp,     &
                                   0.14652511110987343_wp,     &
                                   6.6800942890542642E-002_wp, &
                                   9.7481872093250391E-002_wp, &
                                   0.11259903214901998_wp,     &
                                   0.39346934028735420_wp,     &
                                   0.90858983894062939_wp,     &
                                   0.10882198108584371_wp,     &
                                   0.23781653702703306_wp,     &
                                   0.72306981823295779_wp,     &
                                   2.0202707316911983E-002_wp, &
                                   0.48063938155792130_wp,     &
                                   9.3418177655848922_wp,      &
                                   12.458353319134403_wp,      &
                                   17.392461335221014_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_gamma_pdf(x(i), alpha(i), beta(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_gamma_cdf(x(i), alpha(i), beta(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_gamma_ppf(p(i), alpha(i), beta(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_gamma


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_exp(tol) result(status)

! ==== Description
!! Tests Exponential distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   2.0_wp,   &
                                   3.0_wp,   &
                                   5.0_wp,   &
                                   8.0_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: lambda(n) = &
                                  [1.0_wp,   &
                                   0.5_wp,   &
                                   0.25_wp,  &
                                   1.5_wp,   &
                                   2.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.60653065971263342_wp,     &
                                   0.18393972058572117_wp,     &
                                   0.11809163818525367_wp,     &
                                   8.2962655522175049E-004_wp, &
                                   2.2507034943851823E-007_wp, &
                                   0.39346934028736658_wp,     &
                                   0.63212055882855767_wp,     &
                                   0.52763344725898531_wp,     &
                                   0.99944691562985222_wp,     &
                                   0.99999988746482527_wp,     &
                                   2.0202707317225482E-002_wp, &
                                   0.57536414490565502_wp,     &
                                   2.7725887222372596_wp,      &
                                   0.92419624074575324_wp,     &
                                   2.3025850929605385_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_exp_pdf(x(i), lambda(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_exp_cdf(x(i), lambda(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_exp_ppf(p(i), lambda(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_exp


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_chi2(tol) result(status)

! ==== Description
!! Tests Chi-squared distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   8.0_wp,   &
                                   12.0_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: df(n) = &
                                  [1.0_wp,   &
                                   2.0_wp,   &
                                   4.0_wp,   &
                                   8.0_wp,   &
                                   16.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.43939128946772238_wp,    &
                                   0.18393972058572117_wp,    &
                                   0.10260624827987350_wp,    &
                                   9.7683407406582282E-002_wp,&
                                   6.8838489020562874E-002_wp,&
                                   0.52049987781304596_wp,    &
                                   0.63212055882853913_wp,    &
                                   0.71270250481630781_wp,    &
                                   0.56652987963325563_wp,    &
                                   0.25602023954626635_wp,    &
                                   6.2845016131518605E-004_wp,&
                                   0.57536414490186871_wp,    &
                                   3.3566939800311957_wp,     &
                                   10.040622904161864_wp,     &
                                   31.982384930015542_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_chi2_pdf(x(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_chi2_cdf(x(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_chi2_ppf(p(i), df(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_chi2


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_f(tol) result(status)

! ==== Description
!! Tests F-distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   1.0_wp,   &
                                   2.0_wp,   &
                                   4.0_wp,   &
                                   6.0_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: df1(n) = &
                                  [1.0_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   10.0_wp,  &
                                   20.0_wp]
  real(wp), parameter :: df2(n) = &
                                  [2.0_wp,   &
                                   3.0_wp,   &
                                   6.0_wp,   &
                                   12.0_wp,  &
                                   24.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.35777087639996635_wp,    &
                                   0.27885480092693399_wp,    &
                                   0.16030771434147398_wp,    &
                                   1.2813306262322135E-002_wp,&
                                   4.3616061370696777E-005_wp,&
                                   0.44721359549995787_wp,    &
                                   0.53524199845510989_wp,    &
                                   0.78832567250062413_wp,    &
                                   0.98676257208008011_wp,    &
                                   0.99996975103439700_wp,    &
                                   8.0032012803599173E-004_wp,&
                                   0.31712059283108829_wp,    &
                                   0.97653639711552387_wp,    &
                                   1.4996214812128983_wp,     &
                                   2.7379972346352588_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_f_pdf(x(i), df1(i), df2(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_f_cdf(x(i), df1(i), df2(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_f_ppf(p(i), df1(i), df2(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_f


! ==================================================================== !
! -------------------------------------------------------------------- !
function test_gpd(tol) result(status)

! ==== Description
!! Tests Generalised Pareto distribution pdf, cdf, and ppf and checks against answers.
!! If answers deviate, return test (passed) status as false.

! ==== Declarations
  real(wp), intent(in)   :: tol        !! deviation tolerance
  logical                :: status     !! status (test passed = true)
  real(wp)               :: res        !! result
  integer(i4)            :: i, j
  integer(i4), parameter :: n = 5      !! number of tests per function
  real(wp), parameter :: x(n) = &
                                  [0.5_wp,   &
                                   1.0_wp,   &
                                   2.0_wp,   &
                                   5.0_wp,   &
                                   8.0_wp]
  real(wp), parameter :: p(n) = &
                                  [0.02_wp,  &
                                   0.25_wp,  &
                                   0.5_wp,   &
                                   0.75_wp,  &
                                   0.99_wp]
  real(wp), parameter :: xi(n) = &
                                  [0.5_wp,   &
                                   0.0_wp,   &
                                   -0.2_wp,  &
                                   0.1_wp,   &
                                   0.3_wp]
  real(wp), parameter :: mu(n)    =          &
                                  [0.0_wp,   &
                                   0.0_wp,   &
                                   0.0_wp,   &
                                   1.5_wp,   &
                                   0.0_wp]
  real(wp), parameter :: sigma(n) = &
                                  [1.0_wp,   &
                                   2.0_wp,   &
                                   1.5_wp,   &
                                   0.5_wp,   &
                                   1.0_wp]
  real(wp), parameter :: ans(n*3) = &
                                  [0.51200000000000001_wp,    &
                                   0.30326532985631671_wp,    &
                                   0.19280329218107004_wp,    &
                                   5.8356852566531598E-003_wp,&
                                   4.9765078425504311E-003_wp,&
                                   0.35999999999999999_wp,    &
                                   0.39346934028736658_wp,    &
                                   0.78791637860082298_wp,    &
                                   0.99503966753184481_wp,    &
                                   0.98307987333532854_wp,    &
                                   2.0305089104421636E-002_wp,&
                                   0.57536414490356180_wp,    &
                                   0.97087077527906906_wp,    &
                                   2.2434917749851753_wp,     &
                                   9.9369056851165709_wp]

! ==== Instructions
  j = 0
  status = .true.
  do i = 1, n
     j = j + 1
     res = fsml_gpd_pdf(x(i), xi(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_gpd_cdf(x(i), xi(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do
  do i = 1, n
     j = j + 1
     res = fsml_gpd_ppf(p(i), xi(i), mu(i), sigma(i))
     if (abs( res - ans(j) ) .gt. tol) status = .false.
  end do

end function test_gpd



! ==================================================================== !
! -------------------------------------------------------------------- !
function test_invalid() result(status)

! ==== Description
!! Tests on a subset of functions if invalid arguments return NaN values.

! ==== Declarations
  logical  :: status     !! status (test passed = true)
  real(wp) :: res        !! result

! ==== Instructions

  status = .true.

  ! invalid sigma for normal pdf
  res = fsml_norm_pdf(0.5_wp, 0.5_wp, -0.1_wp)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! invalid p for normal ppf
  res = fsml_norm_ppf(1.1_wp)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! invalid df for t pdf
  res = fsml_t_pdf(0.5_wp, -0.5_wp)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! invalid tail for t cdf
  res = fsml_t_cdf(0.2_wp, 10.2_wp, tail="wrong")
  if (.not. f_utl_is_nan(res)) status = .false.
  ! invalid alpha for gamma pdf
  res = fsml_gamma_pdf(0.2_wp, alpha=-0.001_wp)
  if (.not. f_utl_is_nan(res)) status = .false.
  ! invalid scale for chi2 cdf
  res = fsml_chi2_cdf(0.9_wp, 20.0_wp, scale=-2.3_wp)
  if (.not. f_utl_is_nan(res)) status = .false.

end function test_invalid


end program test_dst

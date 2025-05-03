program example_dst

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Examples for distributions functions (dst module)                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

  use :: fsml
  use :: iso_fortran_env, dp => real64

  implicit none

  print*

  ! ---- Normal Distribution

  print*, "> normal distribution"

  ! standardised normal pdf with x=1
  print*, fsml_norm_pdf(1.0_dp, mu=0.0_dp, sigma=1.0_dp)
  ! 0.24197072451914337

  ! normal distribution pdf with x=0.5, mu=0.4 and sigma=1.2
  print*, fsml_norm_pdf(0.5_dp, mu=0.4_dp, sigma=1.2_dp)
  ! 0.33129955521528498

  ! normal distribution cdf with x=2.3, mu=0.0 and signma=1.0
  print*, fsml_norm_cdf(2.3_dp, mu=0.0_dp, sigma=1.0_dp, tail="left")
  ! 0.98927588997832416

  ! normal distribution ppf with p=0.3, mu=0.0 and signma=1.0
  print*, fsml_norm_ppf(0.3_dp, mu=0.0_dp, sigma=1.0_dp)
  ! -0.52440051270878030

  print*

  ! ---- Student T Distribution

  print*, "> student t distribution"

  ! standard t pdf with x=1.5 and df=10
  print*, fsml_t_pdf(1.5_dp, df=10, mu=0.0_dp, sigma=1.0_dp)
  ! 0.12744479428709160

  ! t distribution pdf with x=0.5, df=200, ~mu=0.8 and ~sigma=1.2
  ! similar to norm due to high df
  print*, fsml_t_pdf(0.5_dp, df=200, mu=0.4_dp, sigma=1.2_dp)
  ! 0.33087996676641318

  ! get confidence interval with t distribution cdf with x=2.3, df=10, ~mu=0.0 and ~signma=1.0
  print*, fsml_t_cdf(2.3_dp, df=10, mu=0.0_dp, sigma=1.0_dp, tail="confidence")
  ! 0.95574568671571991

  ! t distribution ppf with p=0.9, df=15, ~mu=0.0 and ~signma=1.0
  print*, fsml_t_ppf(0.9_dp, df=15, mu=0.0_dp, sigma=1.0_dp)
  ! 1.3406056078565598

  print*

  ! ---- Gamma Distribution

  print*, "> gamma distribution"

!   print*, fsml_gamma_pdf(0.2_dp, alpha=1.2_dp, beta=0.6_dp, loc=0.0_dp)
!   print*, fsml_gamma_cdf(0.2_dp, alpha=1.2_dp, beta=0.6_dp, loc=0.0_dp)
!   print*, fsml_gamma_ppf(0.2_dp, alpha=1.2_dp, beta=0.6_dp, loc=0.0_dp)

  print*

  ! ---- Exponential Distribution

  print*, "> exponential distribution"

  print*

  ! ---- Generalised Pareto Distribution

  print*, "> generalised pareto distribution"

  print*

  ! ---- Chi-Squared Distribution

  print*, "> chi-squared distribution"

!   print*, fsml_chi2_pdf(0.5_dp, df=20, loc=0.5_dp, scale=1.0_dp)
!   print*, fsml_chi2_cdf(1.5_dp, df=10, loc=0.0_dp, scale=1.0_dp)
!   print*, fsml_chi2_ppf(0.2_dp, df=10, loc=2.0_dp, scale=1.2_dp)

  print*

end program example_dst

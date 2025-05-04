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

  ! standard gamma pdf with x=2.3, defaults for shape parameters (alpha=1.0, beta=1.0) and loc
  print*, fsml_gamma_pdf(2.3_dp, alpha=1.0_dp, beta=1.0_dp, loc=0.0_dp)
  ! 0.10025884372280375

  ! standard gamma pdf with x=2.3, defaults for shape parameters and location shift loc=0.5
  print*, fsml_gamma_pdf(2.3_dp, alpha=1.0_dp, beta=1.0_dp, loc=0.5_dp)
  ! 0.16529888822158656

  ! gamma cdf with x=1.2, alpha=1.2, beta=0.6 and no location shift, and right end of the tail (survival function)
  print*, fsml_gamma_cdf(1.2_dp, alpha=1.2_dp, beta=0.6_dp, loc=0.0_dp, tail="right")
  ! 0.18230123290900657

  ! gamma ppf for p=0.5, alpha=0.5, beta=1.1 and no location shift
  print*, fsml_gamma_ppf(0.5_dp, alpha=0.5_dp, beta=1.1_dp, loc=0.0_dp)
  ! 0.25021503271636902


  print*

  ! ---- Exponential Distribution

  print*, "> exponential distribution"

  ! standard exponential pdf with x=1.3 and defaults for lambda=1.0 and loc=0.0
  print*, fsml_exp_pdf(1.3_dp, lambda=1.0_dp, loc=0.0_dp)
  ! 0.27253179303401259

  ! exponential pdf with x=1.3 and defaults for lambda=0.2 and loc=1.0
  print*, fsml_exp_pdf(1.3_dp, lambda=0.2_dp, loc=1.0_dp)
  ! 0.18835290671684976

  ! exponential cdf with x=0.1 and defaults for lambda=1.0 and loc=0.0
  print*, fsml_exp_cdf(1.3_dp, lambda=1.0_dp, loc=0.0_dp)
  ! 0.72746820696598746

  ! exponential ppf with p=0.95 and defaults for lambda=1.0 and loc=0.0
  print*, fsml_exp_ppf(0.5_dp, lambda=1.0_dp, loc=0.0_dp)
  ! 0.69314718055931490

  print*

  ! ---- Chi-Squared Distribution

  print*, "> chi-squared distribution"

  ! chi squared pdf with x=20.0, df=20, loc=0.5 and scale=1.0
  print*, fsml_chi2_pdf(20.0_dp, df=20, loc=0.5_dp, scale=1.0_dp)
  ! 6.3955413942221373E-002

  ! chi squared pdf with x=5.1, df=10, loc=0.0 and scale=0.5
  print*, fsml_chi2_pdf(5.1_dp, df=10, loc=0.0_dp, scale=0.5_dp)
  ! 0.17185714984072062

  ! chi squared cdf with x=11.5, df=10, loc=0.0 and scale=1.0
  print*, fsml_chi2_cdf(11.5_dp, df=10, loc=0.0_dp, scale=1.0_dp)
  ! 0.68008856946173257

  ! chi squared ppf with p=0.2, df=10, loc=2.0 and scale=1.2
  print*, fsml_chi2_ppf(0.2_dp, df=10, loc=2.0_dp, scale=1.2_dp)
  ! 9.4148951072402269

  print*

  ! ---- Generalised Pareto Distribution

  print*, "> generalised pareto distribution"

  ! GPD pdf with x=1.1, xi=1.2, mu=0.0 and sigma=1.0
  print*, fsml_gpd_pdf(1.1_dp, xi=1.2_dp, mu=0.0_dp, sigma=1.0_dp)
  ! 0.21376603130006952

  ! GPD pdf with x=1.9, xi=0.2, mu=0.0 and sigma=1.2
  print*, fsml_gpd_pdf(1.9_dp, xi=0.2_dp, mu=0.0_dp, sigma=1.2_dp)
  ! 0.15994243683480086

  ! GPD cdf with x=2.1, xi=2.7, mu=1.2 and sigma=1.0 (standard left tail)
  print*, fsml_gpd_cdf(2.1_dp, xi=2.7_dp, mu=1.2_dp, sigma=1.0_dp, tail="left")
  ! 0.36650539816689109

  ! GPD ppf with p=0.2, xi=0.7, mu=0.0 and sigma=1.0 (standard left tail)
  print*, fsml_gpd_ppf(0.2_dp, xi=0.7_dp, mu=1.0_dp, sigma=1.0_dp)
  ! 1.2415150853975381

  print*


end program example_dst

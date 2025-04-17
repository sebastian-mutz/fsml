module fsml

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Main module; provides interfaces.                                  |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! FSML interface module.

  ! load modules
  use :: fsml_typ
  use :: fsml_sts
  use :: fsml_dst
  use :: fsml_utl

  ! basic options
  implicit none
  private

  ! declare public statistics procedures
  public :: fsml_mean, fsml_var, fsml_std, fsml_cov, fsml_reg, fsml_corr
  ! declare public distribution procedures
  public :: fsml_norm_pdf, fsml_norm_cdf, fsml_norm_ppf
  public :: fsml_t_pdf, fsml_t_cdf, fsml_t_ppf
  public :: fsml_exp_pdf, fsml_exp_cdf, fsml_exp_ppf
  public :: fsml_gpd_pdf, fsml_gpd_cdf, fsml_gpd_ppf
  ! declare public utility procedures
  public :: fsml_read_csv
  ! declare public derived types
  public :: fsml_typ_df

! ==== Interfaces (API)
! API: interfaces to offer simpler, consistent public procedure names

! ---- Basic Statistics

! arithmetic mean
interface fsml_mean
  module procedure f_sts_mean
end interface

! variance
interface fsml_var
  module procedure f_sts_var
end interface

! standard deviation
interface fsml_std
  module procedure f_sts_std
end interface

! covariance
interface fsml_cov
  module procedure f_sts_cov
end interface

! linear trend (regression coefficient)
interface fsml_reg
  module procedure f_sts_reg
end interface

! Pearson correlation coefficient
interface fsml_corr
  module procedure f_sts_corr
end interface

! ---- Statistical Distributions

! normal distribution pdf
interface fsml_norm_pdf
  module procedure f_dst_norm_pdf
end interface

! normal distribution cdf
interface fsml_norm_cdf
  module procedure f_dst_norm_cdf
end interface

! normal distribution ppf
interface fsml_norm_ppf
  module procedure f_dst_norm_ppf
end interface

! t distribution pdf
interface fsml_t_pdf
  module procedure f_dst_t_pdf
end interface

! t distribution cdf
interface fsml_t_cdf
  module procedure f_dst_t_cdf
end interface

! t distribution ppf
interface fsml_t_ppf
  module procedure f_dst_t_ppf
end interface

! exponential distribution pdf
interface fsml_exp_pdf
  module procedure f_dst_exp_pdf
end interface

! exponential distribution cdf
interface fsml_exp_cdf
  module procedure f_dst_exp_cdf
end interface

! exponential distribution ppf
interface fsml_exp_ppf
  module procedure f_dst_exp_ppf
end interface

! generalised pareto distribution pdf
interface fsml_gpd_pdf
  module procedure f_dst_gpd_pdf
end interface

! generalised pareto distribution cdf
interface fsml_gpd_cdf
  module procedure f_dst_gpd_cdf
end interface

! generalised pareto distribution ppf
interface fsml_gpd_ppf
  module procedure f_dst_gpd_ppf
end interface

! ---- Utilities

! read csv file into dataframe
interface fsml_read_csv
  module procedure s_utl_read_csv
end interface

end module fsml

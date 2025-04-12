module fsml

! |--------------------------------------------------------------------|
! | fsml - fortran statistics and machine learning library             |
! |                                                                    |
! | about                                                              |
! | -----                                                              |
! | Main module.                                                       |
! |                                                                    |
! | license : MIT                                                      |
! | author  : Sebastian G. Mutz (sebastian@sebastianmutz.com)          |
! |--------------------------------------------------------------------|

! FORD
!! FSML module.

  ! load modules
  use :: fsml_typ
  use :: fsml_sts
  use :: fsml_dst
  use :: fsml_utl

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: fsml_mean, fsml_var, fsml_std, fsml_cov, fsml_reg, fsml_corr
  public :: fsml_pdf_norm, fsml_cdf_norm, fsml_pdf_t, fsml_cdf_t
  public :: fsml_readcsv
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
interface fsml_pdf_norm
  module procedure f_dst_pdf_norm
end interface

! normal distribution cdf
interface fsml_cdf_norm
  module procedure f_dst_cdf_norm
end interface

! student t distribution pdf
interface fsml_pdf_t
  module procedure f_dst_pdf_t
end interface

! student t distribution cdf
interface fsml_cdf_t
  module procedure f_dst_cdf_t
end interface

! ---- Utilities

! read csv file into dataframe
interface fsml_readcsv
  module procedure s_utl_readcsv
end interface

end module fsml

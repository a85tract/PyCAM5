module cloud_fraction

  use ap_cldfrc_scheme, only: &
       cldfrc_readnl,          &
       cldfrc_register,        &
       cldfrc_init,            &
       cldfrc_getparams,       &
       cldfrc => cldfrc_run,   &
       cldfrc_fice

  implicit none
  private

  public :: cldfrc_readnl
  public :: cldfrc_register
  public :: cldfrc_init
  public :: cldfrc_getparams
  public :: cldfrc
  public :: cldfrc_fice

end module cloud_fraction

! Compatibility facade for the active McICA RRTMG longwave kernel.
module rrtmg_lw_rad

  use ap_rrtmg_lw_scheme, only: rrtmg_lw => rrtmg_lw_run, inatm

  implicit none
  private

  public :: rrtmg_lw
  public :: inatm

end module rrtmg_lw_rad

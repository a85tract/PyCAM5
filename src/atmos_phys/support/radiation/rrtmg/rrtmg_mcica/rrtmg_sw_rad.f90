! Compatibility facade for the active McICA RRTMG shortwave kernel.
module rrtmg_sw_rad

  use ap_rrtmg_sw_scheme, only: rrtmg_sw => rrtmg_sw_run, inatm_sw

  implicit none
  private

  public :: rrtmg_sw
  public :: inatm_sw

end module rrtmg_sw_rad

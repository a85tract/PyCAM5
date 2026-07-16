! CAM compatibility facade for the decoupled RRTMG shortwave driver.
! The process body, initialization routine, and saved solar-band state share
! one scheme-module instance and are re-exported under the historical API.
module radsw

  use ap_rad_rrtmg_sw_scheme, only: radsw_init, &
       rad_rrtmg_sw => rad_rrtmg_sw_run

  implicit none
  private

  public :: radsw_init
  public :: rad_rrtmg_sw

end module radsw

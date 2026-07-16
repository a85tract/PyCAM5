! CAM compatibility facade for the decoupled RRTMG longwave driver.
! The process body, initialization routine, and saved top-level state share
! one scheme-module instance and are re-exported under the historical API.
module radlw

  use ap_rad_rrtmg_lw_scheme, only: radlw_init, &
       rad_rrtmg_lw => rad_rrtmg_lw_run

  implicit none
  private

  public :: radlw_init
  public :: rad_rrtmg_lw

end module radlw

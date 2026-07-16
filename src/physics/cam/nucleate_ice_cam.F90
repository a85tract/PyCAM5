module nucleate_ice_cam

  use ap_nucleate_ice_cam_calc_scheme, only: &
       use_preexisting_ice,                    &
       nucleate_ice_cam_readnl,                &
       nucleate_ice_cam_register,              &
       nucleate_ice_cam_init,                  &
       nucleate_ice_cam_calc => nucleate_ice_cam_calc_run

  implicit none
  private

  public :: use_preexisting_ice
  public :: nucleate_ice_cam_readnl
  public :: nucleate_ice_cam_register
  public :: nucleate_ice_cam_init
  public :: nucleate_ice_cam_calc

end module nucleate_ice_cam

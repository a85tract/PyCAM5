module dust_sediment_mod

  use ap_dust_sediment_tend_scheme, only: dust_sediment_vel, &
       dust_sediment_tend => dust_sediment_tend_run

  implicit none
  private

  public :: dust_sediment_vel
  public :: dust_sediment_tend

end module dust_sediment_mod

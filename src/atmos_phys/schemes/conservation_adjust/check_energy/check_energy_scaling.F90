module ap_check_energy_scaling_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: check_energy_scaling_run

contains

  !> \section arg_table_check_energy_scaling_run Argument Table
  !! \htmlinclude check_energy_scaling_run.html
  !!
  subroutine check_energy_scaling_run(ncol, pver, scaling_dycore, errmsg, errflg)
    integer,  intent(in)  :: ncol, pver
    real(r8), intent(out) :: scaling_dycore(:,:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! PI-atm does not have the newer CAM cp/cv energy-formula conversion.
    ! Identity is the exact compatibility representation of the old behavior.
    scaling_dycore(:ncol,:pver) = 1._r8
  end subroutine check_energy_scaling_run

end module ap_check_energy_scaling_scheme

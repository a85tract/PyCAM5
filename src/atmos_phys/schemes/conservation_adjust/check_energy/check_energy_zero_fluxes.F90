module ap_check_energy_zero_fluxes_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: check_energy_zero_fluxes_run

contains

  !> \section arg_table_check_energy_zero_fluxes_run Argument Table
  !! \htmlinclude check_energy_zero_fluxes_run.html
  !!
  subroutine check_energy_zero_fluxes_run(nvalues, zero_flux, errmsg, errflg)
    integer,  intent(in)  :: nvalues
    real(r8), intent(out) :: zero_flux(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0
    zero_flux(:nvalues) = 0._r8
  end subroutine check_energy_zero_fluxes_run

end module ap_check_energy_zero_fluxes_scheme

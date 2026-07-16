module ap_check_energy_save_teout_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: check_energy_save_teout_run

contains

  !> \section arg_table_check_energy_save_teout_run Argument Table
  !! \htmlinclude check_energy_save_teout_run.html
  !!
  subroutine check_energy_save_teout_run(nvalues, te_cur, teout, errmsg, errflg)
    integer,  intent(in)  :: nvalues
    real(r8), intent(in)  :: te_cur(:)
    real(r8), intent(out) :: teout(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0
    teout(:nvalues) = te_cur(:nvalues)
  end subroutine check_energy_save_teout_run

end module ap_check_energy_save_teout_scheme

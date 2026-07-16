module ap_check_energy_gmean_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: check_energy_gmean_run

contains

  !> \section arg_table_check_energy_gmean_run Argument Table
  !! \htmlinclude check_energy_gmean_run.html
  !!
  subroutine check_energy_gmean_run(has_chunks, dtime, gravit, &
       teinp_local, teout_local, psurf_local, ptopb, &
       teinp_glob, teout_glob, psurf_glob, ptopb_glob, tedif_glob, heat_glob, &
       errmsg, errflg)
    logical,  intent(in)    :: has_chunks
    real(r8), intent(in)    :: dtime, gravit
    real(r8), intent(in)    :: teinp_local, teout_local, psurf_local, ptopb
    real(r8), intent(inout) :: teinp_glob, teout_glob, psurf_glob, ptopb_glob
    real(r8), intent(inout) :: tedif_glob, heat_glob
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! The CAM adapter performs the MPI/global reduction.  These scalar
    ! operations are unchanged from the PI-atm check_energy_gmean routine.
    if (has_chunks) then
       teinp_glob = teinp_local
       teout_glob = teout_local
       psurf_glob = psurf_local
       ptopb_glob = ptopb

       tedif_glob = teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)
    else
       heat_glob = 0._r8
    end if
  end subroutine check_energy_gmean_run

end module ap_check_energy_gmean_scheme

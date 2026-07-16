! CAM compatibility adapter for the decoupled modal aerosol nucleation
! process.  The public CAM entry points are unchanged.
module modal_aero_newnuc

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use ppgrid, only : pcols, pver
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_newnuc_scheme, only : modal_aero_newnuc_sub_run, &
       modal_aero_newnuc_init

  implicit none
  private

  integer, parameter :: pcnstxx = gas_pcnst

  public :: modal_aero_newnuc_sub, modal_aero_newnuc_init

contains

  subroutine modal_aero_newnuc_sub(                              &
       lchnk, ncol, nstep, loffset, deltat, t, pmid, pdel, zm,   &
       pblh, qv, cld, q, del_h2so4_gasprod, del_h2so4_aeruptk)

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: zm(pcols,pver)
    real(r8), intent(in) :: pblh(pcols)
    real(r8), intent(in) :: qv(pcols,pver)
    real(r8), intent(in) :: cld(ncol,pver)
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: del_h2so4_gasprod(ncol,pver)
    real(r8), intent(in) :: del_h2so4_aeruptk(ncol,pver)

    call t_startf('ap_modal_aero_newnuc_sub_run')
    call modal_aero_newnuc_sub_run(                              &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, loffset,      &
         deltat, t, pmid, pdel, zm,                              &
         pblh, qv, cld, q, del_h2so4_gasprod, del_h2so4_aeruptk)
    call t_stopf('ap_modal_aero_newnuc_sub_run')

  end subroutine modal_aero_newnuc_sub

end module modal_aero_newnuc

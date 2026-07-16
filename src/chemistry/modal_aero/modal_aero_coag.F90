! CAM compatibility adapter for the decoupled modal aerosol coagulation
! process.  Initialization and process tables remain owned by one scheme
! module instance and are re-exported through the historical CAM module.
module modal_aero_coag

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use ppgrid, only : pcols, pver
  use modal_aero_data, only : ntot_amode
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_coag_scheme, only : modal_aero_coag_sub_run, &
       modal_aero_coag_init, pair_option_acoag, maxpair_acoag,   &
       maxspec_acoag, npair_acoag, modefrm_acoag, modetoo_acoag, &
       modetooeff_acoag, nspecfrm_acoag, lspecfrm_acoag,         &
       lspectoo_acoag

  implicit none
  private

  integer, parameter :: pcnstxx = gas_pcnst

  public :: modal_aero_coag_sub, modal_aero_coag_init
  public :: pair_option_acoag, maxpair_acoag, maxspec_acoag
  public :: npair_acoag, modefrm_acoag, modetoo_acoag
  public :: modetooeff_acoag, nspecfrm_acoag
  public :: lspecfrm_acoag, lspectoo_acoag

contains

  subroutine modal_aero_coag_sub(                               &
       lchnk, ncol, nstep, loffset, deltat_main, t, pmid, pdel, &
       q, dgncur_a, dgncur_awet, wetdens_a)

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    real(r8), intent(in) :: deltat_main
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dgncur_a(pcols,pver,ntot_amode)
    real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)
    real(r8), intent(in) :: wetdens_a(pcols,pver,ntot_amode)

    call t_startf('ap_modal_aero_coag_sub_run')
    call modal_aero_coag_sub_run(                               &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, ntot_amode, &
         loffset, deltat_main, t, pmid,                         &
         pdel, q, dgncur_a, dgncur_awet, wetdens_a)
    call t_stopf('ap_modal_aero_coag_sub_run')

  end subroutine modal_aero_coag_sub

end module modal_aero_coag

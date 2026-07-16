! CAM compatibility adapter for the decoupled modal aerosol gas-aerosol
! exchange process.  Initialization and shared process state live with the
! scheme so all users observe the original single module state.
module modal_aero_gasaerexch

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use ppgrid, only : pcols, pver
  use modal_aero_data, only : ntot_amode
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_gasaerexch_scheme, only : &
       modal_aero_gasaerexch_sub_run, modal_aero_gasaerexch_init, &
       maxspec_pcage, npair_pcage, modefrm_pcage, modetoo_pcage, &
       nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage, &
       n_so4_monolayers_pcage, dr_so4_monolayers_pcage, &
       soa_equivso4_factor

  implicit none
  private

  integer, parameter :: pcnstxx = gas_pcnst

  public :: modal_aero_gasaerexch_sub, modal_aero_gasaerexch_init
  public :: maxspec_pcage, npair_pcage, modefrm_pcage, modetoo_pcage
  public :: nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage
  public :: n_so4_monolayers_pcage, dr_so4_monolayers_pcage
  public :: soa_equivso4_factor

contains

  subroutine modal_aero_gasaerexch_sub(                         &
       lchnk, ncol, nstep, loffset, deltat, t, pmid, pdel, qh2o, &
       troplev, q, qqcw, dqdt_other, dqqcwdt_other, dgncur_a,    &
       dgncur_awet, sulfeq)

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    integer, intent(in) :: troplev(pcols)
    real(r8), intent(in) :: deltat
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: qh2o(pcols,pver)
    real(r8), intent(in) :: dgncur_a(pcols,pver,ntot_amode)
    real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)
    real(r8), pointer :: sulfeq(:,:,:)

    call t_startf('ap_modal_aero_gasaerexch_sub_run')
    call modal_aero_gasaerexch_sub_run(                         &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, ntot_amode,  &
         loffset, deltat, t, pmid, pdel,                        &
         qh2o, troplev, q, qqcw, dqdt_other, dqqcwdt_other,     &
         dgncur_a, dgncur_awet, sulfeq)
    call t_stopf('ap_modal_aero_gasaerexch_sub_run')

  end subroutine modal_aero_gasaerexch_sub

end module modal_aero_gasaerexch

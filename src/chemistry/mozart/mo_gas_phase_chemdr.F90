! CAM compatibility adapter for the decoupled gas-phase chemistry driver.
module mo_gas_phase_chemdr

  use shr_kind_mod, only : r8 => shr_kind_r8
  use constituents, only : pcnst
  use ppgrid, only : pcols, pver, pverp
  use physics_buffer, only : physics_buffer_desc
  use perf_mod, only : t_startf, t_stopf
  use gas_phase_chemdr_scheme, only : gas_phase_chemdr_inti, &
       gas_phase_chemdr_run, map2chm

  implicit none
  private

  public :: gas_phase_chemdr
  public :: gas_phase_chemdr_inti
  public :: map2chm

contains

  subroutine gas_phase_chemdr(lchnk, ncol, imozart, q, &
       phis, zm, zi, calday, tfld, pmid, pdel, pint, cldw, troplev, &
       ncldwtr, ufld, vfld, delt, ps, xactive_prates, fsds, ts, asdir, &
       ocnfrac, icefrac, precc, precl, snowhland, ghg_chem, latmapback, &
       chem_name, drydepflx, cflx, qtend, pbuf)

    integer, intent(in) :: lchnk, ncol, imozart
    real(r8), intent(in) :: q(pcols,pver,pcnst)
    real(r8), intent(in) :: phis(pcols)
    real(r8), intent(in) :: zm(pcols,pver)
    real(r8), intent(in) :: zi(pcols,pver+1)
    real(r8), intent(in) :: calday
    real(r8), intent(in) :: tfld(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: pint(pcols,pver+1)
    real(r8), intent(in) :: cldw(pcols,pver)
    integer, intent(in) :: troplev(pcols)
    real(r8), intent(in) :: ncldwtr(pcols,pver)
    real(r8), intent(in) :: ufld(pcols,pver)
    real(r8), intent(in) :: vfld(pcols,pver)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: ps(pcols)
    logical, intent(in) :: xactive_prates
    real(r8), intent(in) :: fsds(pcols)
    real(r8), intent(in) :: ts(pcols)
    real(r8), intent(in) :: asdir(pcols)
    real(r8), intent(in) :: ocnfrac(pcols)
    real(r8), intent(in) :: icefrac(pcols)
    real(r8), intent(in) :: precc(pcols)
    real(r8), intent(in) :: precl(pcols)
    real(r8), intent(in) :: snowhland(pcols)
    logical, intent(in) :: ghg_chem
    integer, intent(in) :: latmapback(pcols)
    character(len=*), intent(in) :: chem_name
    real(r8), intent(out) :: drydepflx(pcols,pcnst)
    real(r8), intent(inout) :: cflx(pcols,pcnst)
    real(r8), intent(inout) :: qtend(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)

    call t_startf('ap_gas_phase_chemdr_run')
    call gas_phase_chemdr_run(lchnk, ncol, pcols, pver, pverp, pcnst, &
         imozart, q, phis, zm, zi, &
         calday, tfld, pmid, pdel, pint, cldw, troplev, ncldwtr, ufld, &
         vfld, delt, ps, xactive_prates, fsds, ts, asdir, ocnfrac, &
         icefrac, precc, precl, snowhland, ghg_chem, latmapback, &
         chem_name, drydepflx, cflx, qtend, pbuf)
    call t_stopf('ap_gas_phase_chemdr_run')

  end subroutine gas_phase_chemdr

end module mo_gas_phase_chemdr

! CAM compatibility adapter for the decoupled Neu wet-deposition process.
module mo_neu_wetdep

  use shr_kind_mod, only : r8 => shr_kind_r8
  use constituents, only : pcnst
  use ppgrid, only : pcols, pver, pverp
  use perf_mod, only : t_startf, t_stopf
  use neu_wetdep_scheme, only : neu_wetdep_init, neu_wetdep_tend_run, &
       do_neu_wetdep

  implicit none
  private

  public :: neu_wetdep_init
  public :: neu_wetdep_tend
  public :: do_neu_wetdep

contains

  subroutine neu_wetdep_tend(lchnk, ncol, mmr, pmid, pdel, zint, tfld, delt, &
       prain, nevapr, cld, cmfdqr, wd_tend)

    integer,  intent(in) :: lchnk, ncol
    real(r8), intent(in) :: mmr(pcols,pver,pcnst)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: zint(pcols,pver+1)
    real(r8), intent(in) :: tfld(pcols,pver)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: prain(ncol,pver)
    real(r8), intent(in) :: nevapr(ncol,pver)
    real(r8), intent(in) :: cld(ncol,pver)
    real(r8), intent(in) :: cmfdqr(ncol,pver)
    real(r8), intent(inout) :: wd_tend(pcols,pver,pcnst)

    call t_startf('ap_neu_wetdep_tend_run')
    call neu_wetdep_tend_run(lchnk, ncol, pcols, pver, pverp, pcnst, mmr, &
         pmid, pdel, zint, tfld, delt, prain, nevapr, cld, cmfdqr, wd_tend)
    call t_stopf('ap_neu_wetdep_tend_run')

  end subroutine neu_wetdep_tend

end module mo_neu_wetdep

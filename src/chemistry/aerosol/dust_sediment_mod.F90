module dust_sediment_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, pverp
  use physconst, only: gravit
  use perf_mod, only: t_startf, t_stopf
  use cam_abortutils, only: endrun
  use ap_dust_sediment_tend_scheme, only: dust_sediment_vel, &
       dust_sediment_tend_run

  implicit none
  private

  public :: dust_sediment_vel
  public :: dust_sediment_tend

contains

  subroutine dust_sediment_tend(ncol, dtime, pint, pmid, pdel, t, &
       dustmr, pvdust, dusttend, sfdust)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: dtime
    real(r8), intent(in) :: pint(pcols,pverp), pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver), t(pcols,pver)
    real(r8), intent(in) :: dustmr(pcols,pver), pvdust(pcols,pverp)
    real(r8), intent(out) :: dusttend(pcols,pver), sfdust(pcols)
    character(len=512) :: errmsg
    integer :: errflg

    call t_startf('ap_dust_sediment_tend_run')
    call dust_sediment_tend_run(ncol, dtime, pint, pmid, pdel, t, &
         dustmr, pvdust, gravit, dusttend, sfdust, errmsg, errflg)
    call t_stopf('ap_dust_sediment_tend_run')
    if (errflg /= 0) call endrun(trim(errmsg))

  end subroutine dust_sediment_tend

end module dust_sediment_mod

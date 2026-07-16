module ap_compute_tms_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: compute_tms_run

  real(r8), parameter :: horomin = 1._r8
  real(r8), parameter :: z0max   = 100._r8
  real(r8), parameter :: dv2min  = 0.01_r8

contains

  !> \section arg_table_compute_tms_run Argument Table
  !! \htmlinclude compute_tms_run.html
  !!
  subroutine compute_tms_run(pcols, pver, ncol, orocnst, z0fac, karman, &
       gravit, rair, u, v, t, pmid, exner, zm, sgh, landfrac, ksrf, &
       taux, tauy, errmsg, errflg)

    integer,  intent(in)  :: pcols
    integer,  intent(in)  :: pver
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: orocnst
    real(r8), intent(in)  :: z0fac
    real(r8), intent(in)  :: karman
    real(r8), intent(in)  :: gravit
    real(r8), intent(in)  :: rair
    real(r8), intent(in)  :: u(pcols,pver)
    real(r8), intent(in)  :: v(pcols,pver)
    real(r8), intent(in)  :: t(pcols,pver)
    real(r8), intent(in)  :: pmid(pcols,pver)
    real(r8), intent(in)  :: exner(pcols,pver)
    real(r8), intent(in)  :: zm(pcols,pver)
    real(r8), intent(in)  :: sgh(pcols)
    real(r8), intent(in)  :: landfrac(pcols)
    real(r8), intent(out) :: ksrf(pcols)
    real(r8), intent(out) :: taux(pcols)
    real(r8), intent(out) :: tauy(pcols)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer  :: i
    integer  :: kb, kt
    real(r8) :: horo
    real(r8) :: z0oro
    real(r8) :: dv2
    real(r8) :: ri
    real(r8) :: stabfri
    real(r8) :: rho
    real(r8) :: cd
    real(r8) :: vmag

    errmsg = ''
    errflg = 0

    do i = 1, ncol
       horo = orocnst * sgh(i)

       if (horo < horomin) then
          ksrf(i) = 0._r8
          taux(i) = 0._r8
          tauy(i) = 0._r8
       else
          z0oro = min(z0fac * horo, z0max)
          cd = (karman / log((zm(i,pver) + z0oro) / z0oro))**2

          kt  = pver - 1
          kb  = pver
          dv2 = max((u(i,kt) - u(i,kb))**2 + &
               (v(i,kt) - v(i,kb))**2, dv2min)

          ri = 2._r8 * gravit * &
               (t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb)) * &
               (zm(i,kt) - zm(i,kb)) / &
               ((t(i,kt) * exner(i,kt) + t(i,kb) * exner(i,kb)) * dv2)

          stabfri = max(0._r8, min(1._r8, 1._r8 - ri))
          cd      = cd * stabfri

          rho     = pmid(i,pver) / (rair * t(i,pver))
          vmag    = sqrt(u(i,pver)**2 + v(i,pver)**2)
          ksrf(i) = rho * cd * vmag * landfrac(i)
          taux(i) = -ksrf(i) * u(i,pver)
          tauy(i) = -ksrf(i) * v(i,pver)
       end if
    end do

  end subroutine compute_tms_run

end module ap_compute_tms_scheme

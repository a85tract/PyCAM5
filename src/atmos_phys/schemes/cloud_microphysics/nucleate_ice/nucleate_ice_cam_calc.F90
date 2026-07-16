module ap_nucleate_ice_cam_calc_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8
  use nucleate_ice, only: nucleati

  implicit none
  private

  public :: nucleate_ice_cam_calc_run

contains

  subroutine nucleate_ice_cam_calc_run( &
       pcols, pver, ncol, top_lev, tmelt, preexisting_ice_enabled, &
       wsubi, t, pmid, relhum, icldm, qc, qi, ni, rho, &
       so4_num, dst_num, soot_num, &
       naai, naai_hom, nihf, niimm, nidep, nimey, &
       fhom, wice, weff, innso4, innbc, inndust, inhet, inhom, &
       infrehom, infrein)

    integer, intent(in) :: pcols, pver, ncol, top_lev
    real(r8), intent(in) :: tmelt
    logical, intent(in) :: preexisting_ice_enabled
    real(r8), intent(in) :: wsubi(pcols,pver), t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver), relhum(pcols,pver)
    real(r8), intent(in) :: icldm(pcols,pver), qc(pcols,pver)
    real(r8), intent(in) :: qi(pcols,pver), ni(pcols,pver)
    real(r8), intent(in) :: rho(pcols,pver)
    real(r8), intent(in) :: so4_num(pcols,pver)
    real(r8), intent(in) :: dst_num(pcols,pver), soot_num(pcols,pver)
    real(r8), intent(out) :: naai(pcols,pver), naai_hom(pcols,pver)
    real(r8), intent(out) :: nihf(pcols,pver), niimm(pcols,pver)
    real(r8), intent(out) :: nidep(pcols,pver), nimey(pcols,pver)
    real(r8), intent(out) :: fhom(pcols,pver), wice(pcols,pver)
    real(r8), intent(out) :: weff(pcols,pver), innso4(pcols,pver)
    real(r8), intent(out) :: innbc(pcols,pver), inndust(pcols,pver)
    real(r8), intent(out) :: inhet(pcols,pver), inhom(pcols,pver)
    real(r8), intent(out) :: infrehom(pcols,pver), infrein(pcols,pver)

    integer :: i, k

    naai = 0._r8
    naai_hom = 0._r8
    nihf = 0._r8
    niimm = 0._r8
    nidep = 0._r8
    nimey = 0._r8
    fhom = 0._r8
    wice = 0._r8
    weff = 0._r8
    innso4 = 0._r8
    innbc = 0._r8
    inndust = 0._r8
    inhet = 0._r8
    inhom = 0._r8
    infrehom = 0._r8
    infrein = 0._r8

    do k = top_lev, pver
       do i = 1, ncol
          if (t(i,k) < tmelt - 5._r8) then
             call nucleati( &
                  wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k), &
                  qc(i,k), qi(i,k), ni(i,k), rho(i,k), &
                  so4_num(i,k), dst_num(i,k), soot_num(i,k), &
                  naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                  wice(i,k), weff(i,k), fhom(i,k))

             naai_hom(i,k) = nihf(i,k)

             nihf(i,k) = nihf(i,k) * rho(i,k)
             niimm(i,k) = niimm(i,k) * rho(i,k)
             nidep(i,k) = nidep(i,k) * rho(i,k)
             nimey(i,k) = nimey(i,k) * rho(i,k)

             if (preexisting_ice_enabled) then
                innso4(i,k) = so4_num(i,k) * 1.e6_r8
                innbc(i,k) = soot_num(i,k) * 1.e6_r8
                inndust(i,k) = dst_num(i,k) * 1.e6_r8
                infrein(i,k) = 1._r8
                inhet(i,k) = niimm(i,k) + nidep(i,k)
                inhom(i,k) = nihf(i,k)
                if (inhom(i,k) > 1.e3_r8) infrehom(i,k) = 1._r8

                if (infrehom(i,k) < 0.5_r8 .and. inhet(i,k) < 1._r8) then
                   innso4(i,k) = 0._r8
                   innbc(i,k) = 0._r8
                   inndust(i,k) = 0._r8
                   infrein(i,k) = 0._r8
                   inhet(i,k) = 0._r8
                   inhom(i,k) = 0._r8
                   infrehom(i,k) = 0._r8
                   wice(i,k) = 0._r8
                   weff(i,k) = 0._r8
                   fhom(i,k) = 0._r8
                end if
             end if
          end if
       end do
    end do

  end subroutine nucleate_ice_cam_calc_run

end module ap_nucleate_ice_cam_calc_scheme

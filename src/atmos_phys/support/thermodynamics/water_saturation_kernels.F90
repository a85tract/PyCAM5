! Intrinsic-only water saturation kernels matching CAM wv_sat_methods.
module ap_water_saturation_kernels

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  integer, parameter, public :: ap_old_goff_gratch_idx = 0
  integer, parameter, public :: ap_goff_gratch_idx = 1
  integer, parameter, public :: ap_murphy_koop_idx = 2
  integer, parameter, public :: ap_bolton_idx = 3

  public :: ap_qsat_water

contains

  elemental subroutine ap_qsat_water(t, p, scheme_idx, epsilo, tmelt, &
       tboil, es, qs)
    real(r8), intent(in) :: t, p
    integer, intent(in) :: scheme_idx
    real(r8), intent(in) :: epsilo, tmelt, tboil
    real(r8), intent(out) :: es, qs
    real(r8) :: omeps

    select case (scheme_idx)
    case (ap_goff_gratch_idx)
       es = 10._r8**(-7.90298_r8*(tboil/t-1._r8) +             &
            5.02808_r8*log10(tboil/t) -                         &
            1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil)) &
            -1._r8) + 8.1328e-3_r8*                            &
            (10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8) +    &
            log10(1013.246_r8))*100._r8
    case (ap_murphy_koop_idx)
       es = exp(54.842763_r8 - (6763.22_r8/t) -                 &
            (4.210_r8*log(t)) + (0.000367_r8*t) +              &
            (tanh(0.0415_r8*(t-218.8_r8)) *                    &
            (53.878_r8-(1331.22_r8/t)-(9.44523_r8*log(t)) +    &
            0.014025_r8*t)))
    case (ap_old_goff_gratch_idx)
       es = old_goff_gratch_svp_water(t, tboil)
    case (ap_bolton_idx)
       es = 611.2_r8*exp((17.67_r8*(t-tmelt))/                  &
            ((t-tmelt)+243.5_r8))
    case default
       es = 0._r8
    end select

    omeps = 1._r8 - epsilo
    if ((p-es) <= 0._r8) then
       qs = 1._r8
    else
       qs = epsilo*es/(p-omeps*es)
    end if
    es = min(es, p)
  end subroutine ap_qsat_water

  elemental function old_goff_gratch_svp_water(t, tboil) result(es)
    real(r8), intent(in) :: t, tboil
    real(r8) :: es
    real(r8) :: ps, e1, e2, f1, f2, f3, f4, f5, f

    ps = 1013.246_r8
    e1 = 11.344_r8*(1._r8-t/tboil)
    e2 = -3.49149_r8*(tboil/t-1._r8)
    f1 = -7.90298_r8*(tboil/t-1._r8)
    f2 = 5.02808_r8*log10(tboil/t)
    f3 = -1.3816_r8*(10._r8**e1-1._r8)/10000000._r8
    f4 = 8.1328_r8*(10._r8**e2-1._r8)/1000._r8
    f5 = log10(ps)
    f = f1+f2+f3+f4+f5
    es = (10._r8**f)*100._r8
  end function old_goff_gratch_svp_water

end module ap_water_saturation_kernels

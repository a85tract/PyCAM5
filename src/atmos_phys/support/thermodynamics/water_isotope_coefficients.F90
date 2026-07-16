module ap_water_isotope_coefficients

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_TKTRIP

  implicit none
  private
  save

  public :: wiso_alpl
  public :: wiso_alpi
  public :: wiso_akel
  public :: wiso_akci

  integer, parameter :: isph2o = 1
  integer, parameter :: isphdo = 3

  real(r8), parameter :: dkfac = 0.58_r8
  real(r8), parameter :: tkini = 253.15_r8
  real(r8), parameter :: fsata = 1.000_r8
  real(r8), parameter :: fsatb = -0.002_r8
  real(r8), parameter :: ssatmx = 2.00_r8
  real(r8), parameter :: fkhum = 0.25_r8
  real(r8), parameter :: tzero = SHR_CONST_TKTRIP

  real(r8), parameter :: difrm(4) = &
       (/ 1._r8, 1._r8, 0.9757_r8, 0.9727_r8 /)
  real(r8), parameter :: alpal(4) = &
       (/ 0._r8, 0._r8, 1158.8e-12_r8, 0.35041e+6_r8 /)
  real(r8), parameter :: alpbl(4) = &
       (/ 0._r8, 0._r8, -1620.1e-9_r8, -1.6664e+3_r8 /)
  real(r8), parameter :: alpcl(4) = &
       (/ 0._r8, 0._r8, 794.84e-6_r8, 6.7123_r8 /)
  real(r8), parameter :: alpdl(4) = &
       (/ 0._r8, 0._r8, -161.04e-3_r8, -7.685e-3_r8 /)
  real(r8), parameter :: alpel(4) = &
       (/ 0._r8, 0._r8, 2.9992e+6_r8, 0._r8 /)
  real(r8), parameter :: alpai(4) = &
       (/ 0._r8, 0._r8, 16289._r8, 0._r8 /)
  real(r8), parameter :: alpbi(4) = &
       (/ 0._r8, 0._r8, 0._r8, 11.839_r8 /)
  real(r8), parameter :: alpci(4) = &
       (/ 0._r8, 0._r8, -9.45e-2_r8, -28.224e-3_r8 /)

contains

  function wiso_alpl(isp, tk)
    integer, intent(in) :: isp
    real(r8), intent(in) :: tk
    real(r8) :: wiso_alpl

    if (isp == isph2o) then
       wiso_alpl = 1._r8
       return
    end if

    if (isp == isphdo) then
       wiso_alpl = exp(alpal(isp)*tk**3 + alpbl(isp)*tk**2 + &
            alpcl(isp)*tk + alpdl(isp) + alpel(isp)/tk**3)
    else
       wiso_alpl = exp(alpal(isp)/tk**3 + alpbl(isp)/tk**2 + &
            alpcl(isp)/tk + alpdl(isp))
    end if
  end function wiso_alpl

  function wiso_alpi(isp, tk)
    integer, intent(in) :: isp
    real(r8), intent(in) :: tk
    real(r8) :: wiso_alpi

    if (isp == isph2o) then
       wiso_alpi = 1._r8
       return
    end if

    wiso_alpi = exp(alpai(isp)/tk**2 + alpbi(isp)/tk + alpci(isp))
  end function wiso_alpi

  function wiso_akel(isp, tk, hum0, alpeq)
    integer, intent(in) :: isp
    real(r8), intent(in) :: tk
    real(r8), intent(in) :: hum0
    real(r8), intent(in) :: alpeq
    real(r8) :: wiso_akel
    real(r8) :: h0
    real(r8) :: heff
    real(r8) :: difrmj
    real(r8) :: dondi

    h0 = min(1.0_r8, hum0)
    difrmj = difrm(isp)
    heff = wiso_heff(h0)
    dondi = (1/difrmj)**dkfac
    wiso_akel = alpeq*heff / (alpeq*dondi*(heff-1._r8) + 1._r8)
  end function wiso_akel

  function wiso_akci(isp, tk, alpeq)
    integer, intent(in) :: isp
    real(r8), intent(in) :: tk
    real(r8), intent(in) :: alpeq
    real(r8) :: wiso_akci
    real(r8) :: sat1
    real(r8) :: difrmj
    real(r8) :: dondi

    if (tk < tkini) then
       sat1 = max(1._r8, wiso_ssatf(tk))
       difrmj = difrm(isp)
       dondi = 1._r8/difrmj
       wiso_akci = alpeq*sat1 / (alpeq*dondi*(sat1-1._r8) + 1._r8)
    else
       wiso_akci = alpeq
    end if
  end function wiso_akci

  function wiso_heff(h0)
    real(r8), intent(in) :: h0
    real(r8) :: wiso_heff

    wiso_heff = min(1.0_r8, fkhum*h0 + 1.0_r8-fkhum)
  end function wiso_heff

  function wiso_ssatf(tk)
    real(r8), intent(in) :: tk
    real(r8) :: wiso_ssatf

    wiso_ssatf = fsata + fsatb*(tk-tzero)
    wiso_ssatf = max(wiso_ssatf, 1.0_r8)
    wiso_ssatf = min(wiso_ssatf, ssatmx)
  end function wiso_ssatf

end module ap_water_isotope_coefficients

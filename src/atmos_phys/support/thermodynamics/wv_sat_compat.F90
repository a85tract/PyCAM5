module ap_wv_sat_compat
  use shr_kind_mod, only: r8 => shr_kind_r8
  use wv_sat_methods, only: wv_sat_qsat_water

  implicit none
  private

  public :: qsat_water

contains

elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)
  real(r8), intent(in) :: t
  real(r8), intent(in) :: p
  real(r8), intent(out) :: es
  real(r8), intent(out) :: qs
  real(r8), intent(out), optional :: gam
  real(r8), intent(out), optional :: dqsdt
  real(r8), intent(out), optional :: enthalpy

  call wv_sat_qsat_water(t, p, es, qs)

  ! ZM only requests es and qs.  Keep the legacy qsat_water call shape so
  ! that removing CAM's wv_saturation module does not alter the science
  ! kernel's floating-point compilation boundary.
  if (present(gam)) gam = 0._r8
  if (present(dqsdt)) dqsdt = 0._r8
  if (present(enthalpy)) enthalpy = 0._r8
end subroutine qsat_water

end module ap_wv_sat_compat

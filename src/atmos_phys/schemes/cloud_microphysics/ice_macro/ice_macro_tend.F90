module ap_ice_macro_tend_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: ice_macro_tend_run

contains

! Saturation adjustment for ice
! Add ice mass if supersaturated
!> \section arg_table_ice_macro_tend_run Argument Table
!! \htmlinclude ice_macro_tend_run.html
elemental subroutine ice_macro_tend_run(naai,t,qsi,qv,qi,ni,xxls,deltat, &
     stend,qvtend,qitend,nitend)

  real(r8), intent(in)  :: naai   !Activated number of ice nuclei
  real(r8), intent(in)  :: t      !temperature (k)
  real(r8), intent(in)  :: qsi    !saturation specific humidity over ice
  real(r8), intent(in)  :: qv     !water vapor mixing ratio
  real(r8), intent(in)  :: qi     !ice mixing ratio
  real(r8), intent(in)  :: ni     !ice number concentration
  real(r8), intent(in)  :: xxls   !latent heat of sublimation
  real(r8), intent(in)  :: deltat !timestep
  real(r8), intent(out) :: stend  ! 'temperature' tendency
  real(r8), intent(out) :: qvtend !vapor tendency
  real(r8), intent(out) :: qitend !ice mass tendency
  real(r8), intent(out) :: nitend !ice number tendency

  real(r8) :: tau
  logical  :: tau_constant

  tau_constant = .true.

  stend = 0._r8
  qvtend = 0._r8
  qitend = 0._r8
  nitend = 0._r8

  if (naai.gt.1.e-18_r8.and.qv.gt.qsi) then

     !optional timescale on condensation
     !tau in sections. Try 300. or tau = f(T): 300s  t> 268, 1800s for t<238
     !
     if (.not. tau_constant) then
        if( t.gt. 268.15_r8 ) then
           tau = 300.0_r8
        elseif(t.lt.238.15_r8 ) then
           tau = 1800._r8
        else
           tau = 300._r8 + (1800._r8 - 300._r8) * ( 268.15_r8 - t ) / 30._r8
        endif
     else
         tau = 300._r8
     end if

     qitend = (qv-qsi)/deltat !* exp(-tau/deltat)
     qvtend = 0._r8 - qitend
     stend  = qitend * xxls    ! moist static energy tend...[J/kg/s] !

     ! kg(h2o)/kg(air)/s * J/kg(h2o)  = J/kg(air)/s (=W/kg)
     ! if ice exists (more than 1 L-1) and there is condensation, do not add to number (= growth), else, add 10um ice

     if (ni.lt.1.e3_r8.and.(qi+qitend*deltat).gt.1e-18_r8) then
        nitend = nitend + 3._r8 * qitend/(4._r8*3.14_r8* 10.e-6_r8**3*997._r8)
     endif

  endif

end subroutine ice_macro_tend_run

end module ap_ice_macro_tend_scheme

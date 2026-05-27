module wv_saturation

!--------------------------------------------------------------------!
! Module Overview:                                                   !
!                                                                    !
! This module provides an interface to wv_sat_methods, providing     !
! saturation vapor pressure and related calculations to CAM.         !
!                                                                    !
! The original wv_saturation codes were introduced by J. J. Hack,    !
! February 1990. The code has been extensively rewritten since then, !
! including a total refactoring in Summer 2012.                      !
!                                                                    !
!--------------------------------------------------------------------!
! Methods:                                                           !
!                                                                    !
! Pure water/ice saturation vapor pressures are calculated on the    !
! fly, with the specific method determined by a runtime option.      !
! Mixed phase SVP is interpolated from the internal table, estbl,    !
! which is created during initialization.                            !
!                                                                    !
! The default method for calculating SVP is determined by a namelist !
! option, and used whenever svp_water/ice or qsat are called.        !
!                                                                    !
!--------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8
use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
use physconst,    only: epsilo, &
                        latvap, &
                        latice, &
                        rh2o,   &
                        cpair,  &
                        tmelt,  &
                        h2otrip

use wv_sat_methods, only: &
     svp_to_qsat => wv_sat_svp_to_qsat

implicit none
private
save

! Public interfaces
! Namelist, initialization, finalization
public wv_sat_readnl
public wv_sat_init
public wv_sat_final

! Saturation vapor pressure calculations
public svp_water
public svp_ice
  
! Mixed phase (water + ice) saturation vapor pressure table lookup
public estblf

public svp_to_qsat

! Subroutines that return both SVP and humidity
! Optional arguments do temperature derivatives
public qsat           ! Mixed phase
public qsat_water     ! SVP over water only
public qsat_ice       ! SVP over ice only

! Wet bulb temperature solver
public :: findsp_vc, findsp

! Data

! This value is slightly high, but it seems to be the value for the
! steam point of water originally (and most frequently) used in the
! Goff & Gratch scheme.
real(r8), parameter :: tboil = 373.16_r8

! Table of saturation vapor pressure values (estbl) from tmin to
! tmax+1 Kelvin, in one degree increments.  ttrice defines the
! transition region, estbl contains a combination of ice & water
! values.
! Make these public parameters in case another module wants to see the
! extent of the table.
  real(r8), public, parameter :: tmin = 127.16_r8
  real(r8), public, parameter :: tmax = 375.16_r8

  real(r8), parameter :: ttrice = 20.00_r8  ! transition range from es over H2O to es over ice

  integer :: plenest                             ! length of estbl
  real(r8), allocatable, target :: estbl(:)      ! table values of saturation vapor pressure

  real(r8) :: omeps      ! 1.0_r8 - epsilo

  real(r8) :: c3         ! parameter used by findsp
  logical :: use_native_wv_saturation_impl = .false.
  logical :: wv_saturation_impl_selected = .false.
  logical :: wv_saturation_proof_written = .false.
  logical :: wv_sat_readnl_logged = .false.
  logical :: wv_sat_init_logged = .false.
  logical :: wv_sat_final_logged = .false.
  logical :: tq_enthalpy_logged = .false.
  logical :: findsp_logged = .false.
  logical :: findsp_vc_logged = .false.

  interface
     function wv_sat_readnl_codon() result(out_c) bind(c, name="wv_sat_readnl_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t) :: out_c
     end function wv_sat_readnl_codon

     function wv_saturation_value_codon(value_c) result(value_out) &
          bind(c, name="wv_saturation_value_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: value_c
       real(c_double) :: value_out
     end function wv_saturation_value_codon

     pure function wv_saturation_limit_es_codon(es_c, p_c) result(es_out) &
          bind(c, name="wv_saturation_limit_es_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: es_c, p_c
       real(c_double) :: es_out
     end function wv_saturation_limit_es_codon

     pure function wv_saturation_tq_enthalpy_codon(cpair_c, t_c, q_c, hltalt_c) result(enthalpy_out) &
          bind(c, name="wv_saturation_tq_enthalpy_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: cpair_c, t_c, q_c, hltalt_c
       real(c_double) :: enthalpy_out
     end function wv_saturation_tq_enthalpy_codon

     function wv_saturation_touch_codon(stage_c) result(stage_out) bind(c, name="wv_saturation_touch_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function wv_saturation_touch_codon

     pure function wv_saturation_no_ip_hltalt_codon(t_c, tmelt_c, latvap_c) result(hltalt_out) &
          bind(c, name="wv_saturation_no_ip_hltalt_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: t_c, tmelt_c, latvap_c
       real(c_double) :: hltalt_out
     end function wv_saturation_no_ip_hltalt_codon

     pure subroutine wv_saturation_calc_hltalt_codon(t_c, tmelt_c, ttrice_c, latvap_c, latice_c, &
          pcf1_c, pcf2_c, pcf3_c, pcf4_c, pcf5_c, hltalt_c, tterm_c) &
          bind(c, name="wv_saturation_calc_hltalt_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: t_c, tmelt_c, ttrice_c, latvap_c, latice_c
       real(c_double), value :: pcf1_c, pcf2_c, pcf3_c, pcf4_c, pcf5_c
       real(c_double), intent(out) :: hltalt_c, tterm_c
     end subroutine wv_saturation_calc_hltalt_codon

     pure function wv_saturation_deriv_dqsdt_codon(t_c, p_c, es_c, qs_c, hltalt_c, tterm_c, &
          rh2o_c, omeps_c) result(dqsdt_out) bind(c, name="wv_saturation_deriv_dqsdt_codon")
       use iso_c_binding, only: c_double
       real(c_double), value :: t_c, p_c, es_c, qs_c, hltalt_c, tterm_c, rh2o_c, omeps_c
       real(c_double) :: dqsdt_out
     end function wv_saturation_deriv_dqsdt_codon

     pure function svp_water_codon(t_c, idx_c) result(es_c) bind(c, name="svp_water_codon")
       use iso_c_binding, only: c_double, c_int64_t
       real(c_double), value :: t_c
       integer(c_int64_t), value :: idx_c
       real(c_double) :: es_c
     end function svp_water_codon

     pure function svp_ice_codon(t_c, idx_c) result(es_c) bind(c, name="svp_ice_codon")
       use iso_c_binding, only: c_double, c_int64_t
       real(c_double), value :: t_c
       integer(c_int64_t), value :: idx_c
       real(c_double) :: es_c
     end function svp_ice_codon

     pure function svp_trans_codon(t_c, idx_c, tmelt_c, ttrice_c) result(es_c) &
          bind(c, name="svp_trans_codon")
       use iso_c_binding, only: c_double, c_int64_t
       real(c_double), value :: t_c, tmelt_c, ttrice_c
       integer(c_int64_t), value :: idx_c
       real(c_double) :: es_c
     end function svp_trans_codon

     pure function estblf_codon(t_c, tmin_c, tmax_c, estbl_p, plenest_c) result(es_c) &
          bind(c, name="estblf_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       real(c_double), value :: t_c, tmin_c, tmax_c
       type(c_ptr), value :: estbl_p
       integer(c_int64_t), value :: plenest_c
       real(c_double) :: es_c
     end function estblf_codon

     pure subroutine wv_saturation_findsp_codon(q_c, t_c, p_c, use_ice_c, idx_c, &
          epsilo_c, omeps_c, cpair_c, c3_c, tmelt_c, ttrice_c, latvap_c, latice_c, &
          rh2o_c, pcf1_c, pcf2_c, pcf3_c, pcf4_c, pcf5_c, tmin_c, tmax_c, estbl_p, &
          plenest_c, tsp_p, qsp_p, status_p) bind(c, name="wv_saturation_findsp_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       real(c_double), value :: q_c, t_c, p_c
       integer(c_int64_t), value :: use_ice_c, idx_c, plenest_c
       real(c_double), value :: epsilo_c, omeps_c, cpair_c, c3_c, tmelt_c, ttrice_c
       real(c_double), value :: latvap_c, latice_c, rh2o_c
       real(c_double), value :: pcf1_c, pcf2_c, pcf3_c, pcf4_c, pcf5_c, tmin_c, tmax_c
       type(c_ptr), value :: estbl_p, tsp_p, qsp_p, status_p
     end subroutine wv_saturation_findsp_codon

     subroutine wv_sat_init_table_codon(plenest_c, tmin_c, tmelt_c, ttrice_c, idx_c, estbl_p) &
          bind(c, name="wv_sat_init_table_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: plenest_c, idx_c
       real(c_double), value :: tmin_c, tmelt_c, ttrice_c
       type(c_ptr), value :: estbl_p
     end subroutine wv_sat_init_table_codon

     subroutine findsp_vc_codon(n_c, use_ice_c, idx_c, epsilo_c, omeps_c, cpair_c, c3_c, &
          tmelt_c, ttrice_c, latvap_c, latice_c, rh2o_c, pcf1_c, pcf2_c, pcf3_c, pcf4_c, &
          pcf5_c, tmin_c, tmax_c, estbl_p, plenest_c, q_p, t_p, p_p, tsp_p, qsp_p, status_p) &
          bind(c, name="findsp_vc_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: n_c, use_ice_c, idx_c, plenest_c
       real(c_double), value :: epsilo_c, omeps_c, cpair_c, c3_c, tmelt_c, ttrice_c
       real(c_double), value :: latvap_c, latice_c, rh2o_c
       real(c_double), value :: pcf1_c, pcf2_c, pcf3_c, pcf4_c, pcf5_c, tmin_c, tmax_c
       type(c_ptr), value :: estbl_p, q_p, t_p, p_p, tsp_p, qsp_p, status_p
     end subroutine findsp_vc_codon
  end interface

  ! Set coefficients for polynomial approximation of difference
  ! between saturation vapor press over water and saturation pressure
  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  ! valid in the range -40 < t < 0 (degrees C).
  real(r8) :: pcf(5) = (/ &
       5.04469588506e-01_r8, &
       -5.47288442819e+00_r8, &
       -3.67471858735e-01_r8, &
       -8.95963532403e-03_r8, &
       -7.78053686625e-05_r8 /)

!   --- Degree 6 approximation ---
!  real(r8) :: pcf(6) = (/ &
!       7.63285250063e-02, &
!       5.86048427932e+00, &
!       4.38660831780e-01, &
!       1.37898276415e-02, &
!       2.14444472424e-04, &
!       1.36639103771e-06 /)

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

subroutine wv_saturation_select_impl()

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wv_saturation_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('WV_SATURATION_IMPL', value=impl_name, length=n, &
       status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wv_saturation_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wv_saturation_impl = .false.
  end if

  wv_saturation_impl_selected = .true.

  if (masterproc) then
     if (use_native_wv_saturation_impl) then
        write(iulog,*) 'wv_saturation implementation = native'
     else
        write(iulog,*) 'wv_saturation implementation = codon'
     end if
  end if

end subroutine wv_saturation_select_impl

subroutine wv_saturation_proof_once()

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  if (wv_saturation_proof_written) return
  wv_saturation_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'wv_saturation entered (init scalar helpers = codon)'
  end if

end subroutine wv_saturation_proof_once

real(r8) function wv_saturation_value(value) result(out)

  use iso_c_binding, only: c_double

  real(r8), intent(in) :: value

  call wv_saturation_select_impl()

  if (use_native_wv_saturation_impl) then
     out = value
     return
  end if

  call wv_saturation_proof_once()
  out = real(wv_saturation_value_codon(real(value, c_double)), r8)

end function wv_saturation_value

subroutine wv_saturation_log_direct(logged, proof_line)

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
  end if

end subroutine wv_saturation_log_direct

subroutine wv_sat_readnl(nlfile)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Get runtime options for wv_saturation.                         !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_get_scheme_idx, &
                            wv_sat_valid_idx, &
                            wv_sat_set_default

  use spmd_utils,      only: masterproc
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand
  use cam_abortutils,  only: endrun

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   
  ! Local variables
  integer :: unitn, ierr
  integer(c_int64_t) :: out_c
  logical :: group_found

  character(len=32) :: wv_sat_scheme = "GoffGratch"

  character(len=*), parameter :: subname = 'wv_sat_readnl'

  namelist /wv_sat_nl/ wv_sat_scheme
  !-----------------------------------------------------------------------------

  call wv_saturation_select_impl()
  if (.not. use_native_wv_saturation_impl) then
     group_found = .false.
     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'wv_sat_nl', status=ierr)
        group_found = ierr == 0
        close(unitn)
        call freeunit(unitn)
     end if
#ifdef SPMD
     call mpibcast(group_found, 1, mpilog, 0, mpicom)
#endif
     if (.not. group_found) then
        call wv_saturation_proof_once()
        out_c = wv_sat_readnl_codon()
        if (out_c == 0_c_int64_t) then
           call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
           return
        end if
        if (.not. wv_sat_set_default(wv_sat_scheme)) then
           call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
           return
        end if
        call wv_saturation_log_direct(wv_sat_readnl_logged, 'wv_sat_readnl direct = codon')
        return
     end if
  end if

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'wv_sat_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, wv_sat_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
           return
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  call mpibcast(wv_sat_scheme, len(wv_sat_scheme) , mpichar, 0, mpicom)
#endif

  if (.not. wv_sat_set_default(wv_sat_scheme)) then
     call endrun('wv_sat_readnl :: Invalid wv_sat_scheme.')
     return
  end if

end subroutine wv_sat_readnl

subroutine wv_sat_init
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Initialize module (e.g. setting parameters, initializing the   !
  !   SVP lookup table).                                             !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_methods_init, &
                            wv_sat_get_scheme_idx, &
                            wv_sat_valid_idx, &
                            wv_sat_get_default_idx
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use shr_assert_mod, only: shr_assert_in_domain
  use error_messages, only: handle_errmsg

  integer :: status

  ! For wv_sat_methods error reporting.
  character(len=256) :: errstring

  ! For generating internal SVP table.
  real(r8) :: t         ! Temperature
  integer  :: i         ! Increment counter

  ! Precalculated because so frequently used.
  call wv_saturation_select_impl()
  if (use_native_wv_saturation_impl) then
     omeps = 1.0_r8 - epsilo
  else
     call wv_saturation_proof_once()
     omeps = real(wv_saturation_value_codon(real(1.0_r8 - epsilo, c_double)), r8)
     call wv_saturation_log_direct(wv_sat_init_logged, 'wv_sat_init direct = codon')
  end if

  ! Transition range method is only valid for transition temperatures at:
  ! -40 deg C < T < 0 deg C
  call shr_assert_in_domain(ttrice, ge=0._r8, le=40._r8, varname="ttrice",&
       msg="wv_sat_init: Invalid transition temperature range.")

! This parameter uses a hardcoded 287.04_r8?
  if (use_native_wv_saturation_impl) then
     c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair
  else
     c3 = real(wv_saturation_value_codon(real(287.04_r8*(7.5_r8*log(10._r8))/cpair, c_double)), r8)
  end if

! Init "methods" module containing actual SVP formulae.

  call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, ttrice, &
       epsilo, errstring)

  call handle_errmsg(errstring, subname="wv_sat_methods_init")

  ! Add two to make the table slightly too big, just in case.
  plenest = ceiling(tmax-tmin) + 2

  ! Allocate SVP table.
  allocate(estbl(plenest), stat=status)
  if (status /= 0) then 
     call endrun('wv_sat_init :: ERROR allocating saturation vapor pressure table')
     return
  end if

  if (use_native_wv_saturation_impl) then
     do i = 1, plenest
       estbl(i) = svp_trans(tmin + real(i-1,r8))
     end do
  else
     call wv_sat_init_table_codon(int(plenest, c_int64_t), real(tmin, c_double), &
          real(tmelt, c_double), real(ttrice, c_double), int(wv_sat_get_default_idx(), c_int64_t), &
          c_loc(estbl(1)))
  end if

  if (masterproc) then
     write(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
  end if

end subroutine wv_sat_init

subroutine wv_sat_final
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Deallocate global variables in module.                         !
  !------------------------------------------------------------------!
  use cam_abortutils, only: endrun

  integer :: status
  integer(c_int64_t) :: touch_c

  call wv_saturation_select_impl()
  if (.not. use_native_wv_saturation_impl) then
     call wv_saturation_proof_once()
     touch_c = wv_saturation_touch_codon(3_c_int64_t)
     if (touch_c == 3_c_int64_t) then
        call wv_saturation_log_direct(wv_sat_final_logged, &
             'wv_sat_final direct = codon; finalization selector/touch direct = codon; native deallocate/error island')
     end if
  end if

  if (allocated(estbl)) then

     deallocate(estbl, stat=status)

     if (status /= 0) then
        call endrun('wv_sat_final :: ERROR deallocating table')
        return
     end if

  end if

end subroutine wv_sat_final

!---------------------------------------------------------------------
! DEFAULT SVP FUNCTIONS
!---------------------------------------------------------------------

! Compute saturation vapor pressure over water
elemental function svp_water(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_get_default_idx

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = real(svp_water_codon(real(t, c_double), &
       int(wv_sat_get_default_idx(), c_int64_t)), r8)

end function svp_water

! Compute saturation vapor pressure over ice
elemental function svp_ice(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_get_default_idx

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = real(svp_ice_codon(real(t, c_double), &
       int(wv_sat_get_default_idx(), c_int64_t)), r8)

end function svp_ice

! Compute saturation vapor pressure with an ice-water transition
elemental function svp_trans(t) result(es)

  use wv_sat_methods, only: &
       wv_sat_get_default_idx

  real(r8), intent(in) :: t ! Temperature (K)
  real(r8) :: es            ! SVP (Pa)

  es = real(svp_trans_codon(real(t, c_double), int(wv_sat_get_default_idx(), c_int64_t), &
       real(tmelt, c_double), real(ttrice, c_double)), r8)

end function svp_trans

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Does linear interpolation from nearest values found
! in the table (estbl).
elemental function estblf(t) result(es)

  real(r8), intent(in) :: t ! Temperature 
  real(r8) :: es            ! SVP (Pa)

  es = real(estblf_codon(real(t, c_double), real(tmin, c_double), real(tmax, c_double), &
       c_loc(estbl(1)), int(plenest, c_int64_t)), r8)

end function estblf

! Get enthalpy based only on temperature
! and specific humidity.
elemental function tq_enthalpy(t, q, hltalt) result(enthalpy)

  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: q      ! Specific humidity
  real(r8), intent(in) :: hltalt ! Modified hlat for T derivatives

  real(r8) :: enthalpy

  enthalpy = real(wv_saturation_tq_enthalpy_codon(real(cpair, c_double), real(t, c_double), &
       real(q, c_double), real(hltalt, c_double)), r8)
  
end function tq_enthalpy

!---------------------------------------------------------------------
! LATENT HEAT OF VAPORIZATION CORRECTIONS
!---------------------------------------------------------------------

elemental subroutine no_ip_hltalt(t, hltalt)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of pure liquid water at  !
  !   a given temperature.                                           !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), target, intent(out) :: hltalt  ! Appropriately modified hlat

  hltalt = real(wv_saturation_no_ip_hltalt_codon(real(t, c_double), real(tmelt, c_double), &
       real(latvap, c_double)), r8)

end subroutine no_ip_hltalt

elemental subroutine calc_hltalt(t, hltalt, tterm)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate latent heat of vaporization of water at a given      !
  !   temperature, taking into account the ice phase if temperature  !
  !   is below freezing.                                             !
  !   Optional argument also calculates a term used to calculate     !
  !   d(es)/dT within the water-ice transition range.                !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t        ! Temperature
  ! Outputs
  real(r8), intent(out) :: hltalt  ! Appropriately modified hlat
  ! Term to account for d(es)/dT in transition region.
  real(r8), target, intent(out), optional :: tterm

  ! Local variables
  real(r8) :: tterm_loc

  call wv_saturation_calc_hltalt_codon(real(t, c_double), real(tmelt, c_double), &
       real(ttrice, c_double), real(latvap, c_double), real(latice, c_double), &
       real(pcf(1), c_double), real(pcf(2), c_double), real(pcf(3), c_double), &
       real(pcf(4), c_double), real(pcf(5), c_double), hltalt, tterm_loc)

  if (present(tterm)) tterm = tterm_loc

end subroutine calc_hltalt

!---------------------------------------------------------------------
! OPTIONAL OUTPUTS
!---------------------------------------------------------------------

! Temperature derivative outputs, for qsat_*
elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, &
     gam, dqsdt)

  ! Inputs
  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: p      ! Pressure
  real(r8), intent(in) :: es     ! Saturation vapor pressure
  real(r8), intent(in) :: qs     ! Saturation specific humidity
  real(r8), intent(in) :: hltalt ! Modified latent heat
  real(r8), intent(in) :: tterm  ! Extra term for d(es)/dT in
                                 ! transition region.

  ! Outputs
  real(r8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)

  ! Local variables
  real(r8) :: desdt        ! d(es)/dt
  real(r8) :: dqsdt_loc    ! local copy of dqsdt

  dqsdt_loc = real(wv_saturation_deriv_dqsdt_codon(real(t, c_double), real(p, c_double), &
       real(es, c_double), real(qs, c_double), real(hltalt, c_double), real(tterm, c_double), &
       real(rh2o, c_double), real(omeps, c_double)), r8)

  if (present(dqsdt)) dqsdt = dqsdt_loc
  if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)

end subroutine deriv_outputs

!---------------------------------------------------------------------
! QSAT (SPECIFIC HUMIDITY) PROCEDURES
!---------------------------------------------------------------------

elemental subroutine qsat(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Look up and return saturation vapor pressure from precomputed  !
  !   table, then calculate and return saturation specific humidity. !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives
  real(r8) :: tterm        ! Account for d(es)/dT in transition region

  es = estblf(t)

  qs = svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = real(wv_saturation_limit_es_codon(real(es, c_double), real(p, c_double)), r8)

  ! Calculate optional arguments.
  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call calc_hltalt(t, hltalt, tterm)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     call deriv_outputs(t, p, es, qs, hltalt, tterm, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat

elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_qsat_water

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_water(t, p, es, qs)
  es = real(wv_saturation_limit_es_codon(real(es, c_double), real(p, c_double)), r8)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! "generalized" analytic expression for t derivative of es
     ! accurate to within 1 percent for 173.16 < t < 373.16
     call no_ip_hltalt(t, hltalt)

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_water

elemental subroutine qsat_ice(t, p, es, qs, gam, dqsdt, enthalpy)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  use wv_sat_methods, only: wv_sat_qsat_ice

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q

  ! Local variables
  real(r8) :: hltalt       ! Modified latent heat for T derivatives

  call wv_sat_qsat_ice(t, p, es, qs)
  es = real(wv_saturation_limit_es_codon(real(es, c_double), real(p, c_double)), r8)

  if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

     ! For pure ice, just add latent heats.
     hltalt = latvap + latice

     if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

     ! For pure water/ice transition term is 0.
     call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
          gam=gam, dqsdt=dqsdt)

  end if

end subroutine qsat_ice

!---------------------------------------------------------------------
! FINDSP (WET BULB TEMPERATURE) PROCEDURES
!---------------------------------------------------------------------

subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)

  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun
  use wv_sat_methods,   only: wv_sat_get_default_idx

  ! Wrapper for findsp which is 1D and handles the output status.
  ! Changing findsp to elemental restricted debugging output.
  ! If that output is needed again, it's preferable *not* to copy findsp,
  ! but to change the existing version.

  ! input arguments
  real(r8), target, intent(in) :: q(:)        ! water vapor (kg/kg)
  real(r8), target, intent(in) :: t(:)        ! temperature (K)
  real(r8), target, intent(in) :: p(:)        ! pressure    (Pa)
  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations

  ! output arguments
  real(r8), target, intent(out) :: tsp(:)     ! saturation temp (K)
  real(r8), target, intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)

  integer, target :: status(size(q))   ! flag representing state of output
                               ! 0 => Successful convergence
                               ! 1 => No calculation done: pressure or specific
                               !      humidity not within usable range
                               ! 2 => Run failed to converge
                               ! 4 => Temperature fell below minimum
                               ! 8 => Enthalpy not conserved

  integer :: n, i

  n = size(q)

  call wv_saturation_select_impl()
  if (use_native_wv_saturation_impl) then
     call findsp(q, t, p, use_ice, tsp, qsp, status)
  else
     do i = 1,n
        call findsp(q(i), t(i), p(i), use_ice, tsp(i), qsp(i), status(i))
     end do
     call wv_saturation_log_direct(findsp_logged, &
          'findsp direct = codon; scalar wet-bulb kernel entered from findsp_vc')
     call wv_saturation_log_direct(findsp_vc_logged, &
          'findsp_vc direct = codon; scalar findsp loop/error reporting native island')
  end if

  ! Currently, only 2 and 8 seem to be treated as fatal errors.
  do i = 1,n
     if (status(i) == 2) then
        write(iulog,*) ' findsp not converging at i = ', i
        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
        call endrun ('wv_saturation::FINDSP -- not converging')
     else if (status(i) == 8) then
        write(iulog,*) ' the enthalpy is not conserved at i = ', i
        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
        call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
     endif
  end do

end subroutine findsp_vc

elemental subroutine findsp (q, t, p, use_ice, tsp, qsp, status)
  use wv_sat_methods, only: wv_sat_get_default_idx
!----------------------------------------------------------------------- 
! 
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!

  real(r8), intent(in) :: q        ! water vapor (kg/kg)
  real(r8), intent(in) :: t        ! temperature (K)
  real(r8), intent(in) :: p        ! pressure    (Pa)
  logical,  intent(in) :: use_ice  ! flag to include ice phase in calculations
!
! output arguments
!
  real(r8), target, intent(out) :: tsp      ! saturation temp (K)
  real(r8), target, intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
  integer,  target, intent(out) :: status   ! flag representing state of output
                                    ! 0 => Successful convergence
                                    ! 1 => No calculation done: pressure or specific
                                    !      humidity not within usable range
                                    ! 2 => Run failed to converge
                                    ! 4 => Temperature fell below minimum
                                    ! 8 => Enthalpy not conserved
!
! local variables
!
  integer, parameter :: iter = 8    ! max number of times to iterate the calculation
  integer :: l                      ! iterator

  real(r8) es                   ! sat. vapor pressure
  real(r8) gam                  ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
  real(r8) dgdt                 ! work variable
  real(r8) g                    ! work variable
  real(r8) hltalt               ! lat. heat. of vap.
  real(r8) qs                   ! spec. hum. of water vapor

! work variables
  real(r8) t1, q1, dt, dq
  real(r8) qvd
  real(r8) r1b, c1, c2
  real(r8), parameter :: dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
  real(r8), parameter :: dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
  real(r8) enin, enout

  if (.not. use_native_wv_saturation_impl) then
     call wv_saturation_findsp_codon(real(q, c_double), real(t, c_double), real(p, c_double), &
          int(merge(1, 0, use_ice), c_int64_t), int(wv_sat_get_default_idx(), c_int64_t), &
          real(epsilo, c_double), real(omeps, c_double), real(cpair, c_double), real(c3, c_double), &
          real(tmelt, c_double), real(ttrice, c_double), real(latvap, c_double), real(latice, c_double), &
          real(rh2o, c_double), real(pcf(1), c_double), real(pcf(2), c_double), &
          real(pcf(3), c_double), real(pcf(4), c_double), real(pcf(5), c_double), &
          real(tmin, c_double), real(tmax, c_double), c_loc(estbl(1)), int(plenest, c_int64_t), &
          c_loc(tsp), c_loc(qsp), c_loc(status))
     return
  end if

  ! Saturation specific humidity at this temperature
  if (use_ice) then
     call qsat(t, p, es, qs)
  else
     call qsat_water(t, p, es, qs)
  end if

  ! make sure a meaningful calculation is possible
  if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 &
       .or. t < tmin .or. t > tmax) then
     status = 1
     ! Keep initial parameters when conditions aren't suitable
     tsp = t
     qsp = q
     enin = 1._r8
     enout = 1._r8

     return
  end if

  ! Prepare to iterate
  status = 2

  ! Get initial enthalpy
  if (use_ice) then
     call calc_hltalt(t,hltalt)
  else
     call no_ip_hltalt(t,hltalt)
  end if
  enin = tq_enthalpy(t, q, hltalt)

  ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
  c1 = hltalt*c3
  c2 = (t + 36._r8)**2
  r1b = c2/(c2 + c1*qs)
  qvd = r1b * (q - qs)
  tsp = t + ((hltalt/cpair)*qvd)

  ! Generate qsp, gam, and enout from tsp.
  if (use_ice) then
     call qsat(tsp, p, es, qsp, gam=gam, enthalpy=enout)
  else
     call qsat_water(tsp, p, es, qsp, gam=gam, enthalpy=enout)
  end if

  ! iterate on first guess
  do l = 1, iter

     g = enin - enout
     dgdt = -cpair * (1 + gam)

     ! New tsp
     t1 = tsp - g/dgdt
     dt = abs(t1 - tsp)/t1
     tsp = t1

     ! bail out if past end of temperature range
     if ( tsp < tmin ) then
        tsp = tmin
        ! Get latent heat and set qsp to a value
        ! that preserves enthalpy.
        if (use_ice) then
           call calc_hltalt(tsp,hltalt)
        else
           call no_ip_hltalt(tsp,hltalt)
        end if
        qsp = (enin - cpair*tsp)/hltalt
        enout = tq_enthalpy(tsp, qsp, hltalt)
        status = 4
        exit
     end if

     ! Re-generate qsp, gam, and enout from new tsp.
     if (use_ice) then
        call qsat(tsp, p, es, q1, gam=gam, enthalpy=enout)
     else
        call qsat_water(tsp, p, es, q1, gam=gam, enthalpy=enout)
     end if
     dq = abs(q1 - qsp)/max(q1,1.e-12_r8)
     qsp = q1

     ! if converged at this point, exclude it from more iterations
     if (dt < dttol .and. dq < dqtol) then
        status = 0
        exit
     endif
  end do

  ! Test for enthalpy conservation
  if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8

end subroutine findsp

end module wv_saturation

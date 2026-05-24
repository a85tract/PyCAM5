module wv_sat_methods

! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, <constants>, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.

use cam_logfile, only: iulog
use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
use spmd_utils, only: masterproc

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

real(r8), target :: tmelt   ! Melting point of water at 1 atm (K)
real(r8), target :: h2otrip ! Triple point temperature of water (K)
real(r8), target :: tboil   ! Boiling point of water at 1 atm (K)

real(r8), target :: ttrice  ! Ice-water transition range

real(r8), target :: epsilo  ! Ice-water transition range
real(r8), target :: omeps   ! 1._r8 - epsilo

! Indices representing individual schemes
integer, parameter :: Invalid_idx = -1
integer, parameter :: OldGoffGratch_idx = 0
integer, parameter :: GoffGratch_idx = 1
integer, parameter :: MurphyKoop_idx = 2
integer, parameter :: Bolton_idx = 3

! Index representing the current default scheme.
integer, parameter :: initial_default_idx = GoffGratch_idx
integer, target :: default_idx = initial_default_idx
logical :: use_native_wv_sat_methods_impl = .false.
logical :: wv_sat_methods_impl_selected = .false.
logical :: wv_sat_methods_proof_written = .false.

interface
  function wv_sat_methods_value_codon(value_c) result(value_out) &
       bind(c, name="wv_sat_methods_value_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: value_c
    real(c_double) :: value_out
  end function wv_sat_methods_value_codon

  function wv_sat_methods_omeps_codon(epsilo_c) result(omeps_c) &
       bind(c, name="wv_sat_methods_omeps_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: epsilo_c
    real(c_double) :: omeps_c
  end function wv_sat_methods_omeps_codon

  subroutine wv_sat_methods_init_codon(tmelt_in_c, h2otrip_in_c, tboil_in_c, &
       ttrice_in_c, epsilo_in_c, tmelt_p, h2otrip_p, tboil_p, ttrice_p, epsilo_p, omeps_p) &
       bind(c, name="wv_sat_methods_init_codon")
    use iso_c_binding, only: c_double, c_ptr
    real(c_double), value :: tmelt_in_c, h2otrip_in_c, tboil_in_c
    real(c_double), value :: ttrice_in_c, epsilo_in_c
    type(c_ptr), value :: tmelt_p, h2otrip_p, tboil_p, ttrice_p, epsilo_p, omeps_p
  end subroutine wv_sat_methods_init_codon

  pure function wv_sat_valid_idx_codon(idx_c) result(status_c) &
       bind(c, name="wv_sat_valid_idx_codon")
    use iso_c_binding, only: c_int64_t
    integer(c_int64_t), value :: idx_c
    integer(c_int64_t) :: status_c
  end function wv_sat_valid_idx_codon

  pure function wv_sat_svp_to_qsat_codon(es_c, p_c, epsilo_c, omeps_c) result(qs_c) &
       bind(c, name="wv_sat_svp_to_qsat_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: es_c, p_c, epsilo_c, omeps_c
    real(c_double) :: qs_c
  end function wv_sat_svp_to_qsat_codon

  pure function wv_sat_get_scheme_idx_codon(name_len_c, name_ascii_p) result(idx_c) &
       bind(c, name="wv_sat_get_scheme_idx_codon")
    use iso_c_binding, only: c_int64_t, c_ptr
    integer(c_int64_t), value :: name_len_c
    type(c_ptr), value :: name_ascii_p
    integer(c_int64_t) :: idx_c
  end function wv_sat_get_scheme_idx_codon

  function wv_sat_set_default_codon(tmp_idx_c, default_idx_p) result(status_c) &
       bind(c, name="wv_sat_set_default_codon")
    use iso_c_binding, only: c_int64_t, c_ptr
    integer(c_int64_t), value :: tmp_idx_c
    type(c_ptr), value :: default_idx_p
    integer(c_int64_t) :: status_c
  end function wv_sat_set_default_codon

  pure subroutine wv_sat_qsat_water_codon(t_c, p_c, idx_c, epsilo_c, omeps_c, es_p, qs_p) &
       bind(c, name="wv_sat_qsat_water_codon")
    use iso_c_binding, only: c_double, c_int64_t, c_ptr
    real(c_double), value :: t_c, p_c, epsilo_c, omeps_c
    integer(c_int64_t), value :: idx_c
    type(c_ptr), value :: es_p, qs_p
  end subroutine wv_sat_qsat_water_codon

  pure function wv_sat_svp_water_codon(t_c, idx_c) result(es_c) &
       bind(c, name="wv_sat_svp_water_codon")
    use iso_c_binding, only: c_double, c_int64_t
    real(c_double), value :: t_c
    integer(c_int64_t), value :: idx_c
    real(c_double) :: es_c
  end function wv_sat_svp_water_codon

  pure function wv_sat_svp_ice_codon(t_c, idx_c) result(es_c) &
       bind(c, name="wv_sat_svp_ice_codon")
    use iso_c_binding, only: c_double, c_int64_t
    real(c_double), value :: t_c
    integer(c_int64_t), value :: idx_c
    real(c_double) :: es_c
  end function wv_sat_svp_ice_codon

  pure function wv_sat_svp_trans_codon(t_c, idx_c, tmelt_c, ttrice_c) result(es_c) &
       bind(c, name="wv_sat_svp_trans_codon")
    use iso_c_binding, only: c_double, c_int64_t
    real(c_double), value :: t_c, tmelt_c, ttrice_c
    integer(c_int64_t), value :: idx_c
    real(c_double) :: es_c
  end function wv_sat_svp_trans_codon

  pure function goffgratch_svp_water_codon(t_c) result(es_c) &
       bind(c, name="goffgratch_svp_water_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function goffgratch_svp_water_codon

  pure function goffgratch_svp_ice_codon(t_c) result(es_c) &
       bind(c, name="goffgratch_svp_ice_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function goffgratch_svp_ice_codon

  pure function murphykoop_svp_water_codon(t_c) result(es_c) &
       bind(c, name="murphykoop_svp_water_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function murphykoop_svp_water_codon

  pure function murphykoop_svp_ice_codon(t_c) result(es_c) &
       bind(c, name="murphykoop_svp_ice_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function murphykoop_svp_ice_codon

  pure function oldgoffgratch_svp_water_codon(t_c) result(es_c) &
       bind(c, name="oldgoffgratch_svp_water_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function oldgoffgratch_svp_water_codon

  pure function oldgoffgratch_svp_ice_codon(t_c) result(es_c) &
       bind(c, name="oldgoffgratch_svp_ice_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function oldgoffgratch_svp_ice_codon

  pure function bolton_svp_water_codon(t_c) result(es_c) &
       bind(c, name="bolton_svp_water_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: t_c
    real(c_double) :: es_c
  end function bolton_svp_water_codon
end interface

public wv_sat_methods_init
public wv_sat_get_scheme_idx
public wv_sat_valid_idx
public wv_sat_get_default_idx

public wv_sat_set_default
public wv_sat_reset_default

public wv_sat_svp_water
public wv_sat_svp_ice
public wv_sat_svp_trans

! pressure -> humidity conversion
public wv_sat_svp_to_qsat

! Combined qsat operations
public wv_sat_qsat_water
public wv_sat_qsat_ice
public wv_sat_qsat_trans

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

subroutine wv_sat_methods_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wv_sat_methods_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('WV_SAT_METHODS_IMPL', value=impl_name, length=n, &
       status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wv_sat_methods_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wv_sat_methods_impl = .false.
  end if

  wv_sat_methods_impl_selected = .true.

  if (masterproc) then
     if (use_native_wv_sat_methods_impl) then
        write(iulog,*) 'wv_sat_methods implementation = native'
     else
        write(iulog,*) 'wv_sat_methods implementation = codon'
     end if
  end if

end subroutine wv_sat_methods_select_impl

subroutine wv_sat_methods_proof_once()

  if (wv_sat_methods_proof_written) return
  wv_sat_methods_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'wv_sat_methods entered (init constants = codon)'
  end if

end subroutine wv_sat_methods_proof_once

real(r8) function wv_sat_methods_value(value) result(out)
  real(r8), intent(in) :: value

  call wv_sat_methods_select_impl()

  if (use_native_wv_sat_methods_impl) then
     out = value
     return
  end if

  call wv_sat_methods_proof_once()
  out = real(wv_sat_methods_value_codon(real(value, c_double)), r8)

end function wv_sat_methods_value

real(r8) function wv_sat_methods_omeps(value) result(out)
  real(r8), intent(in) :: value

  call wv_sat_methods_select_impl()

  if (use_native_wv_sat_methods_impl) then
     out = 1._r8 - value
     return
  end if

  call wv_sat_methods_proof_once()
  out = real(wv_sat_methods_omeps_codon(real(value, c_double)), r8)

end function wv_sat_methods_omeps

! Get physical constants
subroutine wv_sat_methods_init(kind, tmelt_in, h2otrip_in, tboil_in, &
     ttrice_in, epsilo_in, errstring)
  integer, intent(in) :: kind
  real(r8), intent(in) :: tmelt_in
  real(r8), intent(in) :: h2otrip_in
  real(r8), intent(in) :: tboil_in
  real(r8), intent(in) :: ttrice_in
  real(r8), intent(in) :: epsilo_in
  character(len=*), intent(out)  :: errstring

  errstring = ' '

  if (kind /= r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          kind,' was input kind but ',r8,' is internal kind.'
     return
  end if

  if (ttrice_in < 0._r8) then
     write(errstring,*) 'wv_sat_methods_init: ERROR: ', &
          ttrice_in,' was input for ttrice, but negative range is invalid.'
     return
  end if

  call wv_sat_methods_select_impl()

  if (use_native_wv_sat_methods_impl) then
     tmelt = tmelt_in
     h2otrip = h2otrip_in
     tboil = tboil_in
     ttrice = ttrice_in
     epsilo = epsilo_in
     omeps = 1._r8 - epsilo
     return
  end if

  call wv_sat_methods_proof_once()
  call wv_sat_methods_init_codon(real(tmelt_in, c_double), real(h2otrip_in, c_double), &
       real(tboil_in, c_double), real(ttrice_in, c_double), real(epsilo_in, c_double), &
       c_loc(tmelt), c_loc(h2otrip), c_loc(tboil), c_loc(ttrice), c_loc(epsilo), c_loc(omeps))

end subroutine wv_sat_methods_init

! Look up index by name.
pure function wv_sat_get_scheme_idx(name) result(idx)
  character(len=*), intent(in) :: name
  integer :: idx
  integer :: i
  integer(c_int64_t), target :: name_ascii(len(name))

  do i = 1, len(name)
     name_ascii(i) = int(iachar(name(i:i)), c_int64_t)
  end do

  idx = int(wv_sat_get_scheme_idx_codon(int(len(name), c_int64_t), c_loc(name_ascii(1))))

end function wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
pure function wv_sat_valid_idx(idx) result(status)
  integer, intent(in) :: idx
  logical :: status

  status = wv_sat_valid_idx_codon(int(idx, c_int64_t)) /= 0_c_int64_t

end function wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
function wv_sat_set_default(name) result(status)
  character(len=*), intent(in) :: name
  logical :: status

  ! Don't want to overwrite valid default with invalid,
  ! so assign to temporary and check it first.
  integer :: tmp_idx

  tmp_idx = wv_sat_get_scheme_idx(name)

  status = wv_sat_set_default_codon(int(tmp_idx, c_int64_t), c_loc(default_idx)) /= 0_c_int64_t

end function wv_sat_set_default

! Reset default scheme to initial value.
! The same thing can be accomplished with wv_sat_set_default;
! the real reason to provide this routine is to reset the
! module for testing purposes.
subroutine wv_sat_reset_default()

  default_idx = initial_default_idx

end subroutine wv_sat_reset_default

pure function wv_sat_get_default_idx() result(idx)

  integer :: idx

  idx = default_idx

end function wv_sat_get_default_idx

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to range 0-1.
elemental function wv_sat_svp_to_qsat(es, p) result(qs)

  real(r8), intent(in) :: es  ! SVP
  real(r8), intent(in) :: p   ! Current pressure.
  real(r8) :: qs

  qs = real(wv_sat_svp_to_qsat_codon(real(es, c_double), real(p, c_double), &
       real(epsilo, c_double), real(omeps, c_double)), r8)

end function wv_sat_svp_to_qsat

elemental subroutine wv_sat_qsat_water(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), target, intent(out) :: es  ! Saturation vapor pressure
  real(r8), target, intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  call wv_sat_qsat_water_codon(real(t, c_double), real(p, c_double), &
       int(use_idx, c_int64_t), real(epsilo, c_double), real(omeps, c_double), &
       c_loc(es), c_loc(qs))

end subroutine wv_sat_qsat_water

elemental subroutine wv_sat_qsat_ice(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_ice(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_ice

elemental subroutine wv_sat_qsat_trans(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = wv_sat_svp_trans(t, idx)

  qs = wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine wv_sat_qsat_trans

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

elemental function wv_sat_svp_water(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  es = wv_sat_svp_water_codon(real(t, c_double), int(use_idx, c_int64_t))

end function wv_sat_svp_water

elemental function wv_sat_svp_ice(t, idx) result(es)
  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(r8) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  es = wv_sat_svp_ice_codon(real(t, c_double), int(use_idx, c_int64_t))

end function wv_sat_svp_ice

elemental function wv_sat_svp_trans(t, idx) result (es)

  real(r8), intent(in) :: t
  integer,  intent(in), optional :: idx

  real(r8) :: es
  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  es = wv_sat_svp_trans_codon(real(t, c_double), int(use_idx, c_int64_t), &
       real(tmelt, c_double), real(ttrice, c_double))

end function wv_sat_svp_trans

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Goff & Gratch (1946)

elemental function GoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = goffgratch_svp_water_codon(real(t, c_double))

end function GoffGratch_svp_water

elemental function GoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = goffgratch_svp_ice_codon(real(t, c_double))

end function GoffGratch_svp_ice

pure function goffgratch_svp_water_native_cb(t_c) result(es_c) &
     bind(C, name="goffgratch_svp_water_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t

  t = real(t_c, r8)

  ! uncertain below -70 C
  es_c = real(10._r8**(-7.90298_r8*(tboil/t-1._r8)+ &
       5.02808_r8*log10(tboil/t)- &
       1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/tboil))-1._r8)+ &
       8.1328e-3_r8*(10._r8**(-3.49149_r8*(tboil/t-1._r8))-1._r8)+ &
       log10(1013.246_r8))*100._r8, c_double)

end function goffgratch_svp_water_native_cb

pure function goffgratch_svp_ice_native_cb(t_c) result(es_c) &
     bind(C, name="goffgratch_svp_ice_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t

  t = real(t_c, r8)

  ! good down to -100 C
  es_c = real(10._r8**(-9.09718_r8*(h2otrip/t-1._r8)-3.56654_r8* &
       log10(h2otrip/t)+0.876793_r8*(1._r8-t/h2otrip)+ &
       log10(6.1071_r8))*100._r8, c_double)

end function goffgratch_svp_ice_native_cb

! Murphy & Koop (2005)

elemental function MurphyKoop_svp_water(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = murphykoop_svp_water_codon(real(t, c_double))

end function MurphyKoop_svp_water

elemental function MurphyKoop_svp_ice(t) result(es)
  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = murphykoop_svp_ice_codon(real(t, c_double))

end function MurphyKoop_svp_ice

pure function murphykoop_svp_water_native_cb(t_c) result(es_c) &
     bind(C, name="murphykoop_svp_water_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t

  t = real(t_c, r8)

  ! (good for 123 < T < 332 K)
  es_c = real(exp(54.842763_r8 - (6763.22_r8 / t) - (4.210_r8 * log(t)) + &
       (0.000367_r8 * t) + (tanh(0.0415_r8 * (t - 218.8_r8)) * &
       (53.878_r8 - (1331.22_r8 / t) - (9.44523_r8 * log(t)) + &
       0.014025_r8 * t))), c_double)

end function murphykoop_svp_water_native_cb

pure function murphykoop_svp_ice_native_cb(t_c) result(es_c) &
     bind(C, name="murphykoop_svp_ice_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t

  t = real(t_c, r8)

  ! (good down to 110 K)
  es_c = real(exp(9.550426_r8 - (5723.265_r8 / t) + (3.53068_r8 * log(t)) &
       - (0.00728332_r8 * t)), c_double)

end function murphykoop_svp_ice_native_cb

! Old CAM implementation, also labelled Goff & Gratch (1946)

! The water formula differs only due to compiler-dependent order of
! operations, so differences are roundoff level, usually 0.

! The ice formula gives fairly close answers to the current
! implementation, but has been rearranged, and uses the
! 1 atm melting point of water as the triple point.
! Differences are thus small but above roundoff.

! A curious fact: although using the melting point of water was
! probably a mistake, it mildly improves accuracy for ice svp,
! since it compensates for a systematic error in Goff & Gratch.

elemental function OldGoffGratch_svp_water(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es

  es = oldgoffgratch_svp_water_codon(real(t, c_double))

end function OldGoffGratch_svp_water

elemental function OldGoffGratch_svp_ice(t) result(es)
  real(r8), intent(in) :: t
  real(r8) :: es

  es = oldgoffgratch_svp_ice_codon(real(t, c_double))

end function OldGoffGratch_svp_ice

pure function oldgoffgratch_svp_water_native_cb(t_c) result(es_c) &
     bind(C, name="oldgoffgratch_svp_water_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t
  real(r8) :: ps, e1, e2, f1, f2, f3, f4, f5, f

  t = real(t_c, r8)

  ps = 1013.246_r8
  e1 = 11.344_r8*(1.0_r8 - t/tboil)
  e2 = -3.49149_r8*(tboil/t - 1.0_r8)
  f1 = -7.90298_r8*(tboil/t - 1.0_r8)
  f2 = 5.02808_r8*log10(tboil/t)
  f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
  f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
  f5 = log10(ps)
  f  = f1 + f2 + f3 + f4 + f5

  es_c = real((10.0_r8**f)*100.0_r8, c_double)

end function oldgoffgratch_svp_water_native_cb

pure function oldgoffgratch_svp_ice_native_cb(t_c) result(es_c) &
     bind(C, name="oldgoffgratch_svp_ice_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t
  real(r8) :: term1, term2, term3

  t = real(t_c, r8)

  term1 = 2.01889049_r8/(tmelt/t)
  term2 = 3.56654_r8*log(tmelt/t)
  term3 = 20.947031_r8*(tmelt/t)

  es_c = real(575.185606e10_r8*exp(-(term1 + term2 + term3)), c_double)

end function oldgoffgratch_svp_ice_native_cb

! Bolton (1980)
! zm_conv deep convection scheme contained this SVP calculation.
! It appears to be from D. Bolton, 1980, Monthly Weather Review.
! Unlike the other schemes, no distinct ice formula is associated
! with it. (However, a Bolton ice formula exists in CLUBB.)

! The original formula used degrees C, but this function
! takes Kelvin and internally converts.

elemental function Bolton_svp_water(t) result(es)
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  real(r8), intent(in) :: t  ! Temperature in Kelvin
  real(r8) :: es             ! SVP in Pa

  es = bolton_svp_water_codon(real(t, c_double))

end function Bolton_svp_water

pure function bolton_svp_water_native_cb(t_c) result(es_c) &
     bind(C, name="bolton_svp_water_native_cb")
  real(c_double), value, intent(in) :: t_c
  real(c_double) :: es_c
  real(r8) :: t
  real(r8),parameter :: c1 = 611.2_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8

  t = real(t_c, r8)

  es_c = real(c1*exp( (c2*(t - tmelt))/((t - tmelt)+c3) ), c_double)

end function bolton_svp_water_native_cb

end module wv_sat_methods

module pbl_utils
!-----------------------------------------------------------------------!
! Module to hold PBL-related subprograms that may be used with multiple !
! different vertical diffusion schemes.                                 !
!                                                                       !
! Public subroutines:                                                   !
!
!     calc_obklen                                                       !
!                                                                       !
!------------------ History --------------------------------------------!
! Created: Apr. 2012, by S. Santos                                      !
!-----------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8
use cam_logfile, only: iulog
use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
use spmd_utils, only: masterproc

implicit none
private

! Public Procedures
!----------------------------------------------------------------------!
! Excepting the initialization procedure, these are elemental
! procedures, so they can accept scalars or any dimension of array as
! arguments, as long as all arguments have the same number of
! elements.
public pbl_utils_init
public pbl_utils_log_pure_codon_counts
public calc_ustar
public calc_obklen
public virtem
public compute_radf

real(r8), parameter :: ustar_min = 0.01_r8

real(r8) :: g         ! acceleration of gravity
real(r8) :: vk        ! Von Karman's constant
real(r8) :: cpair     ! specific heat of dry air
real(r8) :: rair      ! gas constant for dry air
real(r8) :: zvir      ! rh2o/rair - 1
logical :: use_native_pbl_utils_impl = .false.
logical :: pbl_utils_impl_selected = .false.
logical :: pbl_utils_proof_written = .false.
logical :: pbl_utils_init_logged = .false.
logical :: calc_ustar_logged = .false.
logical :: calc_obklen_logged = .false.
logical :: virtem_logged = .false.
logical :: use_native_compute_radf_impl = .false.
logical :: compute_radf_impl_selected = .false.
logical :: compute_radf_logged = .false.

interface
  function pbl_utils_value_codon(value_c) result(value_out) &
       bind(c, name="pbl_utils_value_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: value_c
    real(c_double) :: value_out
  end function pbl_utils_value_codon
  function pbl_utils_init_codon(value_c) result(value_out) &
       bind(c, name="pbl_utils_init_codon")
    use iso_c_binding, only: c_double
    real(c_double), value :: value_c
    real(c_double) :: value_out
  end function pbl_utils_init_codon
  function physpkg_pure_counter_codon(which_c) result(count_c) &
       bind(c, name="physpkg_pure_counter_codon")
    use iso_c_binding, only: c_int64_t
    integer(c_int64_t), value :: which_c
    integer(c_int64_t) :: count_c
  end function physpkg_pure_counter_codon
  pure function calc_ustar_rrho_codon(rair_c, t_c, pmid_c) result(rrho_out) &
       bind(c, name="calc_ustar_rrho_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: rair_c, t_c, pmid_c
    real(c_double) :: rrho_out
  end function calc_ustar_rrho_codon
  pure function calc_ustar_codon(taux_c, tauy_c, rrho_c, ustar_min_c) result(ustar_out) &
       bind(c, name="calc_ustar_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: taux_c, tauy_c, rrho_c, ustar_min_c
    real(c_double) :: ustar_out
  end function calc_ustar_codon
  pure function calc_obklen_khfs_codon(shflx_c, rrho_c, cpair_c) result(khfs_out) &
       bind(c, name="calc_obklen_khfs_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: shflx_c, rrho_c, cpair_c
    real(c_double) :: khfs_out
  end function calc_obklen_khfs_codon
  pure function calc_obklen_kqfs_codon(qflx_c, rrho_c) result(kqfs_out) &
       bind(c, name="calc_obklen_kqfs_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: qflx_c, rrho_c
    real(c_double) :: kqfs_out
  end function calc_obklen_kqfs_codon
  pure function calc_obklen_kbfs_codon(khfs_c, zvir_c, ths_c, kqfs_c) result(kbfs_out) &
       bind(c, name="calc_obklen_kbfs_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: khfs_c, zvir_c, ths_c, kqfs_c
    real(c_double) :: kbfs_out
  end function calc_obklen_kbfs_codon
  pure function calc_obklen_codon(thvs_c, ustar_c, g_c, vk_c, kbfs_c) result(obklen_out) &
       bind(c, name="calc_obklen_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: thvs_c, ustar_c, g_c, vk_c, kbfs_c
    real(c_double) :: obklen_out
  end function calc_obklen_codon
  pure function virtem_codon(t_c, q_c, zvir_c) result(value_out) &
       bind(c, name="virtem_codon")
    use iso_c_binding, only: c_double
    real(c_double), intent(in), value :: t_c, q_c, zvir_c
    real(c_double) :: value_out
  end function virtem_codon
  subroutine compute_radf_codon(i_c, pcols_c, pver_c, ncvmax_c, radf_mode_c, qmin_c, g_c, &
       ncvfin_p, ktop_p, ql_p, pi_p, qrlw_p, cldeff_p, zi_p, chs_p, lwp_CL_p, opt_depth_CL_p, &
       radinvfrac_CL_p, radf_CL_p) bind(c, name="compute_radf_codon")
    use iso_c_binding, only: c_double, c_int64_t, c_ptr
    integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvmax_c, radf_mode_c
    real(c_double), value :: qmin_c, g_c
    type(c_ptr), value :: ncvfin_p, ktop_p, ql_p, pi_p, qrlw_p, cldeff_p, zi_p, chs_p
    type(c_ptr), value :: lwp_CL_p, opt_depth_CL_p, radinvfrac_CL_p, radf_CL_p
  end subroutine compute_radf_codon
end interface

contains

subroutine pbl_utils_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (pbl_utils_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('PBL_UTILS_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_pbl_utils_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_pbl_utils_impl = .false.
  end if

  pbl_utils_impl_selected = .true.

  if (masterproc) then
     if (use_native_pbl_utils_impl) then
        write(iulog,*) 'pbl_utils implementation = native'
     else
        write(iulog,*) 'pbl_utils implementation = codon'
     end if
  end if

end subroutine pbl_utils_select_impl

subroutine pbl_utils_proof_once()

  if (pbl_utils_proof_written) return
  pbl_utils_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'pbl_utils entered (init constants = codon)'
  end if

end subroutine pbl_utils_proof_once

subroutine pbl_utils_log_direct(logged, proof_line)

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
  end if

end subroutine pbl_utils_log_direct

subroutine pbl_utils_log_pure_codon_counts()

  integer(c_int64_t) :: hits

  call pbl_utils_select_impl()
  if (use_native_pbl_utils_impl) return

  hits = physpkg_pure_counter_codon(1_c_int64_t)
  if (hits > 0_c_int64_t) then
     call pbl_utils_log_direct(calc_ustar_logged, 'calc_ustar direct = codon; pure counter proof')
  end if

  hits = physpkg_pure_counter_codon(2_c_int64_t)
  if (hits > 0_c_int64_t) then
     call pbl_utils_log_direct(calc_obklen_logged, 'calc_obklen direct = codon; pure counter proof')
  end if

  hits = physpkg_pure_counter_codon(3_c_int64_t)
  if (hits > 0_c_int64_t) then
     call pbl_utils_log_direct(virtem_logged, 'virtem direct = codon; pure counter proof')
  end if

end subroutine pbl_utils_log_pure_codon_counts

subroutine pbl_utils_compute_radf_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (compute_radf_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('PBL_UTILS_COMPUTE_RADF_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_compute_radf_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_compute_radf_impl = .false.
  end if

  compute_radf_impl_selected = .true.

  if (masterproc) then
     if (use_native_compute_radf_impl) then
        write(iulog,*) 'pbl_utils_compute_radf implementation = native'
     else
        write(iulog,*) 'pbl_utils_compute_radf implementation = codon'
     end if
  end if

end subroutine pbl_utils_compute_radf_select_impl

subroutine pbl_utils_compute_radf_log_direct()

  if (compute_radf_logged) return
  compute_radf_logged = .true.

  if (masterproc) then
     write(iulog,'(A)') 'compute_radf direct = codon'
  end if

end subroutine pbl_utils_compute_radf_log_direct

real(r8) function pbl_utils_value(value) result(out)
  real(r8), intent(in) :: value

  call pbl_utils_select_impl()

  if (use_native_pbl_utils_impl) then
     out = value
     return
  end if

  call pbl_utils_proof_once()
  out = real(pbl_utils_value_codon(real(value, c_double)), r8)

end function pbl_utils_value

subroutine pbl_utils_init(g_in,vk_in,cpair_in,rair_in,zvir_in)

  !-----------------------------------------------------------------------!
  ! Purpose: Set constants to be used in calls to later functions         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: g_in       ! acceleration of gravity
  real(r8), intent(in) :: vk_in      ! Von Karman's constant
  real(r8), intent(in) :: cpair_in   ! specific heat of dry air
  real(r8), intent(in) :: rair_in    ! gas constant for dry air
  real(r8), intent(in) :: zvir_in    ! rh2o/rair - 1

  call pbl_utils_select_impl()

  if (use_native_pbl_utils_impl) then
     g = g_in
     vk = vk_in
     cpair = cpair_in
     rair = rair_in
     zvir = zvir_in
     return
  end if

  call pbl_utils_proof_once()
  g = real(pbl_utils_init_codon(real(g_in, c_double)), r8)
  vk = real(pbl_utils_init_codon(real(vk_in, c_double)), r8)
  cpair = real(pbl_utils_init_codon(real(cpair_in, c_double)), r8)
  rair = real(pbl_utils_init_codon(real(rair_in, c_double)), r8)
  zvir = real(pbl_utils_init_codon(real(zvir_in, c_double)), r8)

  if (.not. use_native_pbl_utils_impl) then
     call pbl_utils_log_direct(pbl_utils_init_logged, 'pbl_utils_init direct = codon')
  end if

end subroutine pbl_utils_init

elemental subroutine calc_ustar( t,    pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t         ! surface temperature
  real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

  real(r8), intent(out) :: rrho     ! 1./bottom level density
  real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

  if (use_native_pbl_utils_impl) then
     rrho = rair * t / pmid
     ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
     return
  end if

  rrho = real(calc_ustar_rrho_codon(real(rair, c_double), real(t, c_double), &
       real(pmid, c_double)), r8)
  ustar = real(calc_ustar_codon(real(taux, c_double), real(tauy, c_double), &
       real(rrho, c_double), real(ustar_min, c_double)), r8)

end subroutine calc_ustar

elemental subroutine calc_obklen( ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!

  real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]

  real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen        ! Obukhov length

  if (use_native_pbl_utils_impl) then
     khfs = shflx*rrho/cpair
     kqfs = qflx*rrho
     kbfs = khfs + zvir*ths*kqfs
     obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))
     return
  end if

  khfs = real(calc_obklen_khfs_codon(real(shflx, c_double), real(rrho, c_double), &
       real(cpair, c_double)), r8)
  kqfs = real(calc_obklen_kqfs_codon(real(qflx, c_double), real(rrho, c_double)), r8)
  kbfs = real(calc_obklen_kbfs_codon(real(khfs, c_double), real(zvir, c_double), &
       real(ths, c_double), real(kqfs, c_double)), r8)
  obklen = real(calc_obklen_codon(real(thvs, c_double), real(ustar, c_double), &
       real(g, c_double), real(vk, c_double), real(kbfs, c_double)), r8)

end subroutine calc_obklen

elemental real(r8) function virtem(t,q)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t, q

  if (use_native_pbl_utils_impl) then
     virtem = t * (1.0_r8 + zvir*q)
     return
  end if

  virtem = real(virtem_codon(real(t, c_double), real(q, c_double), real(zvir, c_double)), r8)

end function virtem

subroutine compute_radf( choice_radf, i, pcols, pver, ncvmax, ncvfin, ktop, qmin, &
                         ql, pi, qrlw, g, cldeff, zi, chs, lwp_CL, opt_depth_CL,  &
                         radinvfrac_CL, radf_CL )
  ! -------------------------------------------------------------------------- !
  ! Purpose:                                                                   !
  ! Calculate cloud-top radiative cooling contribution to buoyancy production. !
  ! Here,  'radf' [m2/s3] is additional buoyancy flux at the CL top interface  !
  ! associated with cloud-top LW cooling being mainly concentrated near the CL !
  ! top interface ( just below CL top interface ).  Contribution of SW heating !
  ! within the cloud is not included in this radiative buoyancy production     !
  ! since SW heating is more broadly distributed throughout the CL top layer.  !
  ! -------------------------------------------------------------------------- !

  !-----------------!
  ! Input variables !
  !-----------------!
  character(len=6), intent(in) :: choice_radf  ! Method for calculating radf
  integer,  intent(in)  :: i                   ! Index of current column
  integer,  intent(in)  :: pcols               ! Number of atmospheric columns
  integer,  intent(in)  :: pver                ! Number of atmospheric layers
  integer,  intent(in)  :: ncvmax              ! Max numbers of CLs (perhaps equal to pver)
  integer,  target, intent(in)  :: ncvfin(pcols)       ! Total number of CL in column
  integer,  target, intent(in)  :: ktop(pcols, ncvmax) ! ktop for current column
  real(r8), intent(in)  :: qmin                ! Minimum grid-mean LWC counted as clouds [kg/kg]
  real(r8), target, intent(in)  :: ql(pcols, pver)     ! Liquid water specific humidity [ kg/kg ]
  real(r8), target, intent(in)  :: pi(pcols, pver+1)   ! Interface pressures [ Pa ]
  real(r8), target, intent(in)  :: qrlw(pcols, pver)   ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
  real(r8), intent(in)  :: g                   ! Gravitational acceleration
  real(r8), target, intent(in)  :: cldeff(pcols,pver)  ! Effective Cloud Fraction [fraction]
  real(r8), target, intent(in)  :: zi(pcols, pver+1)   ! Interface heights [ m ]
  real(r8), target, intent(in)  :: chs(pcols, pver+1)  ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces.

  !------------------!
  ! Output variables !
  !------------------!
  real(r8), target, intent(out) :: lwp_CL(ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
  real(r8), target, intent(out) :: opt_depth_CL(ncvmax)   ! Optical depth of the CL top layer
  real(r8), target, intent(out) :: radinvfrac_CL(ncvmax)  ! Fraction of LW radiative cooling confined in the top portion of CL
  real(r8), target, intent(out) :: radf_CL(ncvmax)        ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]

  !-----------------!
  ! Local variables !
  !-----------------!
  integer :: kt, ncv
  integer :: radf_mode
  real(r8) :: lwp, opt_depth, radinvfrac, radf


  !-----------------!
  ! Begin main code !
  !-----------------!
  call pbl_utils_compute_radf_select_impl()

  if (.not. use_native_compute_radf_impl) then
    radf_mode = 2
    if( choice_radf .eq. 'orig' ) then
      radf_mode = 0
    elseif( choice_radf .eq. 'ramp' ) then
      radf_mode = 1
    endif

    call pbl_utils_compute_radf_log_direct()
    call compute_radf_codon(int(i, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(ncvmax, c_int64_t), int(radf_mode, c_int64_t), real(qmin, c_double), real(g, c_double), &
         c_loc(ncvfin), c_loc(ktop), c_loc(ql), c_loc(pi), c_loc(qrlw), c_loc(cldeff), c_loc(zi), c_loc(chs), &
         c_loc(lwp_CL), c_loc(opt_depth_CL), c_loc(radinvfrac_CL), c_loc(radf_CL))
    return
  end if

  lwp_CL        = 0._r8
  opt_depth_CL  = 0._r8
  radinvfrac_CL = 0._r8
  radf_CL       = 0._r8

  ! ---------------------------------------- !
  ! Perform do loop for individual CL regime !
  ! ---------------------------------------- !
  do ncv = 1, ncvfin(i)
    kt = ktop(i,ncv)
    !-----------------------------------------------------!
    ! Compute radf for each CL regime and for each column !
    !-----------------------------------------------------!
    if( choice_radf .eq. 'orig' ) then
      if( ql(i,kt) .gt. qmin .and. ql(i,kt-1) .lt. qmin ) then
        lwp       = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        ! Approximate LW cooling fraction concentrated at the inversion by using
        ! polynomial approx to exact formula 1-2/opt_depth+2/(exp(opt_depth)-1))

        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling = [ W/kg ]
        radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)
        ! We can disable cloud LW cooling contribution to turbulence by uncommenting:
        ! radf = 0._r8
      end if

    elseif( choice_radf .eq. 'ramp' ) then

      lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
      opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
      radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
      radinvfrac  = max(cldeff(i,kt)-cldeff(i,kt-1),0._r8) * radinvfrac
      radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling [W/kg]
      radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)

    elseif( choice_radf .eq. 'maxi' ) then

      ! Radiative flux divergence both in 'kt' and 'kt-1' layers are included
      ! 1. From 'kt' layer
        lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = max( radinvfrac * qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 )
      ! 2. From 'kt-1' layer and add the contribution from 'kt' layer
        lwp         = ql(i,kt-1) * ( pi(i,kt) - pi(i,kt-1) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth) + opt_depth**2 )
        radf        = radf + max( radinvfrac * qrlw(i,kt-1) / ( pi(i,kt-1) - pi(i,kt) ) * ( zi(i,kt-1) - zi(i,kt) ), 0._r8 )
        radf        = max( radf, 0._r8 ) * chs(i,kt)

    endif

    lwp_CL(ncv)        = lwp
    opt_depth_CL(ncv)  = opt_depth
    radinvfrac_CL(ncv) = radinvfrac
    radf_CL(ncv)       = radf
  end do ! ncv = 1, ncvfin(i)
end subroutine compute_radf

end module pbl_utils

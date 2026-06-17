module gw_oro

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use gw_utils, only: r8
use coords_1d, only: Coords1D
use spmd_utils, only: masterproc
use cam_logfile, only: iulog

implicit none
private
save

! Public interface
public :: gw_oro_src

logical :: use_native_gw_oro_src_impl = .false.
logical :: gw_oro_src_impl_selected = .false.
logical :: gw_oro_src_entered_logged = .false.

contains

!==========================================================================

subroutine gw_oro_src(ncol, band, p, &
     u, v, t, sgh, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_common, only: GWBand, pver, rair
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  !-----------------------------------------------------------------------
  ! Selectable wrapper for the active fixed-case orographic source helper.
  ! The Codon path receives only Fortran-owned arrays/workspaces.
  !-----------------------------------------------------------------------
  integer, intent(in) :: ncol
  type(GWBand), intent(in) :: band
  type(Coords1D), target, intent(in) :: p
  real(r8), target, intent(in) :: u(ncol,pver), v(ncol,pver)
  real(r8), target, intent(in) :: t(ncol,pver)
  real(r8), target, intent(in) :: sgh(ncol)
  real(r8), target, intent(in) :: zm(ncol,pver)
  real(r8), target, intent(in) :: nm(ncol,pver)
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)
  real(r8), target, intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  real(r8), target, intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  real(r8), target, intent(out) :: xv(ncol), yv(ncol)
  real(r8), target, intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  integer :: i
  integer(c_int64_t), target :: src_level64(ncol), tend_level64(ncol)
  real(r8), target :: hdsp(ncol), tauoro(ncol), nsrc(ncol), rsrc(ncol)
  real(r8), target :: usrc(ncol), vsrc(ncol), dpsrc(ncol)

  interface
     subroutine gw_oro_src_codon(ncol_c, pver_c, ngwv_c, fcrit2_c, kwv_c, rair_c, &
          p_mid_p, p_del_p, p_ifc_p, u_p, v_p, t_p, sgh_p, zm_p, nm_p, &
          src_level_p, tend_level_p, tau_p, ubm_p, ubi_p, xv_p, yv_p, c_p, &
          hdsp_p, tauoro_p, nsrc_p, rsrc_p, usrc_p, vsrc_p, dpsrc_p) &
          bind(c, name="gw_oro_src_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, ngwv_c
       real(c_double), value :: fcrit2_c, kwv_c, rair_c
       type(c_ptr), value :: p_mid_p, p_del_p, p_ifc_p, u_p, v_p, t_p, sgh_p, zm_p, nm_p
       type(c_ptr), value :: src_level_p, tend_level_p, tau_p, ubm_p, ubi_p, xv_p, yv_p, c_p
       type(c_ptr), value :: hdsp_p, tauoro_p, nsrc_p, rsrc_p, usrc_p, vsrc_p, dpsrc_p
     end subroutine gw_oro_src_codon
  end interface

  call gw_oro_src_select_impl()

  if (use_native_gw_oro_src_impl) then
     call gw_oro_src_native(ncol, band, p, u, v, t, sgh, zm, nm, &
          src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  else
     call gw_oro_src_note_entered()
     call gw_oro_src_codon(int(ncol, c_int64_t), int(pver, c_int64_t), int(band%ngwv, c_int64_t), &
          real(band%fcrit2, c_double), real(band%kwv, c_double), real(rair, c_double), &
          c_loc(p%mid), c_loc(p%del), c_loc(p%ifc), c_loc(u), c_loc(v), c_loc(t), c_loc(sgh), &
          c_loc(zm), c_loc(nm), c_loc(src_level64), c_loc(tend_level64), c_loc(tau), c_loc(ubm), &
          c_loc(ubi), c_loc(xv), c_loc(yv), c_loc(c), c_loc(hdsp), c_loc(tauoro), c_loc(nsrc), &
          c_loc(rsrc), c_loc(usrc), c_loc(vsrc), c_loc(dpsrc))
     do i = 1, ncol
        src_level(i) = int(src_level64(i))
        tend_level(i) = int(tend_level64(i))
     end do
  end if

end subroutine gw_oro_src

!==========================================================================

subroutine gw_oro_src_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_ORO_SRC_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_oro_src_append_proof

!==========================================================================

subroutine gw_oro_src_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_oro_src_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('GW_ORO_SRC_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_oro_src_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_oro_src_impl = .false.
  end if

  gw_oro_src_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_oro_src_impl) then
        write(iulog,*) 'gw_oro_src implementation = native'
        call gw_oro_src_append_proof('gw_oro_src selector entered implementation = native')
     else
        write(iulog,*) 'gw_oro_src implementation = codon'
        call gw_oro_src_append_proof('gw_oro_src selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_oro_src_select_impl

!==========================================================================

subroutine gw_oro_src_note_entered()

  if (gw_oro_src_entered_logged) return
  gw_oro_src_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_oro_src entered (unified orographic-source stage dispatch = codon)'
     call gw_oro_src_append_proof('gw_oro_src entered (unified orographic-source stage dispatch = codon)')
     call flush(iulog)
  end if

end subroutine gw_oro_src_note_entered

!==========================================================================

subroutine gw_oro_src_codon_wrap(ncol_local, pver_local, ngwv_local, fcrit2_local, kwv_local, rair_local, &
     p_mid_local, p_del_local, p_ifc_local, u_local, v_local, t_local, sgh_local, zm_local, nm_local, &
     src_level64_local, tend_level64_local, tau_local, ubm_local, ubi_local, xv_local, yv_local, c_local, &
     hdsp_local, tauoro_local, nsrc_local, rsrc_local, usrc_local, vsrc_local, dpsrc_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local, pver_local, ngwv_local
  real(r8), intent(in) :: fcrit2_local, kwv_local, rair_local
  real(r8), target, intent(in) :: p_mid_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: p_del_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: p_ifc_local(ncol_local,pver_local+1)
  real(r8), target, intent(in) :: u_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: v_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: t_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: sgh_local(ncol_local)
  real(r8), target, intent(in) :: zm_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: nm_local(ncol_local,pver_local)
  integer(c_int64_t), target, intent(inout) :: src_level64_local(ncol_local)
  integer(c_int64_t), target, intent(inout) :: tend_level64_local(ncol_local)
  real(r8), target, intent(inout) :: tau_local(ncol_local,-ngwv_local:ngwv_local,pver_local+1)
  real(r8), target, intent(inout) :: ubm_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: ubi_local(ncol_local,pver_local+1)
  real(r8), target, intent(inout) :: xv_local(ncol_local)
  real(r8), target, intent(inout) :: yv_local(ncol_local)
  real(r8), target, intent(inout) :: c_local(ncol_local,-ngwv_local:ngwv_local)
  real(r8), target, intent(inout) :: hdsp_local(ncol_local), tauoro_local(ncol_local)
  real(r8), target, intent(inout) :: nsrc_local(ncol_local), rsrc_local(ncol_local)
  real(r8), target, intent(inout) :: usrc_local(ncol_local), vsrc_local(ncol_local), dpsrc_local(ncol_local)

  interface
     subroutine gw_oro_src_stage_dispatch_codon(ncol_c, pver_c, ngwv_c, fcrit2_c, kwv_c, rair_c, &
          p_mid_p, p_del_p, p_ifc_p, u_p, v_p, t_p, sgh_p, zm_p, nm_p, &
          src_level_p, tend_level_p, tau_p, ubm_p, ubi_p, xv_p, yv_p, c_p, &
          hdsp_p, tauoro_p, nsrc_p, rsrc_p, usrc_p, vsrc_p, dpsrc_p) &
          bind(c, name="gw_oro_src_stage_dispatch_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, ngwv_c
       real(c_double), value :: fcrit2_c, kwv_c, rair_c
       type(c_ptr), value :: p_mid_p, p_del_p, p_ifc_p, u_p, v_p, t_p, sgh_p, zm_p, nm_p
       type(c_ptr), value :: src_level_p, tend_level_p, tau_p, ubm_p, ubi_p, xv_p, yv_p, c_p
       type(c_ptr), value :: hdsp_p, tauoro_p, nsrc_p, rsrc_p, usrc_p, vsrc_p, dpsrc_p
     end subroutine gw_oro_src_stage_dispatch_codon
  end interface

  call gw_oro_src_stage_dispatch_codon(int(ncol_local, c_int64_t), int(pver_local, c_int64_t), int(ngwv_local, c_int64_t), &
       real(fcrit2_local, c_double), real(kwv_local, c_double), real(rair_local, c_double), &
       c_loc(p_mid_local), c_loc(p_del_local), c_loc(p_ifc_local), c_loc(u_local), c_loc(v_local), c_loc(t_local), &
       c_loc(sgh_local), c_loc(zm_local), c_loc(nm_local), c_loc(src_level64_local), c_loc(tend_level64_local), &
       c_loc(tau_local), c_loc(ubm_local), c_loc(ubi_local), c_loc(xv_local), c_loc(yv_local), c_loc(c_local), &
       c_loc(hdsp_local), c_loc(tauoro_local), c_loc(nsrc_local), c_loc(rsrc_local), c_loc(usrc_local), &
       c_loc(vsrc_local), c_loc(dpsrc_local))

end subroutine gw_oro_src_codon_wrap

!==========================================================================

subroutine gw_oro_src_native(ncol, band, p, &
     u, v, t, sgh, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_common, only: GWBand, pver, rair
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(ncol)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)
  ! Midpoint Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real(r8) :: hdsp(ncol)
  ! Max orographic standard deviation to use.
  real(r8) :: sghmax
  ! c=0 stress from orography.
  real(r8) :: tauoro(ncol)
  ! Averages over source region.
  real(r8) :: nsrc(ncol) ! B-V frequency.
  real(r8) :: rsrc(ncol) ! Density.
  real(r8) :: usrc(ncol) ! Zonal wind.
  real(r8) :: vsrc(ncol) ! Meridional wind.

  ! Difference in interface pressure across source region.
  real(r8) :: dpsrc(ncol)

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves
  real(r8), parameter :: orohmin = 10._r8
  ! min wind speed for orographic waves
  real(r8), parameter :: orovmin = 2._r8

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp = 2.0_r8 * sgh

  k = pver
  src_level = k-1
  rsrc = p%mid(:,k)/(rair*t(:,k)) * p%del(:,k)
  usrc = u(:,k) * p%del(:,k)
  vsrc = v(:,k) * p%del(:,k)
  nsrc = nm(:,k)* p%del(:,k)

  do k = pver-1, 1, -1
     do i = 1, ncol
        if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
           src_level(i) = k-1
           rsrc(i) = rsrc(i) + &
                p%mid(i,k) / (rair*t(i,k)) * p%del(i,k)
           usrc(i) = usrc(i) + u(i,k) * p%del(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * p%del(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* p%del(i,k)
        end if
     end do
     ! Break the loop when all source levels found.
     if (all(src_level >= k)) exit
  end do

  do i = 1, ncol
     dpsrc(i) = p%ifc(i,pver+1) - p%ifc(i,src_level(i)+1)
  end do

  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,pver+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  do i = 1, ncol
     if ((ubi(i,pver+1) > orovmin) .and. (hdsp(i) > orohmin)) then
        sghmax = band%fcrit2 * (ubi(i,pver+1) / nsrc(i))**2
        tauoro(i) = 0.5_r8 * band%kwv * min(hdsp(i)**2, sghmax) * &
             rsrc(i) * nsrc(i) * ubi(i,pver+1)
     else
        tauoro(i) = 0._r8
        src_level(i) = pver
     end if
  end do

  ! Set the phase speeds and wave numbers in the direction of the source
  ! wind. Set the source stress magnitude (positive only, note that the
  ! sign of the stress is the same as (c-u).
  tau = 0._r8
  do k = pver, minval(src_level), -1
     where (src_level <= k) tau(:,0,k+1) = tauoro
  end do

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0._r8

end subroutine gw_oro_src_native

end module gw_oro

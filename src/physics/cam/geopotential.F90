
module geopotential

!---------------------------------------------------------------------------------
! Compute geopotential from temperature or
! compute geopotential and temperature from dry static energy.
!
! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
! The diagonal element is the itegration weight from interface k+1 to midpoint k.
! The offdiagonal element is the weight between interfaces.
! 
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pverp
  use dycore,       only: dycore_is
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc

  implicit none
  private
  save

  public geopotential_dse
  public geopotential_t

  logical :: use_native_geopotential_impl = .false.
  logical :: geopotential_impl_selected = .false.
  logical :: geopotential_proof_written = .false.

  interface
     subroutine geopotential_dse_codon(ncol_c, ld_c, pver_c, pverp_c, fvdyn_c, gravit_c, &
          piln_p, pint_p, pmid_p, pdel_p, rpdel_p, dse_p, q_p, phis_p, rair_p, cpair_p, &
          zvir_p, t_p, zi_p, zm_p) bind(c, name="geopotential_dse_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, ld_c, pver_c, pverp_c, fvdyn_c
       real(c_double), value :: gravit_c
       type(c_ptr), value :: piln_p, pint_p, pmid_p, pdel_p, rpdel_p, dse_p, q_p, phis_p
       type(c_ptr), value :: rair_p, cpair_p, zvir_p, t_p, zi_p, zm_p
     end subroutine geopotential_dse_codon

     subroutine geopotential_t_codon(ncol_c, ld_c, pver_c, pverp_c, fvdyn_c, gravit_c, &
          piln_p, pint_p, pmid_p, pdel_p, rpdel_p, t_p, q_p, rair_p, zvir_p, zi_p, zm_p) &
          bind(c, name="geopotential_t_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, ld_c, pver_c, pverp_c, fvdyn_c
       real(c_double), value :: gravit_c
       type(c_ptr), value :: piln_p, pint_p, pmid_p, pdel_p, rpdel_p, t_p, q_p, rair_p
       type(c_ptr), value :: zvir_p, zi_p, zm_p
     end subroutine geopotential_t_codon
  end interface

contains
!===============================================================================
  subroutine geopotential_select_impl()
    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (geopotential_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('GEOPOTENTIAL_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_geopotential_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_geopotential_impl = .false.
    end if

    geopotential_impl_selected = .true.

    if (masterproc) then
       if (use_native_geopotential_impl) then
          write(iulog,*) 'geopotential implementation = native'
       else
          write(iulog,*) 'geopotential implementation = codon'
       end if
    end if
  end subroutine geopotential_select_impl

!===============================================================================
  subroutine geopotential_proof_once()
    if (geopotential_proof_written) return
    geopotential_proof_written = .true.
    if (masterproc) then
       write(iulog,'(A)') 'geopotential helpers entered (dse/t hydrostatic loops = codon)'
    end if
  end subroutine geopotential_proof_once

!===============================================================================
  subroutine geopotential_dse(                                &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol             )
!-----------------------------------------------------------------------
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: piln (:,:)
    real(r8), target, intent(in) :: pmln (:,:)
    real(r8), target, intent(in) :: pint (:,:)
    real(r8), target, intent(in) :: pmid (:,:)
    real(r8), target, intent(in) :: pdel (:,:)
    real(r8), target, intent(in) :: rpdel(:,:)
    real(r8), target, intent(in) :: dse  (:,:)
    real(r8), target, intent(in) :: q    (:,:)
    real(r8), target, intent(in) :: phis (:)
    real(r8), target, intent(in) :: rair (:,:)
    real(r8), intent(in) :: gravit
    real(r8), target, intent(in) :: cpair(:,:)
    real(r8), target, intent(in) :: zvir (:,:)
    real(r8), target, intent(out) :: t(:,:)
    real(r8), target, intent(out) :: zi(:,:)
    real(r8), target, intent(out) :: zm(:,:)

    call geopotential_select_impl()
    if (use_native_geopotential_impl) then
       call geopotential_dse_native(piln, pmln, pint, pmid, pdel, rpdel, dse, q, phis, &
            rair, gravit, cpair, zvir, t, zi, zm, ncol)
    else
       call geopotential_proof_once()
       call geopotential_dse_codon(int(ncol, c_int64_t), int(size(piln, 1), c_int64_t), &
            int(pver, c_int64_t), int(pverp, c_int64_t), &
            merge(1_c_int64_t, 0_c_int64_t, dycore_is('LR')), real(gravit, c_double), &
            c_loc(piln(1,1)), c_loc(pint(1,1)), c_loc(pmid(1,1)), c_loc(pdel(1,1)), &
            c_loc(rpdel(1,1)), c_loc(dse(1,1)), c_loc(q(1,1)), c_loc(phis(1)), &
            c_loc(rair(1,1)), c_loc(cpair(1,1)), c_loc(zvir(1,1)), c_loc(t(1,1)), &
            c_loc(zi(1,1)), c_loc(zm(1,1)))
    end if

    return
  end subroutine geopotential_dse

!===============================================================================
  subroutine geopotential_dse_native(                         &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol             )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the temperature  and geopotential height (above the surface) at the
! midpoints and interfaces from the input dry static energy and pressures.
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
! Input arguments
    integer, intent(in) :: ncol                  ! Number of longitudes

    ! rair, and cpair are passed in as slices of rank 3 arrays allocated
    ! at runtime. Don't specify size to avoid temporary copy.
    real(r8), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(r8), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(r8), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(r8), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(r8), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(r8), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(r8), intent(in) :: dse  (:,:)    ! (pcols,pver)  - dry static energy
    real(r8), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(r8), intent(in) :: phis (:)      ! (pcols)       - surface geopotential
    real(r8), intent(in) :: rair (:,:)    !               - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: cpair(:,:)    !               - specific heat at constant p for dry air
    real(r8), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: t(:,:)       ! (pcols,pver)  - temperature
    real(r8), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(r8), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------------------
!
    logical  :: fvdyn                   ! finite volume dynamics
    integer  :: i,k                     ! Lon, level, level indices
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
!
!----------------------------------------------------------------------------------
    rog(:ncol,:) = rair(:ncol,:) / gravit

! Set dynamics flag
    fvdyn = dycore_is ('LR')

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

! Compute the virtual temperature, zi, zm from bottom up
! Note, zi(i,k) is the interface above zm(i,k)
    do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics
       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

! Now compute tv, t, zm, zi
       do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k)
          tv      = (dse(i,k) - phis(i) - gravit*zi(i,k+1)) / ((cpair(i,k) / tvfac) + &
	                                                               rair(i,k)*hkk(i))

          t (i,k) = tv / tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_dse_native

!===============================================================================
  subroutine geopotential_t(                                 &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   )
!-----------------------------------------------------------------------
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: piln (:,:)
    real(r8), target, intent(in) :: pmln (:,:)
    real(r8), target, intent(in) :: pint (:,:)
    real(r8), target, intent(in) :: pmid (:,:)
    real(r8), target, intent(in) :: pdel (:,:)
    real(r8), target, intent(in) :: rpdel(:,:)
    real(r8), target, intent(in) :: t    (:,:)
    real(r8), target, intent(in) :: q    (:,:)
    real(r8), target, intent(in) :: rair (:,:)
    real(r8), intent(in) :: gravit
    real(r8), target, intent(in) :: zvir (:,:)
    real(r8), target, intent(out) :: zi(:,:)
    real(r8), target, intent(out) :: zm(:,:)

    call geopotential_select_impl()
    if (use_native_geopotential_impl) then
       call geopotential_t_native(piln, pmln, pint, pmid, pdel, rpdel, t, q, rair, gravit, &
            zvir, zi, zm, ncol)
    else
       call geopotential_proof_once()
       call geopotential_t_codon(int(ncol, c_int64_t), int(size(piln, 1), c_int64_t), &
            int(pver, c_int64_t), int(pverp, c_int64_t), &
            merge(1_c_int64_t, 0_c_int64_t, dycore_is('LR')), real(gravit, c_double), &
            c_loc(piln(1,1)), c_loc(pint(1,1)), c_loc(pmid(1,1)), c_loc(pdel(1,1)), &
            c_loc(rpdel(1,1)), c_loc(t(1,1)), c_loc(q(1,1)), c_loc(rair(1,1)), &
            c_loc(zvir(1,1)), c_loc(zi(1,1)), c_loc(zm(1,1)))
    end if

    return
  end subroutine geopotential_t

!===============================================================================
  subroutine geopotential_t_native(                          &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------

use ppgrid, only : pcols

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol                  ! Number of longitudes

    real(r8), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(r8), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(r8), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(r8), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(r8), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(r8), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(r8), intent(in) :: t    (:,:)    ! (pcols,pver)  - temperature
    real(r8), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(r8), intent(in) :: rair (:,:)    ! (pcols,pver)  - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(r8), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    logical  :: fvdyn                   ! finite volume dynamics
    integer  :: i,k                     ! Lon, level indices
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
!
!-----------------------------------------------------------------------
!
    rog(:ncol,:) = rair(:ncol,:) / gravit

! Set dynamics flag

    fvdyn = dycore_is ('LR')

! The surface height is zero by definition.

    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)

    do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics

       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

! Now compute tv, zm, zi

       do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t_native
end module geopotential

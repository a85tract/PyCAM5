module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
! May 2014  F Vitt -- included contributions from MAM4 aerosols and added soa_a2 to the ocphiwet fluxes
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t
use constituents,     only: pcnst, cnst_get_ind
use ppgrid,           only: pcols
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog
use spmd_utils,       only: masterproc
use iso_c_binding,    only: c_double, c_int64_t, c_loc, c_ptr

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep,             &
   set_srf_wetdep_codon_direct

! Private module data
integer :: idx_bc1  = -1
integer :: idx_pom1 = -1
integer :: idx_soa1 = -1
integer :: idx_soa2 = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43 = -1
integer :: idx_bc4  = -1
integer :: idx_pom4 = -1

logical :: bin_fluxes = .false.

logical :: initialized = .false.
logical :: use_native_set_srf_wetdep_impl = .false.
logical :: set_srf_wetdep_impl_selected = .false.
logical :: use_native_set_srf_drydep_impl = .false.
logical :: set_srf_drydep_impl_selected = .false.

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init(bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx, &
                            dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx,bc4_ndx,pom4_ndx)

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! can be initialized with user specified indices
   ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified

   integer, optional, intent(in) :: bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx,dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx
   integer, optional, intent(in) :: bc4_ndx,pom4_ndx

   ! if already initialized abort the run
   if (initialized) then
     call endrun('modal_aero_deposition_init is already initialized')
   endif

   if (present(bc1_ndx)) then
      idx_bc1  = bc1_ndx
   else
      call cnst_get_ind('bc_a1',  idx_bc1)
   endif
   if (present(pom1_ndx)) then
      idx_pom1 = pom1_ndx
   else
      call cnst_get_ind('pom_a1', idx_pom1)
   endif
   if (present(soa1_ndx)) then
      idx_soa1 = soa1_ndx
   else
      call cnst_get_ind('soa_a1', idx_soa1)
   endif
   if (present(soa2_ndx)) then
      idx_soa2 = soa2_ndx
   else
      call cnst_get_ind('soa_a2', idx_soa2)
   endif
   if (present(dst1_ndx)) then
      idx_dst1 = dst1_ndx
   else
      call cnst_get_ind('dst_a1', idx_dst1,abort=.false.)
   endif
   if (present(dst3_ndx)) then
      idx_dst3 = dst3_ndx
   else
      call cnst_get_ind('dst_a3', idx_dst3,abort=.false.)
   endif
   if (present(ncl3_ndx)) then
      idx_ncl3 = ncl3_ndx
   else
      call cnst_get_ind('ncl_a3', idx_ncl3,abort=.false.)
   endif
   if (present(so43_ndx)) then
      idx_so43 = so43_ndx
   else
      call cnst_get_ind('so4_a3', idx_so43,abort=.false.)
   endif
   if (present(bc4_ndx)) then
      idx_bc4 = bc4_ndx
   else
      call cnst_get_ind('bc_a4', idx_bc4,abort=.false.)
   endif
   if (present(pom4_ndx)) then
      idx_pom4 = pom4_ndx
   else
      call cnst_get_ind('pom_a4', idx_pom4,abort=.false.)   
   endif

!  for 7 mode bin_fluxes will be false
   bin_fluxes = idx_dst1>0 .and. idx_dst3>0 .and.idx_ncl3>0 .and. idx_so43>0
   initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (set_srf_wetdep_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('SET_SRF_WETDEP_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_set_srf_wetdep_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_set_srf_wetdep_impl = .false.
   end if

   set_srf_wetdep_impl_selected = .true.

   if (masterproc) then
      if (use_native_set_srf_wetdep_impl) then
         write(iulog,*) 'set_srf_wetdep implementation = native'
      else
         write(iulog,*) 'set_srf_wetdep implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine set_srf_wetdep_select_impl

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), target, intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), target, intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), target, intent(inout) :: cam_out ! cam export state

   call set_srf_wetdep_select_impl()

   if (.not.bin_fluxes) return

   if (use_native_set_srf_wetdep_impl) then
      call set_srf_wetdep_native(aerdepwetis, aerdepwetcw, cam_out)
      return
   end if

   call set_srf_wetdep_codon_invoke(aerdepwetis, aerdepwetcw, cam_out)

end subroutine set_srf_wetdep

!==============================================================================
subroutine set_srf_wetdep_codon_direct(aerdepwetis, aerdepwetcw, cam_out)

! Direct Codon entry used by active shell callers that already selected Codon.

   real(r8), target, intent(in) :: aerdepwetis(:,:)
   real(r8), target, intent(in) :: aerdepwetcw(:,:)
   type(cam_out_t), target, intent(inout) :: cam_out

   call set_srf_wetdep_codon_invoke(aerdepwetis, aerdepwetcw, cam_out)

end subroutine set_srf_wetdep_codon_direct

!==============================================================================
subroutine set_srf_wetdep_codon_invoke(aerdepwetis, aerdepwetcw, cam_out)

   real(r8), target, intent(in) :: aerdepwetis(:,:)
   real(r8), target, intent(in) :: aerdepwetcw(:,:)
   type(cam_out_t), target, intent(inout) :: cam_out

   integer :: ncol

   interface
      subroutine set_srf_wetdep_codon(ncol_c, pcols_c, idx_bc1_c, idx_bc4_c, idx_pom1_c, idx_pom4_c, idx_soa1_c, &
           idx_soa2_c, idx_dst1_c, idx_dst3_c, aerdepwetis_p, aerdepwetcw_p, bcphiwet_p, ocphiwet_p, dstwet1_p, &
           dstwet2_p, dstwet3_p, dstwet4_p) bind(c, name="set_srf_wetdep_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, idx_bc1_c, idx_bc4_c, idx_pom1_c, idx_pom4_c, idx_soa1_c
         integer(c_int64_t), value :: idx_soa2_c, idx_dst1_c, idx_dst3_c
         type(c_ptr), value :: aerdepwetis_p, aerdepwetcw_p, bcphiwet_p, ocphiwet_p, dstwet1_p, dstwet2_p
         type(c_ptr), value :: dstwet3_p, dstwet4_p
      end subroutine set_srf_wetdep_codon
   end interface

   if (.not.bin_fluxes) return

   ncol = cam_out%ncol

   call set_srf_wetdep_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(idx_bc1, c_int64_t), int(idx_bc4, c_int64_t), &
        int(idx_pom1, c_int64_t), int(idx_pom4, c_int64_t), int(idx_soa1, c_int64_t), int(idx_soa2, c_int64_t), &
        int(idx_dst1, c_int64_t), int(idx_dst3, c_int64_t), c_loc(aerdepwetis(1,1)), c_loc(aerdepwetcw(1,1)), &
        c_loc(cam_out%bcphiwet(1)), c_loc(cam_out%ocphiwet(1)), c_loc(cam_out%dstwet1(1)), c_loc(cam_out%dstwet2(1)), &
        c_loc(cam_out%dstwet3(1)), c_loc(cam_out%dstwet4(1)) &
   )

end subroutine set_srf_wetdep_codon_invoke

!==============================================================================
subroutine set_srf_wetdep_native(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   ncol = cam_out%ncol

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      if (idx_bc1>0) &
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc1)+aerdepwetcw(i,idx_bc1))
      if (idx_bc4>0) &
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc4)+aerdepwetcw(i,idx_bc4))

      ! organic carbon fluxes
      if (idx_soa1>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_soa1)+aerdepwetcw(i,idx_soa1))
      if (idx_soa2>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_soa2)+aerdepwetcw(i,idx_soa2))
      if (idx_pom1>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom1)+aerdepwetcw(i,idx_pom1))
      if (idx_pom4>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom4)+aerdepwetcw(i,idx_pom4))

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstwet1(i) = -(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))
      
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = -(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3))
      cam_out%dstwet4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
   enddo

end subroutine set_srf_wetdep_native

!==============================================================================

subroutine set_srf_drydep_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (set_srf_drydep_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('SET_SRF_DRYDEP_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_set_srf_drydep_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_set_srf_drydep_impl = .false.
   end if

   set_srf_drydep_impl_selected = .true.

   if (masterproc) then
      if (use_native_set_srf_drydep_impl) then
         write(iulog,*) 'set_srf_drydep implementation = native'
      else
         write(iulog,*) 'set_srf_drydep implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine set_srf_drydep_select_impl

!==============================================================================

subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

   real(r8), target, intent(in)    :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), target, intent(in)    :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), target, intent(inout) :: cam_out    ! cam export state

   integer :: ncol

   interface
      subroutine set_srf_drydep_codon(ncol_c, pcols_c, idx_bc1_c, idx_bc4_c, idx_pom1_c, idx_pom4_c, idx_soa1_c, &
           idx_soa2_c, idx_dst1_c, idx_dst3_c, aerdepdryis_p, aerdepdrycw_p, bcphidry_p, bcphodry_p, ocphidry_p, &
           ocphodry_p, dstdry1_p, dstdry2_p, dstdry3_p, dstdry4_p) bind(c, name="set_srf_drydep_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, idx_bc1_c, idx_bc4_c, idx_pom1_c, idx_pom4_c, idx_soa1_c
         integer(c_int64_t), value :: idx_soa2_c, idx_dst1_c, idx_dst3_c
         type(c_ptr), value :: aerdepdryis_p, aerdepdrycw_p, bcphidry_p, bcphodry_p, ocphidry_p, ocphodry_p
         type(c_ptr), value :: dstdry1_p, dstdry2_p, dstdry3_p, dstdry4_p
      end subroutine set_srf_drydep_codon
   end interface

   call set_srf_drydep_select_impl()

   if (.not.bin_fluxes) return

   if (use_native_set_srf_drydep_impl) then
      call set_srf_drydep_native(aerdepdryis, aerdepdrycw, cam_out)
      return
   end if

   ncol = cam_out%ncol

   call set_srf_drydep_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(idx_bc1, c_int64_t), int(idx_bc4, c_int64_t), &
        int(idx_pom1, c_int64_t), int(idx_pom4, c_int64_t), int(idx_soa1, c_int64_t), int(idx_soa2, c_int64_t), &
        int(idx_dst1, c_int64_t), int(idx_dst3, c_int64_t), c_loc(aerdepdryis(1,1)), c_loc(aerdepdrycw(1,1)), &
        c_loc(cam_out%bcphidry(1)), c_loc(cam_out%bcphodry(1)), c_loc(cam_out%ocphidry(1)), c_loc(cam_out%ocphodry(1)), &
        c_loc(cam_out%dstdry1(1)), c_loc(cam_out%dstdry2(1)), c_loc(cam_out%dstdry3(1)), c_loc(cam_out%dstdry4(1)) &
   )

end subroutine set_srf_drydep

!==============================================================================

subroutine set_srf_drydep_native(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------
   ncol = cam_out%ncol

   cam_out%bcphidry(:) = 0._r8
   cam_out%bcphodry(:) = 0._r8
   cam_out%ocphidry(:) = 0._r8
   cam_out%ocphodry(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      if (idx_bc1>0) &
           cam_out%bcphidry(i) = cam_out%bcphidry(i) + aerdepdryis(i,idx_bc1)+aerdepdrycw(i,idx_bc1)
      if (idx_bc4>0) &
           cam_out%bcphodry(i) = cam_out%bcphodry(i) + aerdepdryis(i,idx_bc4)+aerdepdrycw(i,idx_bc4)

      ! organic carbon fluxes
      if (idx_pom1>0) &
           cam_out%ocphidry(i) = cam_out%ocphidry(i) + aerdepdryis(i,idx_pom1)+aerdepdrycw(i,idx_pom1)
      if (idx_pom4>0) &
           cam_out%ocphodry(i) = cam_out%ocphodry(i) + aerdepdryis(i,idx_pom4)+aerdepdrycw(i,idx_pom4)
      if (idx_soa1>0) &
           cam_out%ocphidry(i) = cam_out%ocphidry(i) + aerdepdryis(i,idx_soa1)+aerdepdrycw(i,idx_soa1)
      if (idx_soa2>0) &
           cam_out%ocphodry(i) = cam_out%ocphodry(i) + aerdepdryis(i,idx_soa2)+aerdepdrycw(i,idx_soa2)

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstdry1(i) = aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1)
      
      ! Two options for partitioning deposition into bins 2-4:
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)
      cam_out%dstdry4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
   enddo

end subroutine set_srf_drydep_native


!==============================================================================

end module modal_aero_deposition

module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM driver layer for aerosol activation processes.
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                  microphysics and macrophysics and concentrate aerosol information here
!                B. Eaton, Sep 2014 - Refactored to move CAM interface code into the CAM
!                  interface modules and preserve just the driver layer functionality here.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc

use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl, hetfrz_classnuc_cam_register, hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_save_cbaero, hetfrz_classnuc_cam_calc

use cam_history,      only: addfld, phys_decomp, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,       only: endrun
use iso_c_binding,    only: c_double, c_int64_t

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

! Private module data

character(len=16)   :: eddy_scheme
logical             :: micro_do_icesupersat

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8) :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: cldo_idx = -1
integer :: dgnumwet_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4

! modal aerosols
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: npccn_idx, rndst_idx, nacon_idx

logical  :: separate_dust = .false.
logical  :: use_native_run_impl = .false.
logical  :: run_impl_selected = .false.
logical  :: run_linear_entered_logged = .false.
logical  :: run_npccn_copy_entered_logged = .false.

interface
   function microp_aero_register_codon(flag_c) result(out_c) bind(c, name="microp_aero_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_aero_register_codon
   function microp_aero_init_codon(flag_c) result(out_c) bind(c, name="microp_aero_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_aero_init_codon
   function microp_aero_readnl_codon(value_c) result(out_c) bind(c, name="microp_aero_readnl_codon")
      use iso_c_binding, only: c_double
      real(c_double), value :: value_c
      real(c_double) :: out_c
   end function microp_aero_readnl_codon
end interface

!=========================================================================================
contains
!=========================================================================================

subroutine microp_aero_register
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Register pbuf fields for aerosols needed by microphysics
   ! 
   ! Author: Cheryl Craig October 2012
   ! 
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   if (microp_aero_register_codon(1_c_int64_t) == 0_c_int64_t) return

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)

   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)
 
   call nucleate_ice_cam_register()
   call hetfrz_classnuc_cam_register()

end subroutine microp_aero_register

!=========================================================================================

subroutine microp_aero_init

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Initialize constants for aerosols needed by microphysics
   ! 
   ! Author: Andrew Gettelman May 2010
   ! 
   !-----------------------------------------------------------------------

   ! local variables
   integer  :: iaer, ierr
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   integer(c_int64_t) :: active_c
   !-----------------------------------------------------------------------
   active_c = microp_aero_init_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme,  &
                     history_amwg_out         = history_amwg, &
                     micro_do_icesupersat_out = micro_do_icesupersat)

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')   
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2_nadv')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx      = pbuf_get_index('AST')

   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      call ndrop_init()

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

   else

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer) )

         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
      end do

      call ndrop_bam_init()

   end if

   call addfld('LCLOUD', ' ', pver, 'A', 'Liquid cloud fraction used in stratus activation', phys_decomp)

   call addfld('WSUB     ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity'                   ,phys_decomp)
   call addfld('WSUBI    ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity for ice'           ,phys_decomp)

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

   call nucleate_ice_cam_init(mincld, bulk_scale)
   call hetfrz_classnuc_cam_init(mincld)

end subroutine microp_aero_init

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variable
   call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
#endif

   ! set local variables
   bulk_scale = microp_aero_readnl_codon(real(microp_aero_bulk_scale, c_double))

   call nucleate_ice_cam_readnl(nlfile)
   call hetfrz_classnuc_cam_readnl(nlfile)

end subroutine microp_aero_readnl

!=========================================================================================

subroutine microp_aero_run_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (run_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MICROP_AERO_RUN_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_run_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_run_impl = .false.
   end if

   run_impl_selected = .true.

   if (masterproc) then
      if (use_native_run_impl) then
         write(iulog,*) 'microp_aero_run implementation = native'
      else
         write(iulog,*) 'microp_aero_run implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine microp_aero_run_select_impl

!=========================================================================================

subroutine microp_aero_run_linear_log_entered()

   if (run_linear_entered_logged) return
   run_linear_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'microp_aero_linear entered (rho/diag_TKE/lcloud/modal contact helpers = codon)'
      call flush(iulog)
   end if

end subroutine microp_aero_run_linear_log_entered

!=========================================================================================

subroutine microp_aero_run_npccn_copy_log_entered()

   if (run_npccn_copy_entered_logged) return
   run_npccn_copy_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'microp_aero_npccn_copy entered (dropmixnuc NPCCN transfer = codon)'
      call flush(iulog)
   end if

end subroutine microp_aero_run_npccn_copy_log_entered

!=========================================================================================

subroutine microp_aero_run ( &
   state, ptend, deltatin, pbuf)

   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   call microp_aero_run_select_impl()
   call microp_aero_run_codon(state, ptend, deltatin, pbuf, use_native_run_impl)

end subroutine microp_aero_run

!=========================================================================================

subroutine microp_aero_run_codon ( &
   state, ptend, deltatin, pbuf, use_native_impl)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   logical,                     intent(in)    :: use_native_impl

   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m
   integer :: itim_old
   integer :: nmodes

   real(r8), pointer :: ast(:,:)        

   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)

   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8), target :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8), target :: lcldm(pcols,pver)   ! liq cloud fraction

   real(r8), target :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8), target :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: qcld                ! total cloud water
   real(r8), target :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc           ! variables for modal scheme.

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8), target :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8), target :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)
   real(r8) :: nucboast

   real(r8) :: wght

   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number

   interface
      subroutine microp_aero_init_fields_codon(ncol_c, pcols_c, pver_c, rn_dst1_c, rn_dst2_c, rn_dst3_c, rn_dst4_c, &
           npccn_p, nacon_p, rndst_p) &
           bind(c, name="microp_aero_init_fields_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: rn_dst1_c, rn_dst2_c, rn_dst3_c, rn_dst4_c
         type(c_ptr), value :: npccn_p, nacon_p, rndst_p
      end subroutine microp_aero_init_fields_codon
      subroutine microp_aero_rho_codon(ncol_c, pcols_c, pver_c, top_lev_c, rair_c, pmid_p, t_p, rho_p) &
           bind(c, name="microp_aero_rho_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: rair_c
         type(c_ptr), value :: pmid_p, t_p, rho_p
      end subroutine microp_aero_rho_codon
      subroutine microp_aero_diag_tke_wsub_codon(ncol_c, pcols_c, pver_c, top_lev_c, use_preexisting_ice_c, &
           tke_p, wsub_p, wsubi_p) bind(c, name="microp_aero_diag_tke_wsub_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, use_preexisting_ice_c
         type(c_ptr), value :: tke_p, wsub_p, wsubi_p
      end subroutine microp_aero_diag_tke_wsub_codon
      subroutine microp_aero_lcldm_codon(ncol_c, pcols_c, pver_c, top_lev_c, mincld_c, ast_p, lcldm_p) &
           bind(c, name="microp_aero_lcldm_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: mincld_c
         type(c_ptr), value :: ast_p, lcldm_p
      end subroutine microp_aero_lcldm_codon
      subroutine microp_aero_modal_lcloud_codon(ncol_c, pcols_c, pver_c, top_lev_c, qsmall_c, qc_p, qi_p, cldn_p, &
           cldo_p, lcldn_p, lcldo_p) bind(c, name="microp_aero_modal_lcloud_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: qsmall_c
         type(c_ptr), value :: qc_p, qi_p, cldn_p, cldo_p, lcldn_p, lcldo_p
      end subroutine microp_aero_modal_lcloud_codon
      subroutine microp_aero_modal_contact_codon(ncol_c, pcols_c, pver_c, top_lev_c, nmodes_c, mode_coarse_dst_idx_c, &
           separate_dust_c, rn_dst3_c, t_p, rho_p, coarse_dust_p, coarse_nacl_p, num_coarse_p, dgnumwet_p, nacon_p, &
           rndst_p) bind(c, name="microp_aero_modal_contact_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, nmodes_c, mode_coarse_dst_idx_c
         integer(c_int64_t), value :: separate_dust_c
         real(c_double), value :: rn_dst3_c
         type(c_ptr), value :: t_p, rho_p, coarse_dust_p, coarse_nacl_p, num_coarse_p, dgnumwet_p, nacon_p, rndst_p
      end subroutine microp_aero_modal_contact_codon
      subroutine microp_aero_npccn_copy_codon(ncol_c, pcols_c, pver_c, nctend_mixnuc_p, npccn_p) &
           bind(c, name="microp_aero_npccn_copy_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: nctend_mixnuc_p, npccn_p
      end subroutine microp_aero_npccn_copy_codon
   end interface
   !-------------------------------------------------------------------------------

   associate( &
      lchnk => state%lchnk,             &
      ncol  => state%ncol,              &
      t     => state%t,                 &
      qc    => state%q(:,:,cldliq_idx), &
      qi    => state%q(:,:,cldice_idx), &
      nc    => state%q(:,:,numliq_idx), &
      pmid  => state%pmid               )

   itim_old = pbuf_old_tim_idx()
   if (micro_do_icesupersat) then
     call pbuf_get_field(pbuf, cldo_idx,     ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))   
   else
     call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   endif

   call pbuf_get_field(pbuf, npccn_idx, npccn)

   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   if (clim_modal_aero) then

      itim_old = pbuf_old_tim_idx()
      
      if (micro_do_icesupersat) then
        call pbuf_get_field(pbuf, cldo_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))        
      else
        call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      endif

      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )


      call rad_cnst_get_info(0, nmodes=nmodes)
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

      allocate(factnum(pcols,pver,nmodes))

   end if

   ! initialize output and time-varying parameters
   if (use_native_impl) then
      npccn(1:ncol,1:pver)    = 0._r8

      nacon(1:ncol,1:pver,:)  = 0._r8

      ! set default or fixed dust bins for contact freezing
      rndst(1:ncol,1:pver,1) = rn_dst1
      rndst(1:ncol,1:pver,2) = rn_dst2
      rndst(1:ncol,1:pver,3) = rn_dst3
      rndst(1:ncol,1:pver,4) = rn_dst4
   else
      call microp_aero_init_fields_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(rn_dst1, c_double), &
           real(rn_dst2, c_double), real(rn_dst3, c_double), real(rn_dst4, c_double), c_loc(npccn(1,1)), &
           c_loc(nacon(1,1,1)), c_loc(rndst(1,1,1)) &
      )
   end if

   ! save copy of cloud borne aerosols for use in heterogeneous freezing
   if (use_hetfrz_classnuc) then
      call hetfrz_classnuc_cam_save_cbaero(state, pbuf)
   end if

   if (use_native_impl) then
      do k = top_lev, pver
         do i = 1, ncol
            rho(i,k) = pmid(i,k)/(rair*t(i,k))
         end do
      end do
   else
      call microp_aero_run_linear_log_entered()
      call microp_aero_rho_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           real(rair, c_double), c_loc(pmid(1,1)), c_loc(t(1,1)), c_loc(rho(1,1)) &
      )
   end if

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state, pbuf, num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)

   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
         naer2(pcols,pver,naer_all), &
         maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)
         
         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity 
   ! Set to be zero at the surface by initialization.

   select case (trim(eddy_scheme))
   case ('diag_TKE')
      call pbuf_get_field(pbuf, tke_idx, tke)
   case ('CLUBB_SGS')
      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))
      allocate(tke(pcols,pverp))
      tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
   end select

   if (.not. use_native_impl .and. trim(eddy_scheme) == 'diag_TKE') then
      call microp_aero_diag_tke_wsub_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           merge(1_c_int64_t, 0_c_int64_t, use_preexisting_ice), c_loc(tke(1,1)), c_loc(wsub(1,1)), &
           c_loc(wsubi(1,1)) &
      )
   else
      ! Set minimum values above top_lev.
      wsub(:ncol,:top_lev-1)  = 0.20_r8
      wsubi(:ncol,:top_lev-1) = 0.001_r8

      do k = top_lev, pver
         do i = 1, ncol

            select case (trim(eddy_scheme))
            case ('diag_TKE', 'CLUBB_SGS')
               wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
               wsub(i,k) = min(wsub(i,k),10._r8)
            case default
               ! get sub-grid vertical velocity from diff coef.
               ! following morrison et al. 2005, JAS
               ! assume mixing length of 30 m
               dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
               ! use maximum sub-grid vertical vel of 10 m/s
               dum = min(dum, 10._r8)
               ! set wsub to value at current vertical level
               wsub(i,k)  = dum
            end select

            if (eddy_scheme == 'CLUBB_SGS') then
               wsubi(i,k) = max(0.2_r8, wsub(i,k))
               wsubi(i,k) = min(wsubi(i,k), 10.0_r8)
            else
               wsubi(i,k) = max(0.001_r8, wsub(i,k))
               if (.not. use_preexisting_ice) then
                  wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
               endif
            endif

            wsub(i,k)  = max(0.20_r8, wsub(i,k))

         end do
      end do
   end if

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI', wsubi, pcols, lchnk)

   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   call nucleate_ice_cam_calc(state, wsubi, pbuf)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! get liquid cloud fraction, check for minimum

   if (use_native_impl) then
      do k = top_lev, pver
         do i = 1, ncol
            lcldm(i,k) = max(ast(i,k), mincld)
         end do
      end do
   else
      call microp_aero_lcldm_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           real(mincld, c_double), c_loc(ast(1,1)), c_loc(lcldm(1,1)) &
      )
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation

   if (clim_modal_aero) then

      ! for modal aerosol

      ! partition cloud fraction into liquid water part
      if (use_native_impl) then
         lcldn = 0._r8
         lcldo = 0._r8
         do k = top_lev, pver
            do i = 1, ncol
               qcld = qc(i,k) + qi(i,k)
               if (qcld > qsmall) then
                  lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
                  lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
               end if
            end do
         end do
      else
         call microp_aero_modal_lcloud_codon( &
              int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
              real(qsmall, c_double), c_loc(qc(1,1)), c_loc(qi(1,1)), c_loc(cldn(1,1)), c_loc(cldo(1,1)), &
              c_loc(lcldn(1,1)), c_loc(lcldo(1,1)) &
         )
      end if

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      call dropmixnuc( &
         state, ptend, deltatin, pbuf, wsub, &
         lcldn, lcldo, nctend_mixnuc, factnum)

      if (use_native_impl) then
         npccn(:ncol,:) = nctend_mixnuc(:ncol,:)
      else
         call microp_aero_run_npccn_copy_log_entered()
         call microp_aero_npccn_copy_codon( &
              int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              c_loc(nctend_mixnuc(1,1)), c_loc(npccn(1,1)) &
         )
      end if

   else

      ! for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend, state%psetcols, 'none')

      do k = top_lev, pver
         do i = 1, ncol

            if (qc(i,k) >= qsmall) then

               ! get droplet activation rate

               call ndrop_bam_run( &
                  wsub(i,k), t(i,k), rho(i,k), naer2(i,k,:), naer_all, &
                  naer_all, maerosol(i,k,:),  &
                  dum2)
               dum = dum2
            else
               dum = 0._r8
            end if

            npccn(i,k) = (dum*lcldm(i,k) - nc(i,k))/deltatin
         end do
      end do

   end if


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics

   if (.not. use_native_impl .and. clim_modal_aero) then
      call microp_aero_modal_contact_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           int(nmodes, c_int64_t), int(mode_coarse_dst_idx, c_int64_t), &
           merge(1_c_int64_t, 0_c_int64_t, separate_dust), real(rn_dst3, c_double), c_loc(t(1,1)), &
           c_loc(rho(1,1)), c_loc(coarse_dust(1,1)), c_loc(coarse_nacl(1,1)), c_loc(num_coarse(1,1)), &
           c_loc(dgnumwet(1,1,1)), c_loc(nacon(1,1,1)), c_loc(rndst(1,1,1)) &
      )
   else
      do k = top_lev, pver
         do i = 1, ncol

            if (t(i,k) < 269.15_r8) then

               if (clim_modal_aero) then

                  ! For modal aerosols:
                  !  use size '3' for dust coarse mode...
                  !  scale by dust fraction in coarse mode
               
                  dmc  = coarse_dust(i,k)
                  ssmc = coarse_nacl(i,k)

                  if ( separate_dust ) then
                     ! 7-mode -- has separate dust and seasalt mode types and no need for weighting
                     wght = 1._r8
                  else
                     ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                     wght = dmc/(ssmc + dmc)
                  endif

                  if (dmc > 0.0_r8) then
                     nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
                  else
                     nacon(i,k,3) = 0._r8
                  end if

                  !also redefine parameters based on size...

                  rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
                  if (rndst(i,k,3) <= 0._r8) then
                     rndst(i,k,3) = rn_dst3
                  end if

               else

                  !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

                  if (idxdst2 > 0) then
                     nacon(i,k,2) = naer2(i,k,idxdst2)
                  end if
                  if (idxdst3 > 0) then
                     nacon(i,k,3) = naer2(i,k,idxdst3)
                  end if
                  if (idxdst4 > 0) then
                     nacon(i,k,4) = naer2(i,k,idxdst4)
                  end if
               end if
            end if
         end do
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)

   if (.not. clim_modal_aero) then

      ! ccn concentration as diagnostic
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

      deallocate( &
         naer2,    &
         maerosol)

   end if

   ! heterogeneous freezing
   if (use_hetfrz_classnuc) then

      call hetfrz_classnuc_cam_calc(state, deltatin, factnum, pbuf)

   end if

   if (clim_modal_aero) then
      deallocate(factnum)
   end if

   end associate

end subroutine microp_aero_run_codon

!=========================================================================================

end module microp_aero

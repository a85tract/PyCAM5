module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: pi, rair, tmelt
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, phys_decomp, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water, svp_water, svp_ice
use shr_spfn_mod,   only: erf => shr_spfn_erf

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use iso_c_binding, only: c_double, c_int64_t

use nucleate_ice,   only: nucleati_init, nucleati


implicit none
private
save

public :: &
   nucleate_ice_cam_readnl,   &
   nucleate_ice_cam_register, &
   nucleate_ice_cam_init,     &
   nucleate_ice_cam_calc
   

! Namelist variables
logical, public, protected :: use_preexisting_ice = .false.
logical                    :: hist_preexisting_ice = .false.
real(r8)                   :: nucleate_ice_subgrid

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction
real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numice_idx = -1

integer :: &
   naai_idx,     &
   naai_hom_idx

integer :: &
   ast_idx   = -1, &
   dgnum_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: clim_modal_aero

integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

logical  :: separate_dust = .false.
real(r8) :: sigmag_aitken

logical :: nucleate_ice_cam_prep_use_native_impl = .false.
logical :: nucleate_ice_cam_prep_impl_selected = .false.
logical :: nucleate_ice_cam_prep_entered_logged = .false.
logical :: nucleate_ice_cam_post_entered_logged = .false.
logical :: nucleate_ice_cam_modal_dust_entered_logged = .false.
logical :: nucleate_ice_cam_modal_so4_entered_logged = .false.
logical :: nucleate_ice_cam_modal_nucleati_entered_logged = .false.
logical :: nucleate_ice_cam_readnl_logged = .false.
logical :: nucleate_ice_cam_register_logged = .false.
logical :: nucleate_ice_cam_init_logged = .false.

interface
   function nucleate_ice_cam_readnl_codon(flag_c) result(out_c) bind(c, name="nucleate_ice_cam_readnl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function nucleate_ice_cam_readnl_codon
   function nucleate_ice_cam_register_codon(flag_c) result(out_c) bind(c, name="nucleate_ice_cam_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function nucleate_ice_cam_register_codon
   function nucleate_ice_cam_init_codon(flag_c) result(out_c) bind(c, name="nucleate_ice_cam_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function nucleate_ice_cam_init_codon
   function nucleate_ice_cam_init_mincld_codon(value_c) result(out_c) bind(c, name="nucleate_ice_cam_init_mincld_codon")
      use iso_c_binding, only: c_double
      real(c_double), value :: value_c
      real(c_double) :: out_c
   end function nucleate_ice_cam_init_mincld_codon
   function nucleate_ice_cam_init_bulk_scale_codon(value_c) result(out_c) &
        bind(c, name="nucleate_ice_cam_init_bulk_scale_codon")
      use iso_c_binding, only: c_double
      real(c_double), value :: value_c
      real(c_double) :: out_c
   end function nucleate_ice_cam_init_bulk_scale_codon
end interface

!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'
  integer(c_int64_t) :: active_c

  namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
       nucleate_ice_subgrid

  !-----------------------------------------------------------------------------
  active_c = nucleate_ice_cam_readnl_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return
  call nucleate_ice_cam_log_direct(nucleate_ice_cam_readnl_logged, &
       'nucleate_ice_cam_readnl direct = codon; namelist/MPI native islands')

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, nucleate_ice_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(use_preexisting_ice,  1, mpilog, 0, mpicom)
  call mpibcast(hist_preexisting_ice, 1, mpilog, 0, mpicom)
  call mpibcast(nucleate_ice_subgrid, 1, mpir8, 0, mpicom)
#endif

end subroutine nucleate_ice_cam_readnl

!================================================================================================

subroutine nucleate_ice_cam_register()

   if (nucleate_ice_cam_register_codon(1_c_int64_t) == 0_c_int64_t) return
   call nucleate_ice_cam_log_direct(nucleate_ice_cam_register_logged, &
        'nucleate_ice_cam_register direct = codon; pbuf_add_field native CAM API island')

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

!================================================================================================

subroutine nucleate_ice_cam_init(mincld_in, bulk_scale_in)

   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: bulk_scale_in

   ! local variables
   integer  :: iaer
   integer  :: m, n, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
   integer(c_int64_t) :: init_touch_c
   !--------------------------------------------------------------------------------------------

   init_touch_c = nucleate_ice_cam_init_codon(1_c_int64_t)
   if (init_touch_c /= 0_c_int64_t) then
      call nucleate_ice_cam_log_direct(nucleate_ice_cam_init_logged, &
           'nucleate_ice_cam_init direct = codon; scalar init direct = codon; registration/nucleati native boundaries')
   end if

   mincld     = nucleate_ice_cam_init_mincld_codon(real(mincld_in, c_double))
   bulk_scale = nucleate_ice_cam_init_bulk_scale_codon(real(bulk_scale_in, c_double))

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   call addfld('NIHF',  '1/m3', pver, 'A', 'Activated Ice Number Concentation due to homogenous freezing',  phys_decomp)
   call addfld('NIDEP', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to deposition nucleation',phys_decomp)
   call addfld('NIIMM', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to immersion freezing',   phys_decomp)
   call addfld('NIMEY', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to meyers deposition',    phys_decomp)

   if (use_preexisting_ice) then
      call addfld('fhom     ', 'fraction', pver, 'A', 'Fraction of cirrus where homogeneous freezing occur'   ,phys_decomp) 
      call addfld ('WICE      ', 'm/s   ', pver, 'A','Vertical velocity Reduction caused by preexisting ice'  ,phys_decomp)
      call addfld ('WEFF      ', 'm/s   ', pver, 'A','Effective Vertical velocity for ice nucleation' ,phys_decomp)
      call addfld ('INnso4    ','1/m3   ', pver, 'A','Number Concentation so4 used for ice_nucleation',phys_decomp)
      call addfld ('INnbc     ','1/m3   ', pver, 'A','Number Concentation bc  used for ice_nucleation',phys_decomp)
      call addfld ('INndust   ','1/m3   ', pver, 'A','Number Concentation dustused for ice_nucleation',phys_decomp)
      call addfld ('INhet     ','1/m3   ', pver, 'A', &
                'contribution for in-cloud ice number density increase by het nucleation in ice cloud',phys_decomp)
      call addfld ('INhom     ','1/m3   ', pver, 'A', &
                'contribution for in-cloud ice number density increase by hom nucleation in ice cloud',phys_decomp)
      call addfld ('INFrehom  ','frequency',pver,'A','hom IN frequency ice cloud',phys_decomp)
      call addfld ('INFreIN   ','frequency',pver,'A','frequency of ice nucleation occur',phys_decomp)

      if (hist_preexisting_ice) then
         call add_default ('WSUBI   ', 1, ' ')  ! addfld/outfld calls are in microp_aero

         call add_default ('fhom    ', 1, ' ') 
         call add_default ('WICE    ', 1, ' ')
         call add_default ('WEFF    ', 1, ' ')
         call add_default ('INnso4  ', 1, ' ')
         call add_default ('INnbc   ', 1, ' ')
         call add_default ('INndust ', 1, ' ')
         call add_default ('INhet   ', 1, ' ')
         call add_default ('INhom   ', 1, ' ')
         call add_default ('INFrehom', 1, ' ')
         call add_default ('INFreIN ', 1, ' ')
      end if
   end if

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   if (clim_modal_aero) then

      dgnum_idx    = pbuf_get_index('DGNUM' )

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
      if (mode_coarse_dst_idx < 0) mode_coarse_dst_idx = mode_coarse_idx
      if (mode_coarse_slt_idx < 0) mode_coarse_slt_idx = mode_coarse_idx

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

      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)

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
         if (trim(aername(iaer)) == 'DUST1') idxdst1 = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
         if (trim(aername(iaer)) == 'BCPHIL') idxbcphi = iaer
      end do
   end if


   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, iulog, pi, &
        mincld, nucleate_ice_subgrid)

   ! get indices for fields in the physics buffer
   ast_idx      = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_prep_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (nucleate_ice_cam_prep_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NUCLEATE_ICE_CAM_PREP_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      nucleate_ice_cam_prep_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      nucleate_ice_cam_prep_use_native_impl = .false.
   end if

   nucleate_ice_cam_prep_impl_selected = .true.

   if (masterproc) then
      if (nucleate_ice_cam_prep_use_native_impl) then
         write(iulog,*) 'nucleate_ice_cam_prep implementation = native'
      else
         write(iulog,*) 'nucleate_ice_cam_prep implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_prep_select_impl

!================================================================================================

subroutine nucleate_ice_cam_prep_log_entered()

   if (nucleate_ice_cam_prep_entered_logged) return
   nucleate_ice_cam_prep_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'nucleate_ice_cam_prep entered ' // &
           '(rho/icecldf/output zero/relhum prep = codon; qsat/nucleati/history output = native)'
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_prep_log_entered

!================================================================================================

subroutine nucleate_ice_cam_post_log_entered()

   if (nucleate_ice_cam_post_entered_logged) return
   nucleate_ice_cam_post_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'nucleate_ice_cam_post entered ' // &
           '(post-nucleati naai_hom/history rho conversion = codon; nucleati/preexisting = native)'
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_post_log_entered

!================================================================================================

subroutine nucleate_ice_cam_modal_dust_log_entered()

   if (nucleate_ice_cam_modal_dust_entered_logged) return
   nucleate_ice_cam_modal_dust_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'nucleate_ice_cam_modal_dust entered ' // &
           '(modal dust number prep = codon; sulfate erf/log and nucleati = native)'
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_modal_dust_log_entered

!================================================================================================

subroutine nucleate_ice_cam_modal_so4_log_entered()

   if (nucleate_ice_cam_modal_so4_entered_logged) return
   nucleate_ice_cam_modal_so4_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'nucleate_ice_cam_modal_so4 entered ' // &
           '(modal sulfate number prep = codon; nucleati = native)'
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_modal_so4_log_entered

!================================================================================================

subroutine nucleate_ice_cam_modal_nucleati_log_entered()

   if (nucleate_ice_cam_modal_nucleati_entered_logged) return
   nucleate_ice_cam_modal_nucleati_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'nucleate_ice_cam_calc direct = codon modal nucleation loop; ' // &
           'qsat/svp and history/native CAM API islands'
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_modal_nucleati_log_entered

!================================================================================================

subroutine nucleate_ice_cam_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine nucleate_ice_cam_log_direct

!================================================================================================

subroutine nucleate_ice_cam_calc( &
   state, wsubi, pbuf)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8), target,            intent(in)    :: wsubi(:,:)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
 
   ! local workspace

   ! naai and naai_hom are the outputs shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)

   integer :: lchnk, ncol
   integer :: itim_old
   integer :: i, k, m

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)

   real(r8), pointer :: num_accum(:,:)  ! number m.r. of accumulation mode
   real(r8), pointer :: num_aitken(:,:) ! number m.r. of aitken mode
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl
   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio
   real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius

   real(r8), pointer :: ast(:,:)
   real(r8), target :: icecldf(pcols,pver)  ! ice cloud fraction

   real(r8), target :: rho(pcols,pver)      ! air density (kg m-3)

   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8), target :: qs(pcols)    ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water
   real(r8), target :: svp_water_tair(pcols,pver)
   real(r8), target :: svp_ice_tair(pcols,pver)

   real(r8), target :: relhum(pcols,pver)  ! relative humidity
   real(r8), target :: icldm(pcols,pver)   ! ice cloud fraction

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)
   real(r8), target :: dst_num_grid(pcols,pver)      ! modal dust aerosol number (#/cm^3)
   real(r8), target :: so4_num_grid(pcols,pver)      ! modal sulfate aerosol number (#/cm^3)
   real(r8) :: wght
   real(r8) :: dmc
   real(r8) :: ssmc

   ! For pre-existing ice
   real(r8), target :: fhom(pcols,pver)    ! how much fraction of cloud can reach Shom
   real(r8), target :: wice(pcols,pver)    ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8), target :: weff(pcols,pver)    ! effective Vertical velocity for ice nucleation (m/s); weff=wsubi-wice
   real(r8), target :: INnso4(pcols,pver)  ! #/m3, so4 aerosol number used for ice nucleation
   real(r8), target :: INnbc(pcols,pver)   ! #/m3, bc aerosol number used for ice nucleation
   real(r8), target :: INndust(pcols,pver) ! #/m3, dust aerosol number used for ice nucleation
   real(r8), target :: INhet(pcols,pver)   ! #/m3, ice number from het freezing
   real(r8), target :: INhom(pcols,pver)   ! #/m3, ice number from hom freezing
   real(r8), target :: INFrehom(pcols,pver) !  hom freezing occurence frequency.  1 occur, 0 not occur.
   real(r8), target :: INFreIN(pcols,pver)  !  ice nucleation occerence frequency.   1 occur, 0 not occur.

   ! history output for ice nucleation
   real(r8), target :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
   real(r8), target :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
   real(r8), target :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
   real(r8), target :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)

   logical :: use_native_prep_impl
   logical :: use_modal_nucleati_batch
   integer(c_int64_t), target :: nucleati_warn

   interface
      subroutine nucleate_ice_cam_rho_codon(ncol_c, pcols_c, pver_c, top_lev_c, rair_c, pmid_p, t_p, rho_p) &
           bind(c, name="nucleate_ice_cam_rho_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: rair_c
         type(c_ptr), value :: pmid_p, t_p, rho_p
      end subroutine nucleate_ice_cam_rho_codon
      subroutine nucleate_ice_cam_icecldf_codon(ncol_c, pcols_c, pver_c, ast_p, icecldf_p) &
           bind(c, name="nucleate_ice_cam_icecldf_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: ast_p, icecldf_p
      end subroutine nucleate_ice_cam_icecldf_codon
      subroutine nucleate_ice_cam_zero_outputs_codon(ncol_c, pcols_c, pver_c, use_preexisting_ice_c, &
           naai_p, naai_hom_p, nihf_p, niimm_p, nidep_p, nimey_p, fhom_p, wice_p, weff_p, innso4_p, &
           innbc_p, inndust_p, inhet_p, inhom_p, infrehom_p, infrein_p) &
           bind(c, name="nucleate_ice_cam_zero_outputs_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, use_preexisting_ice_c
         type(c_ptr), value :: naai_p, naai_hom_p, nihf_p, niimm_p, nidep_p, nimey_p
         type(c_ptr), value :: fhom_p, wice_p, weff_p, innso4_p, innbc_p, inndust_p
         type(c_ptr), value :: inhet_p, inhom_p, infrehom_p, infrein_p
      end subroutine nucleate_ice_cam_zero_outputs_codon
      subroutine nucleate_ice_cam_relhum_codon(ncol_c, pcols_c, k_c, mincld_c, qn_p, qs_p, icecldf_p, relhum_p, &
           icldm_p) bind(c, name="nucleate_ice_cam_relhum_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, k_c
         real(c_double), value :: mincld_c
         type(c_ptr), value :: qn_p, qs_p, icecldf_p, relhum_p, icldm_p
      end subroutine nucleate_ice_cam_relhum_codon
      subroutine nucleate_ice_cam_post_nucleati_codon(ncol_c, pcols_c, pver_c, top_lev_c, tmelt_c, t_p, rho_p, &
           naai_hom_p, nihf_p, niimm_p, nidep_p, nimey_p) bind(c, name="nucleate_ice_cam_post_nucleati_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: tmelt_c
         type(c_ptr), value :: t_p, rho_p, naai_hom_p, nihf_p, niimm_p, nidep_p, nimey_p
      end subroutine nucleate_ice_cam_post_nucleati_codon
      subroutine nucleate_ice_cam_modal_dst_num_codon(ncol_c, pcols_c, pver_c, top_lev_c, separate_dust_c, &
           rho_p, coarse_dust_p, coarse_nacl_p, num_coarse_p, dst_num_p) &
           bind(c, name="nucleate_ice_cam_modal_dst_num_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, separate_dust_c
         type(c_ptr), value :: rho_p, coarse_dust_p, coarse_nacl_p, num_coarse_p, dst_num_p
      end subroutine nucleate_ice_cam_modal_dst_num_codon
      subroutine nucleate_ice_cam_modal_so4_num_codon(ncol_c, pcols_c, pver_c, top_lev_c, mode_aitken_idx_c, &
           tmelt_c, sigmag_aitken_c, t_p, rho_p, num_aitken_p, dgnum_p, so4_num_p) &
           bind(c, name="nucleate_ice_cam_modal_so4_num_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, mode_aitken_idx_c
         real(c_double), value :: tmelt_c, sigmag_aitken_c
         type(c_ptr), value :: t_p, rho_p, num_aitken_p, dgnum_p, so4_num_p
      end subroutine nucleate_ice_cam_modal_so4_num_codon
      subroutine nucleate_ice_cam_modal_nucleati_batch_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
           tmelt_c, use_hetfrz_classnuc_c, mincld_c, subgrid_c, ci_c, shet_c, minweff_c, gamma4_c, pi_c, &
           wsubi_p, t_p, pmid_p, relhum_p, icldm_p, qc_p, qi_p, ni_p, rho_p, so4_num_p, dst_num_p, &
           svp_water_p, svp_ice_p, naai_p, nihf_p, niimm_p, nidep_p, nimey_p, warn_p) &
           bind(c, name="nucleate_ice_cam_modal_nucleati_batch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, use_hetfrz_classnuc_c
         real(c_double), value :: tmelt_c, mincld_c, subgrid_c, ci_c, shet_c, minweff_c, gamma4_c, pi_c
         type(c_ptr), value :: wsubi_p, t_p, pmid_p, relhum_p, icldm_p, qc_p, qi_p, ni_p, rho_p
         type(c_ptr), value :: so4_num_p, dst_num_p, svp_water_p, svp_ice_p
         type(c_ptr), value :: naai_p, nihf_p, niimm_p, nidep_p, nimey_p, warn_p
      end subroutine nucleate_ice_cam_modal_nucleati_batch_codon
   end interface


   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid

   call nucleate_ice_cam_prep_select_impl()
   use_native_prep_impl = nucleate_ice_cam_prep_use_native_impl
   use_modal_nucleati_batch = (.not. use_native_prep_impl) .and. clim_modal_aero .and. (.not. use_preexisting_ice)

   if (use_native_prep_impl) then
      do k = top_lev, pver
         do i = 1, ncol
            rho(i,k) = pmid(i,k)/(rair*t(i,k))
         end do
      end do
   else
      call nucleate_ice_cam_rho_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           real(rair, c_double), c_loc(pmid(1,1)), c_loc(t(1,1)), c_loc(rho(1,1)) &
      )
   end if

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_accum_idx,  'a', state, pbuf, num_accum)
      call rad_cnst_get_mode_num(0, mode_aitken_idx, 'a', state, pbuf, num_aitken)
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

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   if (use_native_prep_impl) then
      icecldf(:ncol,:pver) = ast(:ncol,:pver)
   else
      call nucleate_ice_cam_icecldf_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           c_loc(ast(1,1)), c_loc(icecldf(1,1)) &
      )
   end if

   if (clim_modal_aero) then
      call pbuf_get_field(pbuf, dgnum_idx, dgnum)
   end if

   if ((.not. use_native_prep_impl) .and. clim_modal_aero .and. (.not. use_preexisting_ice)) then
      call nucleate_ice_cam_modal_dst_num_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           int(merge(1, 0, separate_dust), c_int64_t), c_loc(rho(1,1)), c_loc(coarse_dust(1,1)), &
           c_loc(coarse_nacl(1,1)), c_loc(num_coarse(1,1)), c_loc(dst_num_grid(1,1)) &
      )
      call nucleate_ice_cam_modal_dust_log_entered()
      call nucleate_ice_cam_modal_so4_num_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           int(mode_aitken_idx, c_int64_t), real(tmelt, c_double), real(sigmag_aitken, c_double), &
           c_loc(t(1,1)), c_loc(rho(1,1)), c_loc(num_aitken(1,1)), c_loc(dgnum(1,1,1)), c_loc(so4_num_grid(1,1)) &
      )
      call nucleate_ice_cam_modal_so4_log_entered()
   end if

   ! naai and naai_hom are the outputs from this parameterization
   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
   if (use_native_prep_impl) then
      naai(1:ncol,1:pver)     = 0._r8
      naai_hom(1:ncol,1:pver) = 0._r8

      ! initialize history output fields for ice nucleation
      nihf(1:ncol,1:pver)  = 0._r8
      niimm(1:ncol,1:pver) = 0._r8
      nidep(1:ncol,1:pver) = 0._r8
      nimey(1:ncol,1:pver) = 0._r8

      if (use_preexisting_ice) then
         fhom(:,:)     = 0.0_r8
         wice(:,:)     = 0.0_r8
         weff(:,:)     = 0.0_r8
         INnso4(:,:)   = 0.0_r8
         INnbc(:,:)    = 0.0_r8
         INndust(:,:)  = 0.0_r8
         INhet(:,:)    = 0.0_r8
         INhom(:,:)    = 0.0_r8
         INFrehom(:,:) = 0.0_r8
         INFreIN(:,:)  = 0.0_r8
      endif
   else
      call nucleate_ice_cam_zero_outputs_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(merge(1, 0, use_preexisting_ice), c_int64_t), &
           c_loc(naai(1,1)), c_loc(naai_hom(1,1)), c_loc(nihf(1,1)), c_loc(niimm(1,1)), &
           c_loc(nidep(1,1)), c_loc(nimey(1,1)), c_loc(fhom(1,1)), c_loc(wice(1,1)), &
           c_loc(weff(1,1)), c_loc(INnso4(1,1)), c_loc(INnbc(1,1)), c_loc(INndust(1,1)), &
           c_loc(INhet(1,1)), c_loc(INhom(1,1)), c_loc(INFrehom(1,1)), c_loc(INFreIN(1,1)) &
      )
   end if

   do k = top_lev, pver

      ! Get humidity and saturation vapor pressures
      call qsat_water(t(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      if (use_native_prep_impl) then
         do i = 1, ncol

            relhum(i,k) = qn(i,k)/qs(i)

            ! get cloud fraction, check for minimum
            icldm(i,k) = max(icecldf(i,k), mincld)

         end do
      else
         call nucleate_ice_cam_relhum_codon( &
              int(ncol, c_int64_t), int(pcols, c_int64_t), int(k, c_int64_t), real(mincld, c_double), &
              c_loc(qn(1,1)), c_loc(qs(1)), c_loc(icecldf(1,1)), c_loc(relhum(1,1)), c_loc(icldm(1,1)) &
         )
      end if
      if (use_modal_nucleati_batch) then
         do i = 1, ncol
            svp_water_tair(i,k) = svp_water(t(i,k))
            svp_ice_tair(i,k)   = svp_ice(t(i,k))
         end do
      end if
   end do

   if (.not. use_native_prep_impl) then
      call nucleate_ice_cam_prep_log_entered()
   end if


   if (use_modal_nucleati_batch) then
      nucleati_warn = 0_c_int64_t
      call nucleate_ice_cam_modal_nucleati_batch_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           real(tmelt, c_double), int(merge(1, 0, use_hetfrz_classnuc), c_int64_t), &
           real(mincld, c_double), real(nucleate_ice_subgrid, c_double), real(0.5e3_r8*pi/6._r8, c_double), &
           real(1.3_r8, c_double), real(0.001_r8, c_double), real(6.0_r8, c_double), real(pi, c_double), &
           c_loc(wsubi(1,1)), c_loc(t(1,1)), c_loc(pmid(1,1)), c_loc(relhum(1,1)), c_loc(icldm(1,1)), &
           c_loc(qc(1,1)), c_loc(qi(1,1)), c_loc(ni(1,1)), c_loc(rho(1,1)), c_loc(so4_num_grid(1,1)), &
           c_loc(dst_num_grid(1,1)), c_loc(svp_water_tair(1,1)), c_loc(svp_ice_tair(1,1)), c_loc(naai(1,1)), &
           c_loc(nihf(1,1)), c_loc(niimm(1,1)), c_loc(nidep(1,1)), c_loc(nimey(1,1)), c_loc(nucleati_warn) &
      )
      call nucleate_ice_cam_modal_nucleati_log_entered()
      if (nucleati_warn /= 0_c_int64_t) then
         write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      end if
   else
      do k = top_lev, pver
         do i = 1, ncol

            if (t(i,k) < tmelt - 5._r8) then

               ! compute aerosol number for so4, soot, and dust with units #/cm^3
               so4_num  = 0._r8
               soot_num = 0._r8
               dst1_num = 0._r8
               dst2_num = 0._r8
               dst3_num = 0._r8
               dst4_num = 0._r8
               dst_num  = 0._r8

            if (clim_modal_aero) then
               if ((.not. use_native_prep_impl) .and. (.not. use_preexisting_ice)) then
                  !For modal aerosols, assume for the upper troposphere:
                  ! soot = accumulation mode
                  ! sulfate = aiken mode
                  ! dust = coarse mode
                  ! since modal has internal mixtures.
                  soot_num = num_accum(i,k)*rho(i,k)*1.0e-6_r8
                  dst_num  = dst_num_grid(i,k)
               else
                  !For modal aerosols, assume for the upper troposphere:
                  ! soot = accumulation mode
                  ! sulfate = aiken mode
                  ! dust = coarse mode
                  ! since modal has internal mixtures.
                  soot_num = num_accum(i,k)*rho(i,k)*1.0e-6_r8
                  dmc  = coarse_dust(i,k)*rho(i,k)
                  ssmc = coarse_nacl(i,k)*rho(i,k)

                  if (dmc > 0._r8) then
                     if ( separate_dust ) then
                        ! 7-mode -- has separate dust and seasalt mode types and
                        !           no need for weighting
                        wght = 1._r8
                     else
                        ! 3-mode -- needs weighting for dust since dust and seasalt
                        !           are combined in the "coarse" mode type
                        wght = dmc/(ssmc + dmc)
                     endif
                     dst_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
                  else
                     dst_num = 0.0_r8
                  end if
               end if

               if ((.not. use_native_prep_impl) .and. (.not. use_preexisting_ice)) then
                  so4_num = so4_num_grid(i,k)
               else
                  if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
                     if (.not. use_preexisting_ice) then
                        ! only allow so4 with D>0.1 um in ice nucleation
                        so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8 &
                           * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/dgnum(i,k,mode_aitken_idx))/  &
                           (2._r8**0.5_r8*log(sigmag_aitken))))
                     else
                        ! all so4 from aitken
                        so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8
                     end if
                  else
                     so4_num = 0.0_r8
                  end if
                  so4_num = max(0.0_r8, so4_num)
               end if

            else

               if (idxsul > 0) then 
                  so4_num = naer2(i,k,idxsul)/25._r8 *1.0e-6_r8
               end if
               if (idxbcphi > 0) then 
                  soot_num = naer2(i,k,idxbcphi)/25._r8 *1.0e-6_r8
               end if
               if (idxdst1 > 0) then 
                  dst1_num = naer2(i,k,idxdst1)/25._r8 *1.0e-6_r8
               end if
               if (idxdst2 > 0) then 
                  dst2_num = naer2(i,k,idxdst2)/25._r8 *1.0e-6_r8
               end if
               if (idxdst3 > 0) then 
                  dst3_num = naer2(i,k,idxdst3)/25._r8 *1.0e-6_r8
               end if
               if (idxdst4 > 0) then 
                  dst4_num = naer2(i,k,idxdst4)/25._r8 *1.0e-6_r8
               end if
               dst_num = dst1_num + dst2_num + dst3_num + dst4_num

            end if

            ! *** Turn off soot nucleation ***
            soot_num = 0.0_r8

            call nucleati( &
               wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
               qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
               so4_num, dst_num, soot_num,                               &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               wice(i,k), weff(i,k), fhom(i,k))

            if (use_native_prep_impl .or. use_preexisting_ice) then
               naai_hom(i,k) = nihf(i,k)

               ! output activated ice (convert from #/kg -> #/m3)
               nihf(i,k)     = nihf(i,k) *rho(i,k)
               niimm(i,k)    = niimm(i,k)*rho(i,k)
               nidep(i,k)    = nidep(i,k)*rho(i,k)
               nimey(i,k)    = nimey(i,k)*rho(i,k)

               if (use_preexisting_ice) then
                  INnso4(i,k) =so4_num*1e6_r8  ! (convert from #/cm3 -> #/m3)
                  INnbc(i,k)  =soot_num*1e6_r8
                  INndust(i,k)=dst_num*1e6_r8
                  INFreIN(i,k)=1.0_r8          ! 1,ice nucleation occur
                  INhet(i,k) = niimm(i,k) + nidep(i,k)   ! #/m3, nimey not in cirrus
                  INhom(i,k) = nihf(i,k)                 ! #/m3
                  if (INhom(i,k).gt.1e3_r8)   then ! > 1/L
                     INFrehom(i,k)=1.0_r8       ! 1, hom freezing occur
                  endif

                  ! exclude  no ice nucleaton
                  if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then
                     INnso4(i,k) =0.0_r8
                     INnbc(i,k)  =0.0_r8
                     INndust(i,k)=0.0_r8
                     INFreIN(i,k)=0.0_r8
                     INhet(i,k) = 0.0_r8
                     INhom(i,k) = 0.0_r8
                     INFrehom(i,k)=0.0_r8
                     wice(i,k) = 0.0_r8
                     weff(i,k) = 0.0_r8
                     fhom(i,k) = 0.0_r8
                  endif
               end if
            end if

            end if
         end do
      end do
   end if

   if ((.not. use_native_prep_impl) .and. (.not. use_preexisting_ice)) then
      call nucleate_ice_cam_post_nucleati_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           real(tmelt, c_double), c_loc(t(1,1)), c_loc(rho(1,1)), c_loc(naai_hom(1,1)), &
           c_loc(nihf(1,1)), c_loc(niimm(1,1)), c_loc(nidep(1,1)), c_loc(nimey(1,1)) &
      )
      call nucleate_ice_cam_post_log_entered()
   end if

   if (.not. clim_modal_aero) then

      deallocate( &
         naer2,    &
         maerosol)

   end if

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)

   if (use_preexisting_ice) then
      call outfld( 'fhom' , fhom, pcols, lchnk)
      call outfld( 'WICE' , wice, pcols, lchnk)
      call outfld( 'WEFF' , weff, pcols, lchnk)
      call outfld('INnso4  ',INnso4 , pcols,lchnk)
      call outfld('INnbc   ',INnbc  , pcols,lchnk)
      call outfld('INndust ',INndust, pcols,lchnk)
      call outfld('INhet   ',INhet  , pcols,lchnk)
      call outfld('INhom   ',INhom  , pcols,lchnk)
      call outfld('INFrehom',INFrehom,pcols,lchnk)
      call outfld('INFreIN ',INFreIN, pcols,lchnk)
   end if

end subroutine nucleate_ice_cam_calc

!================================================================================================

end module nucleate_ice_cam

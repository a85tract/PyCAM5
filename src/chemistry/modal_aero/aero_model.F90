!===============================================================================
! Modal Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use perf_mod,       only: t_startf, t_stopf
  use camsrfexch,     only: cam_in_t, cam_out_t
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use physics_buffer, only: pbuf_get_field, pbuf_get_index, pbuf_set_field
  use physconst,      only: gravit, rair, rhoh2o
  use mo_util,        only: chemistry_misc_codon_touch
  use spmd_utils,     only: masterproc

  use cam_history,    only: outfld, fieldname_len
  use chem_mods,      only: gas_pcnst, adv_mass
  use mo_tracname,    only: solsym

  use modal_aero_data,only: cnst_name_cw
  use modal_aero_data,only: ntot_amode, modename_amode, nspec_amode_max
  use ref_pres,       only: top_lev => clim_modal_aero_top_lev

  use modal_aero_wateruptake, only: modal_strat_sulfate

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea  ! tropopspheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea ! stratospheric aerosol dry surface area for chemistry

 ! Misc private data 

  ! number of modes
  integer :: nmodes
  integer :: pblh_idx            = 0
  integer :: dgnum_idx           = 0
  integer :: dgnumwet_idx        = 0
  integer :: rate1_cw2pr_st_idx  = 0  

  integer :: wetdens_ap_idx      = 0
  integer :: qaerwat_idx         = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: nevapr_idx          = 0
  integer :: rprddp_idx          = 0 
  integer :: rprdsh_idx          = 0 
  integer :: sulfeq_idx = -1
  integer, parameter :: drydep_mode_phase_nslot = nspec_amode_max + 2
  integer, parameter :: wetdep_mode_phase_nslot = nspec_amode_max + 2

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8) :: scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
  real(r8) :: scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)

  ! for surf_area_dens 
  integer :: num_idx(ntot_amode) = -1
  integer :: index_tot_mass(ntot_amode,10) = -1
  integer :: index_chm_mass(ntot_amode,10) = -1

  integer :: ndx_h2so4
  character(len=fieldname_len) :: dgnum_name(ntot_amode), dgnumwet_name(ntot_amode)

  ! Namelist variables
  character(len=16) :: wetdep_list(pcnst) = ' '
  character(len=16) :: drydep_list(pcnst) = ' '
  real(r8)          :: sol_facti_cloud_borne   = 1._r8
  real(r8)          :: sol_factb_interstitial  = 0.1_r8
  real(r8)          :: sol_factic_interstitial = 0.4_r8

  integer :: ndrydep = 0
  integer,allocatable :: drydep_indices(:)
  integer :: nwetdep = 0
  integer,allocatable :: wetdep_indices(:)
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)

  logical :: modal_accum_coarse_exch = .false.
  logical :: aero_model_drydep_use_native_impl = .false.
  logical :: aero_model_drydep_impl_selected = .false.
  logical :: aero_model_drydep_proof_written = .false.
  logical :: aero_model_drydep_prepare_shell_proof_written = .false.
  logical :: aero_model_drydep_fullshell_wrap_proof_written = .false.
  integer :: aero_model_drydep_branch_mask = 0
  logical :: aero_model_drydep_branch_selected = .false.
  logical :: aero_model_wetdep_use_native_impl = .false.
  logical :: aero_model_wetdep_impl_selected = .false.
  logical :: aero_model_wetdep_proof_written = .false.
  logical :: aero_model_wetdep_wrap_proof_written = .false.
  logical :: aero_model_wetdep_fullshell_wrap_proof_written = .false.
  logical :: aero_model_wetdep_mode_phase_wrap_proof_written = .false.
  logical :: aero_model_wetdep_mode_phase_use_interstitial = .true.
  logical :: aero_model_wetdep_mode_phase_use_cloudborne = .true.
  logical :: aero_model_wetdep_mode_phase_selected = .false.
  logical :: aero_model_gasaerexch_use_native_impl = .false.
  logical :: aero_model_gasaerexch_impl_selected = .false.
  logical :: aero_model_gasaerexch_proof_written = .false.
  logical :: aero_model_gasaerexch_wrap_proof_written = .false.
  logical :: aero_model_gasaerexch_load_snapshot_proof_written = .false.
  logical :: aero_model_gasaerexch_presetsox_proof_written = .false.
  logical :: aero_model_gasaerexch_store_snapshot_proof_written = .false.
  logical :: aero_model_gasaerexch_aq_save_proof_written = .false.
  logical :: aero_model_gasaerexch_column_flux_use_native_impl = .false.
  logical :: aero_model_gasaerexch_column_flux_impl_selected = .false.
  logical :: aero_model_gasaerexch_h2so4_save_use_native_impl = .false.
  logical :: aero_model_gasaerexch_h2so4_save_impl_selected = .false.
  logical :: aero_model_gasaerexch_h2so4_delta_use_native_impl = .false.
  logical :: aero_model_gasaerexch_h2so4_delta_impl_selected = .false.
  logical :: aero_model_gasaerexch_gas_tend_use_native_impl = .false.
  logical :: aero_model_gasaerexch_gas_tend_impl_selected = .false.
  logical :: aero_model_gasaerexch_aq_tend_use_native_impl = .false.
  logical :: aero_model_gasaerexch_aq_tend_impl_selected = .false.
  logical :: aero_model_emissions_use_native_impl = .false.
  logical :: aero_model_emissions_impl_selected = .false.
  logical :: aero_model_emissions_proof_written = .false.
  logical :: aero_model_emissions_wrap_proof_written = .false.
  logical :: aero_model_emissions_dust_stage_proof_written = .false.
  logical :: aero_model_emissions_seasalt_stage_proof_written = .false.
  logical :: aero_model_emissions_all_stage_proof_written = .false.
  logical :: aero_model_emissions_seasalt_wind_use_native_impl = .false.
  logical :: aero_model_emissions_seasalt_wind_impl_selected = .false.
  logical :: aero_model_emissions_accumulate_sflx_use_native_impl = .false.
  logical :: aero_model_emissions_accumulate_sflx_impl_selected = .false.
  logical :: qqcw2vmr_use_native_impl = .false.
  logical :: qqcw2vmr_impl_selected = .false.
  logical :: vmr2qqcw_use_native_impl = .false.
  logical :: vmr2qqcw_impl_selected = .false.
  logical :: aero_model_register_codon_logged = .false.
  logical :: aero_model_strat_surfarea_codon_logged = .false.

contains
  
  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use iso_c_binding,   only: c_int64_t, c_loc

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr, n, status, l, j, code
    character(len=*), parameter :: subname = 'aero_model_readnl'
    character(len=32) :: impl_name
    logical :: use_native
    integer(c_int64_t) :: status_c
    integer(c_int64_t), target :: aer_wetdep_ascii(16,pcnst)
    integer(c_int64_t), target :: aer_drydep_ascii(16,pcnst)
    integer(c_int64_t), target :: wetdep_ascii(16,pcnst)
    integer(c_int64_t), target :: drydep_ascii(16,pcnst)

    ! Namelist variables
    character(len=16) :: aer_wetdep_list(pcnst) = ' '
    character(len=16) :: aer_drydep_list(pcnst) = ' '

    interface
       function aero_model_readnl_codon(pcnst_c, name_len_c, aer_wetdep_list_p, &
            aer_drydep_list_p, wetdep_list_p, drydep_list_p) result(out_c) &
            bind(c, name="aero_model_readnl_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcnst_c, name_len_c
         type(c_ptr), value :: aer_wetdep_list_p, aer_drydep_list_p
         type(c_ptr), value :: wetdep_list_p, drydep_list_p
         integer(c_int64_t) :: out_c
       end function aero_model_readnl_codon
    end interface

    namelist /aerosol_nl/ aer_wetdep_list, aer_drydep_list, sol_facti_cloud_borne, &
       sol_factb_interstitial, sol_factic_interstitial, modal_strat_sulfate, modal_accum_coarse_exch

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
    call mpibcast(sol_factb_interstitial, 1,                        mpir8,   0, mpicom)
    call mpibcast(sol_factic_interstitial, 1,                       mpir8,   0, mpicom)
    call mpibcast(modal_strat_sulfate,     1,                       mpilog,  0, mpicom)
    call mpibcast(modal_accum_coarse_exch, 1,                       mpilog,  0, mpicom)
#endif

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_READNL_IMPL', impl_name, n, status)
    use_native = .false.
    if (status == 0 .and. n > 0) then
       do l = 1, n
          code = iachar(impl_name(l:l))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(l:l) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native = trim(adjustl(impl_name(:n))) == 'native'
    end if

    if (.not. use_native) then
       do l = 1, pcnst
          do j = 1, len(aer_wetdep_list(1))
             aer_wetdep_ascii(j,l) = int(iachar(aer_wetdep_list(l)(j:j)), c_int64_t)
             aer_drydep_ascii(j,l) = int(iachar(aer_drydep_list(l)(j:j)), c_int64_t)
          end do
       end do

       status_c = aero_model_readnl_codon(int(pcnst, c_int64_t), &
            int(len(aer_wetdep_list(1)), c_int64_t), c_loc(aer_wetdep_ascii(1,1)), &
            c_loc(aer_drydep_ascii(1,1)), c_loc(wetdep_ascii(1,1)), c_loc(drydep_ascii(1,1)))
       if (status_c /= 1_c_int64_t) then
          call endrun('aero_model_readnl_codon failed')
       end if

       do l = 1, pcnst
          wetdep_list(l) = ' '
          drydep_list(l) = ' '
          do j = 1, len(wetdep_list(1))
             wetdep_list(l)(j:j) = achar(int(wetdep_ascii(j,l)))
             drydep_list(l)(j:j) = achar(int(drydep_ascii(j,l)))
          end do
       end do

       if (masterproc) then
          write(iulog,'(A)') 'aero_model_readnl implementation = codon'
          write(iulog,'(A)') 'aero_model_readnl direct = codon; namelist I/O and MPI broadcast native boundary'
          call flush(iulog)
       end if
       return
    end if

    wetdep_list = aer_wetdep_list
    drydep_list = aer_drydep_list

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register
    use iso_c_binding, only : c_int64_t
    use modal_aero_initialize_data, only : modal_aero_register

    integer(c_int64_t) :: active_c

    interface
       function aero_model_register_codon(stage_c) result(out_c) bind(c, name="aero_model_register_codon")
         import :: c_int64_t
         integer(c_int64_t), value :: stage_c
         integer(c_int64_t) :: out_c
       end function aero_model_register_codon
    end interface

    active_c = aero_model_register_codon(1_c_int64_t)
    if (.not. aero_model_register_codon_logged) then
       aero_model_register_codon_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'aero_model_register direct = codon; modal_aero_register native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call modal_aero_register()

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, add_default, phys_decomp
    use phys_control,    only: phys_getopts
    use mo_chem_utls,    only: get_rxt_ndx, get_spc_ndx
    use modal_aero_data, only: cnst_name_cw
    use modal_aero_initialize_data, only: modal_aero_initialize
    use rad_constituents,           only: rad_cnst_get_info
    use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    use drydep_mod,      only: inidrydep
    use wetdep,          only: wetdep_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, id
    character(len=20) :: dummy

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies

    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes

    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 
    integer :: errcode
    character(len=fieldname_len) :: field_name

    character(len=32) :: spec_name
    character(len=32) :: spec_type
    character(len=32) :: mode_type
    integer :: nspec

    call chemistry_misc_codon_touch('aero_model_init', 176)

    dgnum_idx    = pbuf_get_index('DGNUM')
    dgnumwet_idx = pbuf_get_index('DGNUMWET')
    fracis_idx   = pbuf_get_index('FRACIS') 
    prain_idx    = pbuf_get_index('PRAIN')  
    nevapr_idx   = pbuf_get_index('NEVAPR') 
    rprddp_idx   = pbuf_get_index('RPRDDP')  
    rprdsh_idx   = pbuf_get_index('RPRDSH')  
    sulfeq_idx   = pbuf_get_index('MAMH2SO4EQ',errcode)
    
    call phys_getopts( history_aerosol_out=history_aerosol )

    call rad_cnst_get_info(0, nmodes=nmodes)

    call modal_aero_initialize(pbuf2d,modal_accum_coarse_exch)
    call modal_aero_bcscavcoef_init()

    call dust_init()
    call seasalt_init()
    call wetdep_init()


    nwetdep = 0
    ndrydep = 0

    count_species: do m = 1,pcnst
       if ( len_trim(wetdep_list(m)) /= 0 ) then
          nwetdep = nwetdep+1
       endif
       if ( len_trim(drydep_list(m)) /= 0 ) then
          ndrydep = ndrydep+1
       endif
    enddo count_species
    
    if (nwetdep>0) &
         allocate(wetdep_indices(nwetdep))
    if (ndrydep>0) &
         allocate(drydep_indices(ndrydep))

    do m = 1,ndrydep
       call cnst_get_ind ( drydep_list(m), id, abort=.false. )
       if (id>0) then
          drydep_indices(m) = id
       else
          call endrun(subrname//': invalid drydep species: '//trim(drydep_list(m)) )
       endif

       if (masterproc) then
          write(iulog,*) subrname//': '//drydep_list(m)//' will have drydep applied'
       endif
    enddo
    do m = 1,nwetdep
       call cnst_get_ind ( wetdep_list(m), id, abort=.false. )
       if (id>0) then
          wetdep_indices(m) = id
       else
          call endrun(subrname//': invalid wetdep species: '//trim(wetdep_list(m)) )
       endif
       
       if (masterproc) then
          write(iulog,*) subrname//': '//wetdep_list(m)//' will have wet removal'
       endif
    enddo

    if (ndrydep>0) then

       call inidrydep(rair, gravit)

       dummy = 'RAM1'
       call addfld (dummy,'frac ',1, 'A','RAM1',phys_decomp)
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = 'airFV'
       call addfld (dummy,'frac ',1, 'A','FV',phys_decomp)
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (dust_active) then
       ! emissions diagnostics ....

       do m = 1, dust_nbin+dust_nnum
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,'kg/m2/s ',1, 'A',trim(dust_names(m))//' dust surface emission',phys_decomp)
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,'kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,'frac ',1, 'A','Soil erodibility factor',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (seasalt_active) then
       
       dummy = 'SSTSFMBL'
       call addfld (dummy,'kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       do m = 1, seasalt_nbin
          dummy = trim(seasalt_names(m)) // 'SF'
          call addfld (dummy,'kg/m2/s ',1, 'A',trim(seasalt_names(m))//' seasalt surface emission',phys_decomp)
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

    endif

    
    ! set flags for drydep tendencies
    drydep_lq(:) = .false.
    do m=1,ndrydep 
       id = drydep_indices(m)
       drydep_lq(id) =  .true.
    enddo

    ! set flags for wetdep tendencies
    wetdep_lq(:) = .false.
    do m=1,nwetdep
       id = wetdep_indices(m)
       wetdep_lq(id) = .true.
    enddo

    wetdens_ap_idx = pbuf_get_index('WETDENS_AP')
    qaerwat_idx    = pbuf_get_index('QAERWAT')
    pblh_idx       = pbuf_get_index('pblh')

    rate1_cw2pr_st_idx  = pbuf_get_index('RATE1_CW2PR_ST') 
    call pbuf_set_field(pbuf2d, rate1_cw2pr_st_idx, 0.0_r8)

    do m = 1,ndrydep
       
       ! units 
       if (drydep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(drydep_list(m))//'DDF',unit_basename//'/m2/s ',   1, 'A', &
            trim(drydep_list(m))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
       call addfld (trim(drydep_list(m))//'TBF',unit_basename//'/m2/s',   1, 'A', &
            trim(drydep_list(m))//' turbulent dry deposition flux',phys_decomp)
       call addfld (trim(drydep_list(m))//'GVF',unit_basename//'/m2/s ',   1, 'A', &
            trim(drydep_list(m))//' gravitational dry deposition flux',phys_decomp)
       call addfld (trim(drydep_list(m))//'DTQ',unit_basename//'/kg/s ',pver, 'A', &
            trim(drydep_list(m))//' dry deposition',phys_decomp)
       call addfld (trim(drydep_list(m))//'DDV','m/s     ',pver, 'A', &
            trim(drydep_list(m))//' deposition velocity',phys_decomp)

       if ( history_aerosol ) then 
          call add_default (trim(drydep_list(m))//'DDF', 1, ' ')
          call add_default (trim(drydep_list(m))//'TBF', 1, ' ')
          call add_default (trim(drydep_list(m))//'GVF', 1, ' ')
       endif

    enddo

    do m = 1,nwetdep
       
       ! units 
       if (wetdep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(wetdep_list(m))//'SFWET',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSIC',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (incloud, convective) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSIS',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (incloud, stratiform) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSBC',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (belowcloud, convective) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSBS',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (belowcloud, stratiform) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'WET',unit_basename//'/kg/s ',pver, 'A','wet deposition tendency',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SIC',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' ic wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SIS',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' is wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SBC',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' bc wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SBS',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' bs wet deposition',phys_decomp)
       
       if ( history_aerosol ) then          
          call add_default (trim(wetdep_list(m))//'SFWET', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSIC', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSIS', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSBC', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSBS', 1, ' ')
       endif

    enddo

    do m = 1,gas_pcnst

       if  ( solsym(m)(1:3) == 'num') then
          unit_basename = ' 1'  ! Units 'kg' or '1' 
       else
          unit_basename = 'kg'  ! Units 'kg' or '1' 
       end if

       call addfld( 'GS_'//trim(solsym(m)), unit_basename//'/m2/s ',1,  'A', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)', phys_decomp)
       call addfld( 'AQ_'//trim(solsym(m)), unit_basename//'/m2/s ',1,  'A', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)', phys_decomp)
       if ( history_aerosol ) then 
          call add_default( 'GS_'//trim(solsym(m)), 1, ' ')
          call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
       endif
    enddo
    do n = 1,pcnst
       if( .not. (cnst_name_cw(n) == ' ') ) then

          if (cnst_name_cw(n)(1:3) == 'num') then
             unit_basename = ' 1'
          else
             unit_basename = 'kg'  
          endif

          call addfld( cnst_name_cw(n),                unit_basename//'/kg ', pver, 'A', &
               trim(cnst_name_cw(n))//' in cloud water',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFWET', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSIC', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, convective) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSIS', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, stratiform) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSBC', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, convective) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSBS', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, stratiform) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'DDF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'TBF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' turbulent dry deposition flux',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'GVF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' gravitational dry deposition flux',phys_decomp)     

          if ( history_aerosol ) then 
             call add_default( cnst_name_cw(n), 1, ' ' )
             call add_default (trim(cnst_name_cw(n))//'GVF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFWET', 1, ' ') 
             call add_default (trim(cnst_name_cw(n))//'TBF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'DDF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSBS', 1, ' ')      
             call add_default (trim(cnst_name_cw(n))//'SFSIC', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSBC', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSIS', 1, ' ')
          endif
       endif
    enddo
    do n=1,ntot_amode
       dgnum_name(n) = ' '
       dgnumwet_name(n) = ' '
       write(dgnum_name(n),fmt='(a,i1)') 'dgnum',n
       write(dgnumwet_name(n),fmt='(a,i1)') 'dgnumwet',n
       call addfld( dgnum_name(n), 'm', pver, 'I', 'Aerosol mode dry diameter', phys_decomp )
       call addfld( dgnumwet_name(n), 'm', pver, 'I', 'Aerosol mode wet diameter', phys_decomp )
       if ( history_aerosol ) then 
          call add_default( dgnum_name(n), 1, ' ' )
          call add_default( dgnumwet_name(n), 1, ' ' )
       endif
      
       if (modal_strat_sulfate) then
          field_name = ' '
          write(field_name,fmt='(a,i1)') 'wtpct_a',n
          call addfld( field_name, '%', pver, 'I', 'Aerosol mode weight percent H2SO4', phys_decomp )
          if ( history_aerosol ) then 
             call add_default (field_name, 0, 'I')
          endif

          field_name = ' '
          write(field_name,fmt='(a,i1)') 'sulfeq_a',n
          call addfld( field_name, 'kg/kg', pver, 'I', 'H2SO4 equilibrium mixing ratio', phys_decomp )
          if ( history_aerosol ) then 
             call add_default (field_name, 0, 'I')
          endif

          field_name = ' '
          write(field_name,fmt='(a,i1)') 'sulden_a',n
          call addfld( field_name, 'g/cm3', pver, 'I', 'Sulfate aerosol particle mass density', phys_decomp )
          if ( history_aerosol ) then 
             call add_default (field_name, 0, 'I')
          endif

       end if
    end do

    ndx_h2so4 = get_spc_ndx('H2SO4')

    ! for aero_model_surfarea called from mo_usrrxt
    do l=1,ntot_amode
       test_name = ' '
       write(test_name,fmt='(a5,i1)') 'num_a',l
       num_idx(l) = get_spc_ndx( trim(test_name) )
       if (num_idx(l) < 0) then
          write(errmes,fmt='(a,i1)') 'usrrxt_inti: cannot find MAM num_idx ',l
          write(iulog,*) errmes
          call endrun(errmes)
       endif
    end do
    
    ! for surf_area_dens 
    ! define indeces associated with the various aerosol types    
    do n = 1,nmodes
       call rad_cnst_get_info(0, n, mode_type=mode_type, nspec=nspec)
       if ( trim(mode_type) /= 'primary_carbon') then ! ignore the primary_carbon mode
          do l = 1, nspec
             call rad_cnst_get_info(0, n, l, spec_type=spec_type, spec_name=spec_name)
             index_tot_mass(n,l) = get_spc_ndx(spec_name)
             if ( trim(spec_type) == 'sulfate'   .or. &
                  trim(spec_type) == 's-organic' .or. &
                  trim(spec_type) == 'black-c'   .or. &
                  trim(spec_type) == 'ammonium') then
                index_chm_mass(n,l) = get_spc_ndx(spec_name)
             endif
          enddo
       endif
    enddo

  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use dust_sediment_mod, only: dust_sediment_tend
    use drydep_mod,        only: d3ddflux, calcram
    use mo_drydep,         only: fraction_landuse
    use modal_aero_data,   only: qqcw_get_field
    use modal_aero_data,   only: qqcw_fill_cptrs
    use modal_aero_data,   only: cnst_name_cw
    use modal_aero_data,   only: alnsg_amode
    use modal_aero_data,   only: sigmag_amode
    use modal_aero_data,   only: nspec_amode
    use modal_aero_data,   only: numptr_amode
    use modal_aero_data,   only: numptrcw_amode
    use modal_aero_data,   only: lmassptr_amode
    use modal_aero_data,   only: lmassptrcw_amode
    use modal_aero_deposition, only: set_srf_drydep
    use iso_c_binding, only: c_double, c_ptr, c_int64_t, c_loc

  ! args 
    type(physics_state), target, intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)          
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt             ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

    integer :: lchnk                   ! chunk identifier
    integer :: ncol                    ! number of atmospheric columns
    integer :: jvlc                    ! index for last dimension of vlc_xxx arrays
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: lspec                   ! index for aerosol number / chem-mass / water-mass
    integer :: m                       ! aerosol mode index
    integer :: mm                      ! tracer index
    integer :: i

    real(r8) :: tvs(pcols,pver)
    real(r8), target :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: sflx(pcols)          ! deposition flux
    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp) ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8), target :: rad_drop(pcols,pver)
    real(r8), target :: dens_drop(pcols,pver)
    real(r8), target :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity
    real(r8) :: vlc_dry_full(pcols,pver,4,ntot_amode)
    real(r8) :: vlc_grv_full(pcols,pver,4,ntot_amode)
    real(r8) :: vlc_trb_full(pcols,4,ntot_amode)
    real(r8), target :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8), target :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)
    type(c_ptr) :: qqcw_ptrs(pcnst)
    integer(c_int64_t) :: drydep_slot_active(drydep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t) :: drydep_slot_mm(drydep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t) :: drydep_slot_jvlc(drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_qqcw_mode_phase(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_diag_ddv(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_diag_dqdt(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_diag_sflx(pcols,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_diag_dep_trb(pcols,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: drydep_diag_dep_grv(pcols,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)
    logical  :: apply_srf_drydep_local
    integer(c_int64_t), target :: branch_mask_c

    interface
       subroutine aero_model_drydep_codon(apply_srf_drydep_c, branch_mask_p, &
            ncol_c, pcols_c, pver_c, pcnst_c, rair_c, rhoh2o_c, state_t_p, state_pmid_p, &
            rho_p, rad_drop_p, dens_drop_p, sg_drop_p, aerdepdryis_p, aerdepdrycw_p) &
            bind(c, name="aero_model_drydep_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: apply_srf_drydep_c
         type(c_ptr), value :: branch_mask_p
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c
         real(c_double), value :: rair_c, rhoh2o_c
         type(c_ptr), value :: state_t_p, state_pmid_p, rho_p, rad_drop_p, dens_drop_p, sg_drop_p
         type(c_ptr), value :: aerdepdryis_p, aerdepdrycw_p
       end subroutine aero_model_drydep_codon
    end interface

    call aero_model_drydep_select_impl()
    if (.not. aero_model_drydep_use_native_impl) then
       apply_srf_drydep_local = .false.
    else
       apply_srf_drydep_local = .not. aerodep_flx_prescribed()
    end if

    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
                  state%pdel(:,pver),fvin,fv)

    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )
 
    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species
    
    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', lq=drydep_lq)

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 

    if (.not. aero_model_drydep_use_native_impl) then
       call aero_model_drydep_codon( &
            merge(1_c_int64_t, 0_c_int64_t, .not. aerodep_flx_prescribed()), c_loc(branch_mask_c), &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
            real(rair, c_double), real(rhoh2o, c_double), c_loc(state%t), c_loc(state%pmid), c_loc(rho), &
            c_loc(rad_drop), c_loc(dens_drop), c_loc(sg_drop), c_loc(aerdepdryis), c_loc(aerdepdrycw) &
       )
       aero_model_drydep_branch_mask = int(branch_mask_c)
       aero_model_drydep_branch_selected = .true.
       apply_srf_drydep_local = iand(aero_model_drydep_branch_mask, 1) /= 0
       if (masterproc .and. .not. aero_model_drydep_prepare_shell_proof_written) then
          write(iulog,'(A)') 'aero_model_drydep init shell entered (unified branch/rho/drop/aerdep stage dispatch = codon)'
          call aero_model_drydep_append_impl_proof('AERO_MODEL_DRYDEP_PROOF_FILE', &
               'aero_model_drydep init shell entered (unified branch/rho/drop/aerdep stage dispatch = codon)')
          aero_model_drydep_prepare_shell_proof_written = .true.
          call flush(iulog)
       end if
    else
       tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
       rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
       rad_drop(:,:) = 5.0e-6_r8
       dens_drop(:,:) = rhoh2o
       sg_drop(:,:) = 1.46_r8
       aerdepdryis(:,:) = 0._r8
       aerdepdrycw(:,:) = 0._r8
    end if

!
! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
!
! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    jvlc = 3
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 0, fraction_landuse(:,:,lchnk))
    jvlc = 4
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, fraction_landuse(:,:,lchnk))

    if (.not. aero_model_drydep_use_native_impl) then
       do m = 1, ntot_amode
          rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,m) * exp(1.5_r8*(alnsg_amode(m)**2))
          dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
          sg_aer(1:ncol,:) = sigmag_amode(m)

          call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv, &
               vlc_dry_full(:,:,1,m), vlc_trb_full(:,1,m), vlc_grv_full(:,:,1,m), &
               rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, fraction_landuse(:,:,lchnk))
          call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv, &
               vlc_dry_full(:,:,2,m), vlc_trb_full(:,2,m), vlc_grv_full(:,:,2,m), &
               rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, fraction_landuse(:,:,lchnk))

          vlc_dry_full(:,:,3,m) = vlc_dry(:,:,3)
          vlc_dry_full(:,:,4,m) = vlc_dry(:,:,4)
          vlc_grv_full(:,:,3,m) = vlc_grv(:,:,3)
          vlc_grv_full(:,:,4,m) = vlc_grv(:,:,4)
          vlc_trb_full(:,3,m) = vlc_trb(:,3)
          vlc_trb_full(:,4,m) = vlc_trb(:,4)
       end do

       call qqcw_fill_cptrs(pbuf, qqcw_ptrs)
       call aero_model_drydep_fullshell_codon_wrap( &
            ncol, dt, state%pint, state%pdel, rho, vlc_dry_full, vlc_trb_full, vlc_grv_full, state%q, ptend%q, &
            qqcw_ptrs, drydep_slot_active, drydep_slot_mm, drydep_slot_jvlc, drydep_qqcw_mode_phase, &
            drydep_diag_ddv, drydep_diag_dqdt, drydep_diag_sflx, drydep_diag_dep_trb, drydep_diag_dep_grv )

       do m = 1, ntot_amode
          do lphase = 1, 2
             do lspec = 1, drydep_mode_phase_nslot
                if (drydep_slot_active(lspec,lphase,m) == 0_c_int64_t) cycle
                mm = int(drydep_slot_mm(lspec,lphase,m))

                if (lphase == 1) then
                   ptend%lq(mm) = .TRUE.
                   call outfld( trim(cnst_name(mm))//'DDV', drydep_diag_ddv(:,:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name(mm))//'DDF', drydep_diag_sflx(:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name(mm))//'TBF', drydep_diag_dep_trb(:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name(mm))//'GVF', drydep_diag_dep_grv(:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name(mm))//'DTQ', ptend%q(:,:,mm), pcols, lchnk )
                   aerdepdryis(:ncol,mm) = drydep_diag_sflx(:ncol,lspec,lphase,m)
                else
                   call outfld( trim(cnst_name_cw(mm))//'DDF', drydep_diag_sflx(:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name_cw(mm))//'TBF', drydep_diag_dep_trb(:,lspec,lphase,m), pcols, lchnk )
                   call outfld( trim(cnst_name_cw(mm))//'GVF', drydep_diag_dep_grv(:,lspec,lphase,m), pcols, lchnk )
                   aerdepdrycw(:ncol,mm) = drydep_diag_sflx(:ncol,lspec,lphase,m)
                end if
             end do
          end do
       end do

       if (apply_srf_drydep_local) then
          call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
       end if
       return
    end if



    do m = 1, ntot_amode   ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial / cloud-borne forms

          if (lphase == 1) then   ! interstial aerosol - calc settling/dep velocities of mode

! rad_aer = volume mean wet radius (m)
! dgncur_awet = geometric mean wet diameter for number distribution (m)
             rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,m)   &
                                 *exp(1.5_r8*(alnsg_amode(m)**2))
! dens_aer(1:ncol,:) = wet density (kg/m3)
             dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
             sg_aer(1:ncol,:) = sigmag_amode(m)

             jvlc = 1
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  & 
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, fraction_landuse(:,:,lchnk))
             jvlc = 2
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  & 
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, fraction_landuse(:,:,lchnk))
          end if

          do lspec = 0, nspec_amode(m)+1   ! loop over number + constituents + water

             if (lspec == 0) then   ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jvlc = 1
                else
                   mm = numptrcw_amode(m)
                   jvlc = 3
                endif
             else if (lspec <= nspec_amode(m)) then   ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jvlc = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jvlc = 4
                endif
             else   ! water mass
!   bypass dry deposition of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
!                  mm = lwaterptr_amode(m)
                   jvlc = 2
                else
                   mm = 0
                   jvlc = 4
                endif
             endif


          if (mm <= 0) cycle

!         if (lphase == 1) then
          if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
             ptend%lq(mm) = .TRUE.

             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             call outfld( trim(cnst_name(mm))//'DDV', pvmzaer(:,2:pverp), pcols, lchnk )

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     state%q(:,:,mm),  pvmzaer,  ptend%q(:,:,mm), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), state%q(:,:,mm), state%pmid, &
                               state%pdel, tvs, sflx, ptend%q(:,:,mm), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             call outfld( trim(cnst_name(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'GVF', dep_grv, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'DTQ', ptend%q(:,:,mm), pcols, lchnk)
             aerdepdryis(:ncol,mm) = sflx(:ncol)

          else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then  ! aerosol water
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     qaerwat(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), qaerwat(:,:,mm), state%pmid, &
                               state%pdel, tvs, sflx, dqdt_tmp(:,:), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             qaerwat(1:ncol,:,mm) = qaerwat(1:ncol,:,mm) + dqdt_tmp(1:ncol,:) * dt

          else  ! lphase == 2
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)
             fldcw => qqcw_get_field(pbuf, mm,lchnk)

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     fldcw(:,:),  pvmzaer,  dqdt_tmp(:,:), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), fldcw(:,:), state%pmid, &
                               state%pdel, tvs, sflx, dqdt_tmp(:,:), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

             call outfld( trim(cnst_name_cw(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name_cw(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name_cw(mm))//'GVF', dep_grv, pcols, lchnk )
             aerdepdrycw(:ncol,mm) = sflx(:ncol)

          endif

          enddo   ! lspec = 0, nspec_amode(m)+1
       enddo   ! lphase = 1, 2
    enddo   ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then 
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (apply_srf_drydep_local) then
       call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    endif

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_drydep_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_DRYDEP_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_drydep_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_drydep_use_native_impl = .false.
    end if

    aero_model_drydep_impl_selected = .true.

    if (masterproc) then
       if (aero_model_drydep_use_native_impl) then
          write(iulog,*) 'aero_model_drydep implementation = native'
       else
          write(iulog,*) 'aero_model_drydep implementation = codon'
          if (.not. aero_model_drydep_proof_written) then
             call aero_model_drydep_append_impl_proof('AERO_MODEL_DRYDEP_PROOF_FILE', &
                  'aero_model_drydep selector entered implementation = codon')
             aero_model_drydep_proof_written = .true.
          end if
       end if
       call flush(iulog)
    end if

  end subroutine aero_model_drydep_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep_select_branches(apply_srf_drydep_in)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    logical, intent(in) :: apply_srf_drydep_in

    integer(c_int64_t), target :: branch_mask_c

    interface
       subroutine aero_model_drydep_select_branches_codon(apply_srf_drydep_c, branch_mask_p) &
            bind(c, name="aero_model_drydep_select_branches_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: apply_srf_drydep_c
         type(c_ptr), value :: branch_mask_p
       end subroutine aero_model_drydep_select_branches_codon
    end interface

    if (aero_model_drydep_branch_selected) return

    branch_mask_c = 0_c_int64_t
    call aero_model_drydep_select_branches_codon( &
         merge(1_c_int64_t, 0_c_int64_t, apply_srf_drydep_in), &
         c_loc(branch_mask_c) &
    )

    aero_model_drydep_branch_mask = int(branch_mask_c)
    aero_model_drydep_branch_selected = .true.

  end subroutine aero_model_drydep_select_branches

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep_append_impl_proof(env_name, proof_line)

    character(len=*), intent(in) :: env_name, proof_line

    character(len=512) :: proof_path
    integer :: status, n, unit_id

    call get_environment_variable(env_name, value=proof_path, length=n, status=status)
    if (status /= 0 .or. n <= 0) return

    open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
         position='append', iostat=status)
    if (status /= 0) return

    write(unit_id,'(A)') trim(proof_line)
    close(unit_id)

  end subroutine aero_model_drydep_append_impl_proof

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep_fullshell_codon_wrap(ncol, dt, pint, pdel, rho, vlc_dry, vlc_trb, vlc_grv, &
                                                    state_q, ptend_q, qqcw_ptrs, slot_active, slot_mm, slot_jvlc, &
                                                    qqcw_mode_phase, diag_ddv, diag_dqdt, diag_sflx, diag_dep_trb, &
                                                    diag_dep_grv)

    use modal_aero_data
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: dt
    real(r8), target, intent(in) :: pint(pcols,pverp), pdel(pcols,pver), rho(pcols,pver)
    real(r8), target, intent(in) :: vlc_dry(pcols,pver,4,ntot_amode), vlc_trb(pcols,4,ntot_amode)
    real(r8), target, intent(in) :: vlc_grv(pcols,pver,4,ntot_amode)
    real(r8), target, intent(in) :: state_q(pcols,pver,pcnst)
    real(r8), target, intent(inout) :: ptend_q(pcols,pver,pcnst)
    type(c_ptr), intent(in) :: qqcw_ptrs(pcnst)
    integer(c_int64_t), target, intent(inout) :: slot_active(drydep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_mm(drydep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_jvlc(drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: qqcw_mode_phase(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_ddv(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_dqdt(pcols,pver,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx(pcols,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_dep_trb(pcols,drydep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_dep_grv(pcols,drydep_mode_phase_nslot,2,ntot_amode)

    integer :: m, lphase, lspec, slot, mm_local, jvlc_local
    real(r8), target :: q_work(pcols,pver), dqdt_work(pcols,pver), pvmzaer_work(pcols,pverp)
    real(r8), target :: fxdust(pcols,pverp), psi(pcols,pverp), fdot(pcols,pverp), xxk(pcols,pver)
    real(r8), target :: sflx_work(pcols), fxdot(pcols), fxdd(pcols), psistar(pcols), xins(pcols)
    real(r8), target :: s(pcols,pverp), sh(pcols,pverp), d(pcols,pverp), dh(pcols,pverp)
    real(r8), target :: e(pcols,pverp), eh(pcols,pverp), ppl(pcols,pverp), ppr(pcols,pverp)
    real(r8), target :: delxh(pcols,pverp)
    real(r8), parameter :: drydep_mxsedfac = 0.99_r8
    integer(c_int64_t), target :: intz(pcols), status_code, fail_i, fail_k
    real(r8), pointer :: fldcw(:,:)
    character(len=192) :: wrap_proof_line

    interface
       subroutine aero_model_drydep_fullshell_codon(ncol_c, pcols_c, pver_c, pverp_c, ntot_amode_c, nslot_max_c, &
            dt_c, gravit_c, mxsedfac_c, pint_p, pdel_p, rho_p, vlc_dry_p, vlc_trb_p, vlc_grv_p, state_q_p, &
            ptend_q_p, qqcw_mode_phase_p, slot_active_p, slot_mm_p, slot_jvlc_p, diag_ddv_p, diag_dqdt_p, &
            diag_sflx_p, diag_dep_trb_p, diag_dep_grv_p, q_work_p, dqdt_work_p, pvmzaer_work_p, sflx_work_p, &
            fxdust_p, psi_p, fdot_p, xxk_p, fxdot_p, fxdd_p, psistar_p, s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p, &
            xins_p, intz_p, status_p, fail_i_p, fail_k_p) bind(c, name="aero_model_drydep_fullshell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, ntot_amode_c, nslot_max_c
         real(c_double), value :: dt_c, gravit_c, mxsedfac_c
         type(c_ptr), value :: pint_p, pdel_p, rho_p, vlc_dry_p, vlc_trb_p, vlc_grv_p, state_q_p, ptend_q_p
         type(c_ptr), value :: qqcw_mode_phase_p, slot_active_p, slot_mm_p, slot_jvlc_p, diag_ddv_p, diag_dqdt_p
         type(c_ptr), value :: diag_sflx_p, diag_dep_trb_p, diag_dep_grv_p, q_work_p, dqdt_work_p, pvmzaer_work_p
         type(c_ptr), value :: sflx_work_p, fxdust_p, psi_p, fdot_p, xxk_p, fxdot_p, fxdd_p, psistar_p, s_p, sh_p, d_p
         type(c_ptr), value :: dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p, xins_p, intz_p, status_p, fail_i_p, fail_k_p
       end subroutine aero_model_drydep_fullshell_codon
    end interface

    slot_active(:,:,:) = 0_c_int64_t
    slot_mm(:,:,:) = 0_c_int64_t
    slot_jvlc(:,:,:) = 0_c_int64_t
    qqcw_mode_phase(:,:,:,:,:) = 0._r8
    diag_ddv(:,:,:,:,:) = 0._r8
    diag_dqdt(:,:,:,:,:) = 0._r8
    diag_sflx(:,:,:,:) = 0._r8
    diag_dep_trb(:,:,:,:) = 0._r8
    diag_dep_grv(:,:,:,:) = 0._r8

    do m = 1, ntot_amode
       do lphase = 1, 2
          do lspec = 0, nspec_amode(m)+1
             slot = lspec + 1
             if (lspec == 0) then
                if (lphase == 1) then
                   mm_local = numptr_amode(m)
                   jvlc_local = 1
                else
                   mm_local = numptrcw_amode(m)
                   jvlc_local = 3
                end if
             else if (lspec <= nspec_amode(m)) then
                if (lphase == 1) then
                   mm_local = lmassptr_amode(lspec,m)
                   jvlc_local = 2
                else
                   mm_local = lmassptrcw_amode(lspec,m)
                   jvlc_local = 4
                end if
             else
                cycle
             end if

             if (mm_local <= 0) cycle

             slot_active(slot,lphase,m) = 1_c_int64_t
             slot_mm(slot,lphase,m) = int(mm_local, c_int64_t)
             slot_jvlc(slot,lphase,m) = int(jvlc_local, c_int64_t)

             if (lphase == 2) then
                call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, mm_local, fldcw)
                qqcw_mode_phase(:ncol,:,slot,lphase,m) = fldcw(:ncol,:)
             end if
          end do
       end do
    end do

    if (masterproc .and. .not. aero_model_drydep_fullshell_wrap_proof_written) then
       wrap_proof_line = 'aero_model_drydep_fullshell_codon_wrap entered (parent shell = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_drydep_append_impl_proof('AERO_MODEL_DRYDEP_PROOF_FILE', trim(wrap_proof_line))
       aero_model_drydep_fullshell_wrap_proof_written = .true.
       call flush(iulog)
    end if

    status_code = 0_c_int64_t
    fail_i = 0_c_int64_t
    fail_k = 0_c_int64_t
    call aero_model_drydep_fullshell_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
         int(ntot_amode, c_int64_t), int(drydep_mode_phase_nslot, c_int64_t), real(dt, c_double), &
         real(gravit, c_double), real(drydep_mxsedfac, c_double), c_loc(pint), c_loc(pdel), c_loc(rho), c_loc(vlc_dry), &
         c_loc(vlc_trb), c_loc(vlc_grv), c_loc(state_q), c_loc(ptend_q), c_loc(qqcw_mode_phase), c_loc(slot_active), &
         c_loc(slot_mm), c_loc(slot_jvlc), c_loc(diag_ddv), c_loc(diag_dqdt), c_loc(diag_sflx), c_loc(diag_dep_trb), &
         c_loc(diag_dep_grv), c_loc(q_work), c_loc(dqdt_work), c_loc(pvmzaer_work), c_loc(sflx_work), c_loc(fxdust), c_loc(psi), &
         c_loc(fdot), c_loc(xxk), c_loc(fxdot), c_loc(fxdd), c_loc(psistar), c_loc(s), c_loc(sh), c_loc(d), &
         c_loc(dh), c_loc(e), c_loc(eh), c_loc(ppl), c_loc(ppr), c_loc(delxh), c_loc(xins), c_loc(intz), &
         c_loc(status_code), c_loc(fail_i), c_loc(fail_k) )

    if (status_code /= 0_c_int64_t) then
       write(iulog,*) 'aero_model_drydep_fullshell_codon_wrap -- interval was not found ', int(fail_i), int(fail_k)
       call endrun('aero_model_drydep_fullshell_codon_wrap -- interval was not found')
    end if

    do m = 1, ntot_amode
       do slot = 1, drydep_mode_phase_nslot
          if (slot_active(slot,2,m) == 0_c_int64_t) cycle
          call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, int(slot_mm(slot,2,m)), fldcw)
          fldcw(:ncol,:) = qqcw_mode_phase(:ncol,:,slot,2,m)
       end do
    end do

  end subroutine aero_model_drydep_fullshell_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_append_impl_proof(env_name, proof_line)

    character(len=*), intent(in) :: env_name, proof_line

    character(len=512) :: proof_path
    integer :: status, n, unit_id

    call get_environment_variable(env_name, value=proof_path, length=n, status=status)
    if (status /= 0 .or. n <= 0) return

    open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
         position='append', iostat=status)
    if (status /= 0) return

    write(unit_id,'(A)') trim(proof_line)
    close(unit_id)

  end subroutine aero_model_wetdep_append_impl_proof

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_wetdep_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_WETDEP_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_wetdep_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_wetdep_use_native_impl = .false.
    end if

    aero_model_wetdep_impl_selected = .true.

    if (masterproc) then
       if (aero_model_wetdep_use_native_impl) then
          write(iulog,*) 'aero_model_wetdep implementation = native'
       else
          write(iulog,*) 'aero_model_wetdep implementation = codon'
          if (.not. aero_model_wetdep_proof_written) then
             call aero_model_wetdep_append_impl_proof('AERO_MODEL_WETDEP_PROOF_FILE', &
                  'aero_model_wetdep selector entered implementation = codon')
             aero_model_wetdep_proof_written = .true.
          end if
       end if
       call flush(iulog)
    end if

  end subroutine aero_model_wetdep_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_select_mode_phase()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_wetdep_mode_phase_selected) return

    impl_name = 'cloudborne'
    call get_environment_variable('AERO_MODEL_WETDEP_MODE_PHASE_DIRECT', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       select case (trim(adjustl(impl_name(:n))))
       case ('both')
          aero_model_wetdep_mode_phase_use_interstitial = .true.
          aero_model_wetdep_mode_phase_use_cloudborne = .true.
       case ('interstitial')
          aero_model_wetdep_mode_phase_use_interstitial = .true.
          aero_model_wetdep_mode_phase_use_cloudborne = .false.
       case ('cloudborne')
          aero_model_wetdep_mode_phase_use_interstitial = .false.
          aero_model_wetdep_mode_phase_use_cloudborne = .true.
       case ('none', 'legacy')
          aero_model_wetdep_mode_phase_use_interstitial = .false.
          aero_model_wetdep_mode_phase_use_cloudborne = .false.
       case default
          aero_model_wetdep_mode_phase_use_interstitial = .true.
          aero_model_wetdep_mode_phase_use_cloudborne = .true.
       end select
    else
       aero_model_wetdep_mode_phase_use_interstitial = .false.
       aero_model_wetdep_mode_phase_use_cloudborne = .true.
    end if

    aero_model_wetdep_mode_phase_selected = .true.

    if (masterproc .and. .not. aero_model_wetdep_use_native_impl) then
       write(iulog,'(A,L1,A,L1)') 'aero_model_wetdep mode_phase direct interstitial=', &
            aero_model_wetdep_mode_phase_use_interstitial, ' cloudborne=', aero_model_wetdep_mode_phase_use_cloudborne
       call flush(iulog)
    end if

  end subroutine aero_model_wetdep_select_mode_phase

  !=============================================================================
  !=============================================================================
  subroutine qqcw2vmr_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (qqcw2vmr_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('QQCW2VMR_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       qqcw2vmr_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       qqcw2vmr_use_native_impl = .false.
    end if

    qqcw2vmr_impl_selected = .true.

    if (masterproc) then
       if (qqcw2vmr_use_native_impl) then
          write(iulog,*) 'qqcw2vmr implementation = native'
       else
          write(iulog,*) 'qqcw2vmr implementation = codon'
       end if
    end if

  end subroutine qqcw2vmr_select_impl

  !=============================================================================
  !=============================================================================
  subroutine vmr2qqcw_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (vmr2qqcw_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('VMR2QQCW_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       vmr2qqcw_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       vmr2qqcw_use_native_impl = .false.
    end if

    vmr2qqcw_impl_selected = .true.

    if (masterproc) then
       if (vmr2qqcw_use_native_impl) then
          write(iulog,*) 'vmr2qqcw implementation = native'
       else
          write(iulog,*) 'vmr2qqcw implementation = codon'
       end if
    end if

  end subroutine vmr2qqcw_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use modal_aero_deposition, only: set_srf_wetdep, set_srf_wetdep_codon_direct
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, wetdep_inputs_set_codon_direct, wetdep_inputs_t
    use modal_aero_data
    use modal_aero_calcsize,   only: modal_aero_calcsize_sub
    use modal_aero_wateruptake,only: modal_aero_wateruptake_dr
    use iso_c_binding,         only: c_ptr, c_int64_t


    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    integer :: m ! tracer index
    logical :: use_mode_phase_direct, use_fullshell_direct

    integer :: lchnk ! chunk identifier
    integer :: ncol ! number of atmospheric columns

    real(r8) :: iscavt(pcols, pver)

    integer :: mm
    integer :: i,k
    integer, parameter :: wetdep_stage_prepare_tracer = 1
    integer, parameter :: wetdep_stage_finish_interstitial = 2
    integer, parameter :: wetdep_stage_update_water = 3
    integer, parameter :: wetdep_stage_finish_cloudborne = 4

    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)

    real(r8) :: sflx(pcols) ! deposition flux
    real(r8) :: sflx_ics(pcols)
    real(r8) :: sflx_iss(pcols)
    real(r8) :: sflx_bcs(pcols)
    real(r8) :: sflx_bss(pcols)

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: lspec ! index for aerosol number / chem-mass / water-mass
    integer :: lcoardust, lcoarnacl ! indices for coarse mode dust and seasalt masses
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode ! rce 2010/05/02
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl ! rce 2010/05/02
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: hygro_sum_old(pcols,pver) ! before removal [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_del(pcols,pver) ! removal change to [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_old_ik, hygro_sum_new_ik
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    real(r8) :: water_old, water_new ! temporary old/new aerosol water mix-rat
    logical  :: isprx(pcols,pver) ! true if precipation
    real(r8) :: aerdepwetis(pcols,pcnst) ! aerosol wet deposition (interstitial)
    real(r8) :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)
    real(r8) :: codon_dummy2d_a(pcols,pver)
    real(r8) :: codon_dummy2d_b(pcols,pver)
    real(r8) :: codon_dummy2d_c(pcols,pver)
    real(r8) :: codon_dummy1d(pcols)
    real(r8), pointer :: fldcw(:,:)
    type(c_ptr) :: qqcw_ptrs(pcnst)
    integer(c_int64_t) :: wetdep_phase_active(2,ntot_amode)
    integer(c_int64_t) :: wetdep_slot_active(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t) :: wetdep_slot_mm(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t) :: wetdep_slot_jnv(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t) :: wetdep_slot_mass_kind(wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_slot_hygro_scale(wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_qqcw_mode_phase(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_dqdt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_icscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_isscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_bcscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_bsscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_sflx(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_sflx_ics(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_sflx_iss(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_sflx_bcs(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8) :: wetdep_diag_sflx_bss(pcols,wetdep_mode_phase_nslot,2,ntot_amode)

    real(r8), pointer :: dgnumwet(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    type(wetdep_inputs_t) :: dep_inputs

    lchnk = state%lchnk
    ncol  = state%ncol

    call aero_model_wetdep_select_impl()
    call aero_model_wetdep_select_mode_phase()

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)
    
    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled
    
    call t_startf('calcsize')
    ! for prognostic modal aerosols the transfer of mass between aitken and accumulation
    ! modes is done in conjunction with the dry radius calculation
    call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
    call t_stopf('calcsize')

    call t_startf('wateruptake')
    call modal_aero_wateruptake_dr(state, pbuf)
    call t_stopf('wateruptake')

    if (nwetdep<1) return

    if (aero_model_wetdep_use_native_impl) then
       call wetdep_inputs_set( state, pbuf, dep_inputs )
    else
       call wetdep_inputs_set_codon_direct( state, pbuf, dep_inputs, prec, isprx )
    end if

    call pbuf_get_field(pbuf, dgnumwet_idx,       dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, qaerwat_idx,        qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, fracis_idx,         fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
    if (.not. aero_model_wetdep_use_native_impl) then
       call qqcw_fill_cptrs(pbuf, qqcw_ptrs)
    end if

    if (aero_model_wetdep_use_native_impl) then
       prec(:ncol)=0._r8
       do k=1,pver
          where (prec(:ncol) >= 1.e-7_r8)
             isprx(:ncol,k) = .true.
          elsewhere
             isprx(:ncol,k) = .false.
          endwhere
          prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
               *state%pdel(:ncol,k)/gravit
       end do
    end if

    ! calculate the mass-weighted sol_factic for coarse mode species
    ! sol_factic_coarse(:,:) = 0.30_r8 ! tuned 1/4
    f_act_conv_coarse_dust = 0.40_r8 ! rce 2010/05/02
    f_act_conv_coarse_nacl = 0.80_r8 ! rce 2010/05/02
    f_act_conv_coarse(:,:) = 0.60_r8 ! rce 2010/05/02
    if (modeptr_coarse > 0) then
       lcoardust = lptr_dust_a_amode(modeptr_coarse)
       lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
       if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
          call aero_model_wetdep_fill_f_act_conv_coarse(ncol, dt, lcoardust, lcoarnacl, state%q, ptend%q, f_act_conv_coarse)
       end if
    end if

    scavcoefnv(:,:,0) = 0.0_r8 ! below-cloud scavcoef = 0.0 for cloud-borne species

    ! The interstitial phase updates qqcw, so the parent-shell batched path is
    ! only order-safe when both phases stay inside the same Codon call.
    use_fullshell_direct = .not. aero_model_wetdep_use_native_impl
    use_fullshell_direct = use_fullshell_direct .and. aero_model_wetdep_mode_phase_use_interstitial
    use_fullshell_direct = use_fullshell_direct .and. aero_model_wetdep_mode_phase_use_cloudborne

    wetdep_phase_active(:,:) = 0_c_int64_t
    if (use_fullshell_direct) then
       call aero_model_wetdep_fullshell_codon_wrap( &
            ncol, dt, state%pmid, state%q(:,:,1), state%pdel, dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
            dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, dep_inputs%evapr, &
            dep_inputs%totcond, dep_inputs%cldvcu, dep_inputs%cldvst, dlf, isprx, dgnumwet, state%q, ptend%q, &
            fracis, f_act_conv_coarse, qqcw_ptrs, wetdep_phase_active, wetdep_slot_active, wetdep_slot_mm, &
            wetdep_slot_jnv, wetdep_slot_mass_kind, wetdep_slot_hygro_scale, wetdep_qqcw_mode_phase, &
            wetdep_diag_dqdt, wetdep_diag_icscavt, wetdep_diag_isscavt, wetdep_diag_bcscavt, wetdep_diag_bsscavt, &
            wetdep_diag_sflx, wetdep_diag_sflx_ics, wetdep_diag_sflx_iss, wetdep_diag_sflx_bcs, &
            wetdep_diag_sflx_bss )
    end if

    do m = 1, ntot_amode ! main loop over aerosol modes

       do lphase = 1, 2 ! loop over interstitial (1) and cloud-borne (2) forms

          if (use_fullshell_direct) then
             use_mode_phase_direct = wetdep_phase_active(lphase,m) /= 0_c_int64_t
          else
             use_mode_phase_direct = .not. aero_model_wetdep_use_native_impl
             if (lphase == 1) use_mode_phase_direct = use_mode_phase_direct .and. aero_model_wetdep_mode_phase_use_interstitial
             if (lphase == 2) use_mode_phase_direct = use_mode_phase_direct .and. aero_model_wetdep_mode_phase_use_cloudborne
          end if

          if (use_mode_phase_direct) then
             if (.not. use_fullshell_direct) then
                call aero_model_wetdep_mode_phase_codon_wrap( &
                     m, lphase, ncol, dt, state%pmid, state%q(:,:,1), state%pdel, dep_inputs%cldt, dep_inputs%cldcu, &
                     dep_inputs%cmfdqr, dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, dep_inputs%cldvcu, dep_inputs%cldvst, dlf, isprx, dgnumwet, &
                     state%q, ptend%q, qaerwat(:,:,m), fracis, f_act_conv_coarse, qqcw_ptrs, wetdep_slot_active(:,lphase,m), &
                     wetdep_slot_mm(:,lphase,m), wetdep_slot_jnv(:,lphase,m), wetdep_slot_mass_kind(:,lphase,m), &
                     wetdep_slot_hygro_scale(:,lphase,m), wetdep_qqcw_mode_phase(:,:,:,lphase,m), q_tmp, iscavt, f_act_conv, &
                     sol_factic, hygro_sum_old, hygro_sum_del, scavcoefnv(:,:,1), scavcoefnv(:,:,2), &
                     wetdep_diag_dqdt(:,:,:,lphase,m), wetdep_diag_icscavt(:,:,:,lphase,m), wetdep_diag_isscavt(:,:,:,lphase,m), &
                     wetdep_diag_bcscavt(:,:,:,lphase,m), wetdep_diag_bsscavt(:,:,:,lphase,m), wetdep_diag_sflx(:,:,lphase,m), &
                     wetdep_diag_sflx_ics(:,:,lphase,m), wetdep_diag_sflx_iss(:,:,lphase,m), wetdep_diag_sflx_bcs(:,:,lphase,m), &
                     wetdep_diag_sflx_bss(:,:,lphase,m) )
             end if

             if (lphase == 1) then
                do lspec = 1, wetdep_mode_phase_nslot
                   if (wetdep_slot_active(lspec,lphase,m) == 0_c_int64_t) cycle
                   mm = int(wetdep_slot_mm(lspec,lphase,m))
                   ptend%lq(mm) = .TRUE.
                   call outfld( trim(cnst_name(mm))//'WET', wetdep_diag_dqdt(:,:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SIC', wetdep_diag_icscavt(:,:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SIS', wetdep_diag_isscavt(:,:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SBC', wetdep_diag_bcscavt(:,:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SBS', wetdep_diag_bsscavt(:,:,lspec,lphase,m), pcols, lchnk)
                   aerdepwetis(:ncol,mm) = wetdep_diag_sflx(:ncol,lspec,lphase,m)
                   call outfld( trim(cnst_name(mm))//'SFWET', wetdep_diag_sflx(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SFSIC', wetdep_diag_sflx_ics(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SFSIS', wetdep_diag_sflx_iss(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SFSBC', wetdep_diag_sflx_bcs(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name(mm))//'SFSBS', wetdep_diag_sflx_bss(:,lspec,lphase,m), pcols, lchnk)
                end do
             else
                do lspec = 1, wetdep_mode_phase_nslot
                   if (wetdep_slot_active(lspec,lphase,m) == 0_c_int64_t) cycle
                   mm = int(wetdep_slot_mm(lspec,lphase,m))
                   if (.not. use_fullshell_direct) then
                      call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, mm, fldcw)
                      fldcw(:ncol,:) = wetdep_qqcw_mode_phase(:ncol,:,lspec,lphase,m)
                   end if
                   aerdepwetcw(:ncol,mm) = wetdep_diag_sflx(:ncol,lspec,lphase,m)
                   call outfld( trim(cnst_name_cw(mm))//'SFWET', wetdep_diag_sflx(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name_cw(mm))//'SFSIC', wetdep_diag_sflx_ics(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name_cw(mm))//'SFSIS', wetdep_diag_sflx_iss(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name_cw(mm))//'SFSBC', wetdep_diag_sflx_bcs(:,lspec,lphase,m), pcols, lchnk)
                   call outfld( trim(cnst_name_cw(mm))//'SFSBS', wetdep_diag_sflx_bss(:,lspec,lphase,m), pcols, lchnk)
                end do
             end if

             cycle
          end if

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !
          if (lphase == 1) then ! interstial aerosol
             hygro_sum_old(:,:) = 0.0_r8
             hygro_sum_del(:,:) = 0.0_r8
             if (aero_model_wetdep_use_native_impl) then
                call modal_aero_bcscavcoef_get( m, ncol, isprx, dgnumwet, &
                     scavcoefnv(:,:,1), scavcoefnv(:,:,2) )
             else
                call aero_model_wetdep_bcscavcoef_codon_wrap( m, ncol, isprx, dgnumwet, &
                     scavcoefnv(:,:,1), scavcoefnv(:,:,2) )
             end if

             sol_factb = sol_factb_interstitial ! all below-cloud scav ON (0.1 "tuning factor")

             sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial

             sol_factic = sol_factic_interstitial

             if (m == modeptr_pcarbon) then
                ! sol_factic = 0.0_r8 ! conv in-cloud scav OFF (0.0 activation fraction)
                f_act_conv = 0.0_r8 ! rce 2010/05/02
             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
                ! sol_factic = 0.2_r8 ! conv in-cloud scav ON (0.5 activation fraction) ! tuned 1/4
                f_act_conv = 0.4_r8 ! rce 2010/05/02
             else
                ! sol_factic = 0.4_r8 ! conv in-cloud scav ON (1.0 activation fraction) ! tuned 1/4
                f_act_conv = 0.8_r8 ! rce 2010/05/02
             end if

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                                   !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          end if
          !
          ! rce 2010/05/03
          ! wetdepa has "sol_fact" parameters:
          ! sol_facti, sol_factic, sol_factb for liquid cloud

          do lspec = 0, nspec_amode(m)+1 ! loop over number + chem constituents + water

             if (lspec == 0) then ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jnv = 1
                else
                   mm = numptrcw_amode(m)
                   jnv = 0
                endif
             else if (lspec <= nspec_amode(m)) then ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jnv = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jnv = 0
                endif
             else ! water mass
                ! bypass wet removal of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
                   ! mm = lwaterptr_amode(m)
                   jnv = 2
                else
                   mm = 0
                   jnv = 0
                endif
             endif

             if (mm <= 0) cycle


             ! set f_act_conv for interstitial (lphase=1) coarse mode species
             ! for the convective in-cloud, we conceptually treat the coarse dust and seasalt
             ! as being externally mixed, and apply f_act_conv = f_act_conv_coarse_dust/nacl to dust/seasalt
             ! number and sulfate are conceptually partitioned to the dust and seasalt
             ! on a mass basis, so the f_act_conv for number and sulfate are
             ! mass-weighted averages of the values used for dust/seasalt
             if ((lphase == 1) .and. (m == modeptr_coarse)) then
                ! sol_factic = sol_factic_coarse
                f_act_conv = f_act_conv_coarse ! rce 2010/05/02
                if (lspec > 0) then
                   if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                      ! sol_factic = 0.2_r8 ! tuned 1/4
                      f_act_conv = f_act_conv_coarse_dust ! rce 2010/05/02
                   else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                      ! sol_factic = 0.4_r8 ! tuned 1/6
                      f_act_conv = f_act_conv_coarse_nacl ! rce 2010/05/02
                   end if
                end if
             end if


             if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
                ptend%lq(mm) = .TRUE.
                dqdt_tmp(:,:) = 0.0_r8

                if (aero_model_wetdep_use_native_impl) then
                   fldcw => qqcw_get_field(pbuf, mm,lchnk)
                   ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                   q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                else
                   call aero_model_wetdep_codon_wrap( &
                        wetdep_stage_prepare_tracer, ncol, dt, 0.0_r8, state%pdel, state%q(:,:,mm), ptend%q(:,:,mm), &
                        q_tmp, dqdt_tmp, sflx, sflx_ics, sflx_iss, sflx_bcs, sflx_bss, hygro_sum_old, hygro_sum_del, &
                        codon_dummy2d_a, codon_dummy2d_c, icscavt, isscavt, bcscavt, bsscavt, codon_dummy1d &
                   )
                   call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, mm, fldcw)
                end if

                if (aero_model_wetdep_use_native_impl) then
                   call wetdepa_v2( state%pmid, state%q(:,:,1), state%pdel, &
                        dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                        dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                        dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                        dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                        dlf, fracis(:,:,mm), sol_factb, ncol, &
                        scavcoefnv(:,:,jnv), &
                        is_strat_cloudborne=.false.,  &
                        qqcw=fldcw,  &
                        f_act_conv=f_act_conv, &
                        icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                        sol_facti_in=sol_facti, sol_factic_in=sol_factic )
                else
                   call aero_model_wetdep_scavenging_codon_wrap( &
                        2, ncol, dt, state%pmid, state%q(:,:,1), state%pdel, dep_inputs%cldt, dep_inputs%cldcu, &
                        dep_inputs%cmfdqr, dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                        dep_inputs%evapr, dep_inputs%totcond, q_tmp, dqdt_tmp, iscavt, dep_inputs%cldvcu, &
                        dep_inputs%cldvst, dlf, fracis(:,:,mm), sol_factb, scavcoefnv(:,:,jnv), fldcw, f_act_conv, &
                        icscavt, isscavt, bcscavt, bsscavt, sol_facti, sol_factic )
                end if

                if (aero_model_wetdep_use_native_impl) then
                   ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)
                else
                   tmpa = 0.0_r8
                   if (lspec > 0) then
                      tmpa = spechygro(lspectype_amode(lspec,m))/ &
                           specdens_amode(lspectype_amode(lspec,m))
                   end if
                   call aero_model_wetdep_codon_wrap( &
                        wetdep_stage_finish_interstitial, ncol, dt, tmpa, state%pdel, state%q(:,:,mm), ptend%q(:,:,mm), &
                        q_tmp, dqdt_tmp, sflx, sflx_ics, sflx_iss, sflx_bcs, sflx_bss, hygro_sum_old, hygro_sum_del, &
                        codon_dummy2d_a, codon_dummy2d_c, icscavt, isscavt, bcscavt, bsscavt, aerdepwetis(:,mm) &
                   )
                end if

                call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

                if (aero_model_wetdep_use_native_impl) then
                   call aero_model_wetdep_column_flux(ncol, dqdt_tmp, state%pdel, sflx)
                   aerdepwetis(:ncol,mm) = sflx(:ncol)

                   call aero_model_wetdep_column_flux(ncol, icscavt, state%pdel, sflx)
                   sflx_ics(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, isscavt, state%pdel, sflx)
                   sflx_iss(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, bcscavt, state%pdel, sflx)
                   sflx_bcs(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, bsscavt, state%pdel, sflx)
                   sflx_bss(:) = sflx(:)

                   if (lspec > 0) then
                      tmpa = spechygro(lspectype_amode(lspec,m))/ &
                           specdens_amode(lspectype_amode(lspec,m))
                      tmpb = tmpa*dt
                      hygro_sum_old(1:ncol,:) = hygro_sum_old(1:ncol,:) &
                           + tmpa*q_tmp(1:ncol,:)
                      hygro_sum_del(1:ncol,:) = hygro_sum_del(1:ncol,:) &
                           + tmpb*dqdt_tmp(1:ncol,:)
                   end if
                end if
                call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SFSIC', sflx_ics, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SFSIS', sflx_iss, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SFSBC', sflx_bcs, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SFSBS', sflx_bss, pcols, lchnk)

             else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then
                ! aerosol water -- because of how wetdepa treats evaporation of stratiform
                ! precip, it is not appropriate to apply wetdepa to aerosol water
                ! instead, "hygro_sum" = [sum of (mass*hygro/dens)] is calculated before and
                ! after wet removal, and new water is calculated using
                ! new_water = old_water*min(10,(hygro_sum_new/hygro_sum_old))
                ! the "min(10,...)" is to avoid potential problems when hygro_sum_old ~= 0
                ! also, individual wet removal terms (ic,is,bc,bs) are not output to history
                ! ptend%lq(mm) = .TRUE.
                ! dqdt_tmp(:,:) = 0.0_r8
                if (aero_model_wetdep_use_native_impl) then
                   do k = 1, pver
                      do i = 1, ncol
                         ! water_old = max( 0.0_r8, state%q(i,k,mm)+ptend%q(i,k,mm)*dt )
                         water_old = max( 0.0_r8, qaerwat(i,k,mm) )
                         hygro_sum_old_ik = max( 0.0_r8, hygro_sum_old(i,k) )
                         hygro_sum_new_ik = max( 0.0_r8, hygro_sum_old_ik+hygro_sum_del(i,k) )
                         if (hygro_sum_new_ik >= 10.0_r8*hygro_sum_old_ik) then
                            water_new = 10.0_r8*water_old
                         else
                            water_new = water_old*(hygro_sum_new_ik/hygro_sum_old_ik)
                         end if
                         ! dqdt_tmp(i,k) = (water_new - water_old)/dt
                         qaerwat(i,k,mm) = water_new
                      end do
                   end do
                else
                   call aero_model_wetdep_codon_wrap( &
                        wetdep_stage_update_water, ncol, dt, 0.0_r8, state%pdel, codon_dummy2d_a, codon_dummy2d_b, &
                        q_tmp, dqdt_tmp, sflx, sflx_ics, sflx_iss, sflx_bcs, sflx_bss, hygro_sum_old, hygro_sum_del, &
                        qaerwat(:,:,mm), codon_dummy2d_c, icscavt, isscavt, bcscavt, bsscavt, codon_dummy1d &
                   )
                end if

                ! ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                ! call outfld( trim(cnst_name(mm))

                ! sflx(:)=0._r8
                ! do k=1,pver
                ! do i=1,ncol
                ! sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                ! enddo
                ! enddo
                ! call outfld( trim(cnst_name(mm))

             else ! lphase == 2
                dqdt_tmp(:,:) = 0.0_r8

                if (aero_model_wetdep_use_native_impl) then
                   fldcw => qqcw_get_field(pbuf, mm,lchnk)
                   call wetdepa_v2(state%pmid, state%q(:,:,1), state%pdel, &
                        dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                        dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                        dep_inputs%evapr, dep_inputs%totcond, fldcw, dt, &
                        dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                        dlf, fracis_cw, sol_factb, ncol, &
                        scavcoefnv(:,:,jnv), &
                        is_strat_cloudborne=.true.,  &
                        icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                        sol_facti_in=sol_facti, sol_factic_in=sol_factic )
                else
                   call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, mm, fldcw)
                   call aero_model_wetdep_scavenging_codon_wrap( &
                        1, ncol, dt, state%pmid, state%q(:,:,1), state%pdel, dep_inputs%cldt, dep_inputs%cldcu, &
                        dep_inputs%cmfdqr, dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                        dep_inputs%evapr, dep_inputs%totcond, fldcw, dqdt_tmp, iscavt, dep_inputs%cldvcu, &
                        dep_inputs%cldvst, dlf, fracis_cw, sol_factb, scavcoefnv(:,:,jnv), codon_dummy2d_a, f_act_conv, &
                        icscavt, isscavt, bcscavt, bsscavt, sol_facti, sol_factic )
                end if

                if (aero_model_wetdep_use_native_impl) then
                   fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt
                else
                   call aero_model_wetdep_codon_wrap( &
                        wetdep_stage_finish_cloudborne, ncol, dt, 0.0_r8, state%pdel, codon_dummy2d_a, codon_dummy2d_b, &
                        q_tmp, dqdt_tmp, sflx, sflx_ics, sflx_iss, sflx_bcs, sflx_bss, hygro_sum_old, hygro_sum_del, &
                        codon_dummy2d_c, fldcw, icscavt, isscavt, bcscavt, bsscavt, aerdepwetcw(:,mm) &
                   )
                end if

                if (aero_model_wetdep_use_native_impl) then
                   call aero_model_wetdep_column_flux(ncol, dqdt_tmp, state%pdel, sflx)
                   aerdepwetcw(:ncol,mm) = sflx(:ncol)

                   call aero_model_wetdep_column_flux(ncol, icscavt, state%pdel, sflx)
                   sflx_ics(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, isscavt, state%pdel, sflx)
                   sflx_iss(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, bcscavt, state%pdel, sflx)
                   sflx_bcs(:) = sflx(:)
                   call aero_model_wetdep_column_flux(ncol, bsscavt, state%pdel, sflx)
                   sflx_bss(:) = sflx(:)
                end if
                call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
                call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx_ics, pcols, lchnk)
                call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx_iss, pcols, lchnk)
                call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx_bcs, pcols, lchnk)
                call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx_bss, pcols, lchnk)

             endif

          enddo ! lspec = 0, nspec_amode(m)+1
       enddo ! lphase = 1, 2
    enddo ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       if (aero_model_wetdep_use_native_impl) then
          call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
       else
          call set_srf_wetdep_codon_direct(aerdepwetis, aerdepwetcw, cam_out)
       end if
    endif

  endsubroutine aero_model_wetdep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, qqcw_index, fldcw)

    use iso_c_binding, only: c_ptr, c_associated, c_f_pointer
    use cam_abortutils, only: endrun

    type(c_ptr), intent(in) :: qqcw_ptrs(pcnst)
    integer, intent(in) :: qqcw_index
    real(r8), pointer, intent(out) :: fldcw(:,:)

    nullify(fldcw)

    if (qqcw_index < 1 .or. qqcw_index > pcnst) then
       call endrun('aero_model_wetdep_resolve_qqcw_ptr: qqcw index out of range')
    end if
    if (.not. c_associated(qqcw_ptrs(qqcw_index))) then
       call endrun('aero_model_wetdep_resolve_qqcw_ptr: unresolved qqcw pointer')
    end if

    call c_f_pointer(qqcw_ptrs(qqcw_index), fldcw, (/pcols, pver/))

  end subroutine aero_model_wetdep_resolve_qqcw_ptr

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_fullshell_codon_wrap(ncol, dt, pmid, q1, pdel, cldt, cldcu, cmfdqr, &
                                                    evapc, conicw, prain, qme, evapr, totcond, cldvcu, cldvst, dlf, &
                                                    isprx, dgnumwet, state_q, ptend_q, fracis_full, f_act_conv_coarse, &
                                                    qqcw_ptrs, phase_active, slot_active, slot_mm, slot_jnv, &
                                                    slot_mass_kind, slot_hygro_scale, qqcw_mode_phase, diag_dqdt, &
                                                    diag_icscavt, diag_isscavt, diag_bcscavt, diag_bsscavt, diag_sflx, &
                                                    diag_sflx_ics, diag_sflx_iss, diag_sflx_bcs, diag_sflx_bss)

    use modal_aero_data
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: dt
    real(r8), target, intent(in) :: pmid(pcols,pver), q1(pcols,pver), pdel(pcols,pver), cldt(pcols,pver), cldcu(pcols,pver)
    real(r8), target, intent(in) :: cmfdqr(pcols,pver), evapc(pcols,pver), conicw(pcols,pver), prain(pcols,pver)
    real(r8), target, intent(in) :: qme(pcols,pver), evapr(pcols,pver), totcond(pcols,pver), cldvcu(pcols,pver)
    real(r8), target, intent(in) :: cldvst(pcols,pver), dlf(pcols,pver), dgnumwet(pcols,pver,ntot_amode)
    logical, intent(in) :: isprx(pcols,pver)
    real(r8), target, intent(in) :: state_q(pcols,pver,pcnst), f_act_conv_coarse(pcols,pver)
    real(r8), target, intent(inout) :: ptend_q(pcols,pver,pcnst), fracis_full(pcols,pver,pcnst)
    type(c_ptr), intent(in) :: qqcw_ptrs(pcnst)
    integer(c_int64_t), target, intent(inout) :: phase_active(2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_active(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_mm(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_jnv(wetdep_mode_phase_nslot,2,ntot_amode)
    integer(c_int64_t), target, intent(inout) :: slot_mass_kind(wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: slot_hygro_scale(wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: qqcw_mode_phase(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_dqdt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_icscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_isscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_bcscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_bsscavt(pcols,pver,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx_ics(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx_iss(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx_bcs(pcols,wetdep_mode_phase_nslot,2,ntot_amode)
    real(r8), target, intent(inout) :: diag_sflx_bss(pcols,wetdep_mode_phase_nslot,2,ntot_amode)

    integer :: i, k, m, lphase, lspec, slot, mm_local, jnv_local
    logical :: use_mode_phase_direct
    real(r8) :: omsm
    real(r8), target :: phase_sol_factb(2,ntot_amode), phase_sol_facti(2,ntot_amode)
    real(r8), target :: phase_sol_factic_scalar(2,ntot_amode), phase_base_f_act_scalar(2,ntot_amode)
    real(r8), target :: phase_dgnum_mode(ntot_amode)
    real(r8), target :: scavimptblnum_all(nimptblgrow_mind:nimptblgrow_maxd,ntot_amode)
    real(r8), target :: scavimptblvol_all(nimptblgrow_mind:nimptblgrow_maxd,ntot_amode)
    integer(c_int64_t), target :: phase_is_coarse_interstitial(2,ntot_amode), isprx_mask(pcols,pver)
    real(r8), target :: q_tmp_work(pcols,pver), hygro_sum_old(pcols,pver), hygro_sum_del(pcols,pver)
    real(r8), target :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver), scavcoef_work(pcols,pver)
    real(r8), target :: iscavt_work(pcols,pver), f_act_conv_work(pcols,pver), sol_factic_work(pcols,pver)
    real(r8), target :: fracis_dummy_work(pcols,pver), wetdep_dblchek_hist(pcols,pver), wetdep_srct_hist(pcols,pver)
    real(r8), target :: wetdep_rat_hist(pcols,pver), wetdep_fracev_hist(pcols,pver)
    real(r8), target :: wetdep_clds(pcols), wetdep_fracev(pcols), wetdep_fracev_cu(pcols), wetdep_fracp(pcols)
    real(r8), target :: wetdep_pdog(pcols), wetdep_rpdog(pcols), wetdep_precabc(pcols), wetdep_precabs(pcols)
    real(r8), target :: wetdep_rat(pcols), wetdep_scavab(pcols), wetdep_scavabc(pcols), wetdep_srcc(pcols)
    real(r8), target :: wetdep_srcs(pcols), wetdep_srct(pcols), wetdep_fins(pcols), wetdep_finc(pcols)
    real(r8), target :: wetdep_conv_scav_ic(pcols), wetdep_conv_scav_bc(pcols), wetdep_st_scav_ic(pcols)
    real(r8), target :: wetdep_st_scav_bc(pcols), wetdep_odds(pcols), wetdep_dblchek(pcols), wetdep_trac_qqcw(pcols)
    real(r8), target :: wetdep_tracer_incu(pcols), wetdep_tracer_mean(pcols)
    real(r8), pointer :: fldcw(:,:)
    character(len=192) :: wrap_proof_line

    interface
       subroutine aero_model_wetdep_fullshell_codon(ncol_c, pcols_c, pver_c, ntot_amode_c, nslot_max_c, &
            nimptblgrow_mind_c, nimptblgrow_maxd_c, dt_c, gravit_c, omsm_c, dlndg_nimptblgrow_c, &
            f_act_conv_coarse_dust_c, f_act_conv_coarse_nacl_c, pmid_p, q1_p, pdel_p, cldt_p, cldcu_p, cmfdqr_p, &
            evapc_p, conicw_p, prain_p, qme_p, evapr_p, totcond_p, cldvcu_p, cldvst_p, dlf_p, phase_active_p, &
            phase_sol_factb_p, phase_sol_facti_p, phase_sol_factic_scalar_p, phase_base_f_act_scalar_p, &
            phase_dgnum_mode_p, phase_is_coarse_interstitial_p, isprx_mask_p, dgnumwet_p, scavimptblnum_all_p, &
            scavimptblvol_all_p, state_q_p, ptend_q_p, fracis_full_p, f_act_conv_coarse_p, qqcw_mode_phase_p, &
            q_tmp_work_p, hygro_sum_old_p, hygro_sum_del_p, scavcoefnum_p, scavcoefvol_p, scavcoef_work_p, &
            iscavt_work_p, f_act_conv_work_p, sol_factic_work_p, slot_active_p, slot_mm_p, slot_jnv_p, &
            slot_mass_kind_p, slot_hygro_scale_p, diag_dqdt_p, diag_icscavt_p, diag_isscavt_p, diag_bcscavt_p, &
            diag_bsscavt_p, diag_sflx_p, diag_sflx_ics_p, diag_sflx_iss_p, diag_sflx_bcs_p, diag_sflx_bss_p, &
            clds_p, fracev_p, fracev_cu_p, fracp_p, pdog_p, rpdog_p, precabc_p, precabs_p, rat_p, scavab_p, &
            scavabc_p, srcc_p, srcs_p, srct_p, fins_p, finc_p, conv_scav_ic_p, conv_scav_bc_p, st_scav_ic_p, &
            st_scav_bc_p, odds_p, dblchek_p, trac_qqcw_p, tracer_incu_p, tracer_mean_p, fracis_dummy_p, &
            dblchek_hist_p, srct_hist_p, rat_hist_p, fracev_hist_p) bind(c, name="aero_model_wetdep_fullshell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ntot_amode_c, nslot_max_c
         integer(c_int64_t), value :: nimptblgrow_mind_c, nimptblgrow_maxd_c
         real(c_double), value :: dt_c, gravit_c, omsm_c, dlndg_nimptblgrow_c
         real(c_double), value :: f_act_conv_coarse_dust_c, f_act_conv_coarse_nacl_c
         type(c_ptr), value :: pmid_p, q1_p, pdel_p, cldt_p, cldcu_p, cmfdqr_p, evapc_p, conicw_p, prain_p, qme_p
         type(c_ptr), value :: evapr_p, totcond_p, cldvcu_p, cldvst_p, dlf_p, phase_active_p, phase_sol_factb_p
         type(c_ptr), value :: phase_sol_facti_p, phase_sol_factic_scalar_p, phase_base_f_act_scalar_p
         type(c_ptr), value :: phase_dgnum_mode_p, phase_is_coarse_interstitial_p, isprx_mask_p, dgnumwet_p
         type(c_ptr), value :: scavimptblnum_all_p, scavimptblvol_all_p, state_q_p, ptend_q_p, fracis_full_p
         type(c_ptr), value :: f_act_conv_coarse_p, qqcw_mode_phase_p, q_tmp_work_p, hygro_sum_old_p
         type(c_ptr), value :: hygro_sum_del_p, scavcoefnum_p, scavcoefvol_p, scavcoef_work_p, iscavt_work_p
         type(c_ptr), value :: f_act_conv_work_p, sol_factic_work_p, slot_active_p, slot_mm_p, slot_jnv_p
         type(c_ptr), value :: slot_mass_kind_p, slot_hygro_scale_p, diag_dqdt_p, diag_icscavt_p, diag_isscavt_p
         type(c_ptr), value :: diag_bcscavt_p, diag_bsscavt_p, diag_sflx_p, diag_sflx_ics_p, diag_sflx_iss_p
         type(c_ptr), value :: diag_sflx_bcs_p, diag_sflx_bss_p, clds_p, fracev_p, fracev_cu_p, fracp_p, pdog_p
         type(c_ptr), value :: rpdog_p, precabc_p, precabs_p, rat_p, scavab_p, scavabc_p, srcc_p, srcs_p, srct_p
         type(c_ptr), value :: fins_p, finc_p, conv_scav_ic_p, conv_scav_bc_p, st_scav_ic_p, st_scav_bc_p, odds_p
         type(c_ptr), value :: dblchek_p, trac_qqcw_p, tracer_incu_p, tracer_mean_p, fracis_dummy_p
         type(c_ptr), value :: dblchek_hist_p, srct_hist_p, rat_hist_p, fracev_hist_p
       end subroutine aero_model_wetdep_fullshell_codon
    end interface

    phase_active(:,:) = 0_c_int64_t
    phase_sol_factb(:,:) = 0.0_r8
    phase_sol_facti(:,:) = 0.0_r8
    phase_sol_factic_scalar(:,:) = 0.0_r8
    phase_base_f_act_scalar(:,:) = 0.0_r8
    phase_dgnum_mode(:) = 0.0_r8
    scavimptblnum_all(:,:) = scavimptblnum(:,:)
    scavimptblvol_all(:,:) = scavimptblvol(:,:)
    phase_is_coarse_interstitial(:,:) = 0_c_int64_t
    slot_active(:,:,:) = 0_c_int64_t
    slot_mm(:,:,:) = 0_c_int64_t
    slot_jnv(:,:,:) = 0_c_int64_t
    slot_mass_kind(:,:,:) = 0_c_int64_t
    slot_hygro_scale(:,:,:) = 0.0_r8
    qqcw_mode_phase(:,:,:,:,:) = 0.0_r8
    diag_dqdt(:,:,:,:,:) = 0.0_r8
    diag_icscavt(:,:,:,:,:) = 0.0_r8
    diag_isscavt(:,:,:,:,:) = 0.0_r8
    diag_bcscavt(:,:,:,:,:) = 0.0_r8
    diag_bsscavt(:,:,:,:,:) = 0.0_r8
    diag_sflx(:,:,:,:) = 0.0_r8
    diag_sflx_ics(:,:,:,:) = 0.0_r8
    diag_sflx_iss(:,:,:,:) = 0.0_r8
    diag_sflx_bcs(:,:,:,:) = 0.0_r8
    diag_sflx_bss(:,:,:,:) = 0.0_r8

    do k = 1, pver
       do i = 1, ncol
          isprx_mask(i,k) = merge(1_c_int64_t, 0_c_int64_t, isprx(i,k))
       end do
    end do

    do m = 1, ntot_amode
       phase_dgnum_mode(m) = dgnum_amode(m)
       do lphase = 1, 2
          use_mode_phase_direct = .false.
          if (lphase == 1) use_mode_phase_direct = aero_model_wetdep_mode_phase_use_interstitial
          if (lphase == 2) use_mode_phase_direct = aero_model_wetdep_mode_phase_use_cloudborne
          if (.not. use_mode_phase_direct) cycle

          phase_active(lphase,m) = 1_c_int64_t
          if (lphase == 1) then
             phase_sol_factb(lphase,m) = sol_factb_interstitial
             phase_sol_facti(lphase,m) = 0.0_r8
             phase_sol_factic_scalar(lphase,m) = sol_factic_interstitial
             if (m == modeptr_pcarbon) then
                phase_base_f_act_scalar(lphase,m) = 0.0_r8
             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
                phase_base_f_act_scalar(lphase,m) = 0.4_r8
             else
                phase_base_f_act_scalar(lphase,m) = 0.8_r8
             end if
             phase_is_coarse_interstitial(lphase,m) = merge(1_c_int64_t, 0_c_int64_t, m == modeptr_coarse)
          else
             phase_sol_factb(lphase,m) = 0.0_r8
             phase_sol_facti(lphase,m) = sol_facti_cloud_borne
             phase_sol_factic_scalar(lphase,m) = 0.0_r8
             phase_base_f_act_scalar(lphase,m) = 0.0_r8
             phase_is_coarse_interstitial(lphase,m) = 0_c_int64_t
          end if

          do lspec = 0, nspec_amode(m)+1
             slot = lspec + 1
             if (lspec == 0) then
                if (lphase == 1) then
                   mm_local = numptr_amode(m)
                   jnv_local = 1
                else
                   mm_local = numptrcw_amode(m)
                   jnv_local = 0
                end if
             else if (lspec <= nspec_amode(m)) then
                if (lphase == 1) then
                   mm_local = lmassptr_amode(lspec,m)
                   jnv_local = 2
                else
                   mm_local = lmassptrcw_amode(lspec,m)
                   jnv_local = 0
                end if
             else
                cycle
             end if

             if (mm_local <= 0) cycle

             slot_active(slot,lphase,m) = 1_c_int64_t
             slot_mm(slot,lphase,m) = int(mm_local, c_int64_t)
             slot_jnv(slot,lphase,m) = int(jnv_local, c_int64_t)

             if ((lphase == 1) .and. (lspec > 0)) then
                slot_hygro_scale(slot,lphase,m) = spechygro(lspectype_amode(lspec,m)) / &
                     specdens_amode(lspectype_amode(lspec,m))
                if (m == modeptr_coarse) then
                   if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                      slot_mass_kind(slot,lphase,m) = 1_c_int64_t
                   else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                      slot_mass_kind(slot,lphase,m) = 2_c_int64_t
                   end if
                end if
             end if

             call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, mm_local, fldcw)
             qqcw_mode_phase(:ncol,:,slot,lphase,m) = fldcw(:ncol,:)
          end do
       end do
    end do

    if (masterproc .and. .not. aero_model_wetdep_fullshell_wrap_proof_written) then
       wrap_proof_line = 'aero_model_wetdep_fullshell_codon_wrap entered (mode_phase parent shell = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_wetdep_append_impl_proof('AERO_MODEL_WETDEP_PROOF_FILE', trim(wrap_proof_line))
       aero_model_wetdep_fullshell_wrap_proof_written = .true.
       call flush(iulog)
    end if

    omsm = 1._r8 - 2*epsilon(1._r8)
    call aero_model_wetdep_fullshell_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(ntot_amode, c_int64_t), &
         int(wetdep_mode_phase_nslot, c_int64_t), int(nimptblgrow_mind, c_int64_t), int(nimptblgrow_maxd, c_int64_t), &
         real(dt, c_double), real(gravit, c_double), real(omsm, c_double), real(dlndg_nimptblgrow, c_double), &
         real(0.40_r8, c_double), real(0.80_r8, c_double), c_loc(pmid), c_loc(q1), c_loc(pdel), c_loc(cldt), &
         c_loc(cldcu), c_loc(cmfdqr), c_loc(evapc), c_loc(conicw), c_loc(prain), c_loc(qme), c_loc(evapr), &
         c_loc(totcond), c_loc(cldvcu), c_loc(cldvst), c_loc(dlf), c_loc(phase_active), c_loc(phase_sol_factb), &
         c_loc(phase_sol_facti), c_loc(phase_sol_factic_scalar), c_loc(phase_base_f_act_scalar), c_loc(phase_dgnum_mode), &
         c_loc(phase_is_coarse_interstitial), c_loc(isprx_mask), c_loc(dgnumwet), c_loc(scavimptblnum_all(nimptblgrow_mind,1)), &
         c_loc(scavimptblvol_all(nimptblgrow_mind,1)), c_loc(state_q), c_loc(ptend_q), c_loc(fracis_full), &
         c_loc(f_act_conv_coarse), c_loc(qqcw_mode_phase), c_loc(q_tmp_work), c_loc(hygro_sum_old), c_loc(hygro_sum_del), &
         c_loc(scavcoefnum), c_loc(scavcoefvol), c_loc(scavcoef_work), c_loc(iscavt_work), c_loc(f_act_conv_work), &
         c_loc(sol_factic_work), c_loc(slot_active), c_loc(slot_mm), c_loc(slot_jnv), c_loc(slot_mass_kind), &
         c_loc(slot_hygro_scale), c_loc(diag_dqdt), c_loc(diag_icscavt), c_loc(diag_isscavt), c_loc(diag_bcscavt), &
         c_loc(diag_bsscavt), c_loc(diag_sflx), c_loc(diag_sflx_ics), c_loc(diag_sflx_iss), c_loc(diag_sflx_bcs), &
         c_loc(diag_sflx_bss), c_loc(wetdep_clds), c_loc(wetdep_fracev), c_loc(wetdep_fracev_cu), c_loc(wetdep_fracp), &
         c_loc(wetdep_pdog), c_loc(wetdep_rpdog), c_loc(wetdep_precabc), c_loc(wetdep_precabs), c_loc(wetdep_rat), &
         c_loc(wetdep_scavab), c_loc(wetdep_scavabc), c_loc(wetdep_srcc), c_loc(wetdep_srcs), c_loc(wetdep_srct), &
         c_loc(wetdep_fins), c_loc(wetdep_finc), c_loc(wetdep_conv_scav_ic), c_loc(wetdep_conv_scav_bc), &
         c_loc(wetdep_st_scav_ic), c_loc(wetdep_st_scav_bc), c_loc(wetdep_odds), c_loc(wetdep_dblchek), &
         c_loc(wetdep_trac_qqcw), c_loc(wetdep_tracer_incu), c_loc(wetdep_tracer_mean), c_loc(fracis_dummy_work), &
         c_loc(wetdep_dblchek_hist), c_loc(wetdep_srct_hist), c_loc(wetdep_rat_hist), c_loc(wetdep_fracev_hist) &
    )

    do m = 1, ntot_amode
       if (phase_active(2,m) == 0_c_int64_t) cycle
       do slot = 1, wetdep_mode_phase_nslot
          if (slot_active(slot,2,m) == 0_c_int64_t) cycle
          call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, int(slot_mm(slot,2,m)), fldcw)
          fldcw(:ncol,:) = qqcw_mode_phase(:ncol,:,slot,2,m)
       end do
    end do

  end subroutine aero_model_wetdep_fullshell_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_mode_phase_codon_wrap(m, lphase, ncol, dt, pmid, q1, pdel, cldt, cldcu, cmfdqr, &
                                                     evapc, conicw, prain, qme, evapr, totcond, cldvcu, cldvst, dlf, &
                                                     isprx, dgnumwet, state_q, ptend_q, qaerwat_mode, fracis_full, &
                                                     f_act_conv_coarse, qqcw_ptrs, slot_active, slot_mm, slot_jnv, &
                                                     slot_mass_kind, slot_hygro_scale, qqcw_mode_phase, q_tmp_work, &
                                                     iscavt_work, f_act_conv_work, sol_factic_work, hygro_sum_old, &
                                                     hygro_sum_del, scavcoefnum, scavcoefvol, diag_dqdt, diag_icscavt, &
                                                     diag_isscavt, diag_bcscavt, diag_bsscavt, diag_sflx, diag_sflx_ics, &
                                                     diag_sflx_iss, diag_sflx_bcs, diag_sflx_bss)

    use modal_aero_data
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: m, lphase, ncol
    real(r8), intent(in) :: dt
    real(r8), target, intent(in) :: pmid(pcols,pver), q1(pcols,pver), pdel(pcols,pver), cldt(pcols,pver), cldcu(pcols,pver)
    real(r8), target, intent(in) :: cmfdqr(pcols,pver), evapc(pcols,pver), conicw(pcols,pver), prain(pcols,pver)
    real(r8), target, intent(in) :: qme(pcols,pver), evapr(pcols,pver), totcond(pcols,pver), cldvcu(pcols,pver)
    real(r8), target, intent(in) :: cldvst(pcols,pver), dlf(pcols,pver), dgnumwet(pcols,pver,ntot_amode)
    logical, intent(in) :: isprx(pcols,pver)
    real(r8), target, intent(in) :: state_q(pcols,pver,pcnst), f_act_conv_coarse(pcols,pver)
    real(r8), target, intent(inout) :: ptend_q(pcols,pver,pcnst), qaerwat_mode(pcols,pver), fracis_full(pcols,pver,pcnst)
    type(c_ptr), intent(in) :: qqcw_ptrs(pcnst)
    integer(c_int64_t), target, intent(inout) :: slot_active(wetdep_mode_phase_nslot), slot_mm(wetdep_mode_phase_nslot)
    integer(c_int64_t), target, intent(inout) :: slot_jnv(wetdep_mode_phase_nslot), slot_mass_kind(wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: slot_hygro_scale(wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: qqcw_mode_phase(pcols,pver,wetdep_mode_phase_nslot), q_tmp_work(pcols,pver)
    real(r8), target, intent(inout) :: iscavt_work(pcols,pver), f_act_conv_work(pcols,pver), sol_factic_work(pcols,pver)
    real(r8), target, intent(inout) :: hygro_sum_old(pcols,pver), hygro_sum_del(pcols,pver)
    real(r8), target, intent(inout) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)
    real(r8), target, intent(inout) :: diag_dqdt(pcols,pver,wetdep_mode_phase_nslot), diag_icscavt(pcols,pver,wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: diag_isscavt(pcols,pver,wetdep_mode_phase_nslot), diag_bcscavt(pcols,pver,wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: diag_bsscavt(pcols,pver,wetdep_mode_phase_nslot), diag_sflx(pcols,wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: diag_sflx_ics(pcols,wetdep_mode_phase_nslot), diag_sflx_iss(pcols,wetdep_mode_phase_nslot)
    real(r8), target, intent(inout) :: diag_sflx_bcs(pcols,wetdep_mode_phase_nslot), diag_sflx_bss(pcols,wetdep_mode_phase_nslot)

    integer :: i, k, lspec, slot, mm_local, jnv_local
    real(r8) :: dgnum_mode, sol_factb_local, sol_facti_local, sol_factic_scalar, base_f_act_scalar, omsm
    integer(c_int64_t) :: is_coarse_interstitial
    integer(c_int64_t), target :: isprx_mask(pcols,pver)
    real(r8), target :: scavimptblnum_mode(nimptblgrow_mind:nimptblgrow_maxd)
    real(r8), target :: scavimptblvol_mode(nimptblgrow_mind:nimptblgrow_maxd)
    real(r8), target :: scavcoef_work(pcols,pver)
    real(r8), target :: fracis_dummy_work(pcols,pver)
    real(r8), target :: wetdep_dblchek_hist(pcols,pver), wetdep_srct_hist(pcols,pver)
    real(r8), target :: wetdep_rat_hist(pcols,pver), wetdep_fracev_hist(pcols,pver)
    real(r8), target :: wetdep_clds(pcols), wetdep_fracev(pcols), wetdep_fracev_cu(pcols), wetdep_fracp(pcols)
    real(r8), target :: wetdep_pdog(pcols), wetdep_rpdog(pcols), wetdep_precabc(pcols), wetdep_precabs(pcols)
    real(r8), target :: wetdep_rat(pcols), wetdep_scavab(pcols), wetdep_scavabc(pcols), wetdep_srcc(pcols)
    real(r8), target :: wetdep_srcs(pcols), wetdep_srct(pcols), wetdep_fins(pcols), wetdep_finc(pcols)
    real(r8), target :: wetdep_conv_scav_ic(pcols), wetdep_conv_scav_bc(pcols), wetdep_st_scav_ic(pcols)
    real(r8), target :: wetdep_st_scav_bc(pcols), wetdep_odds(pcols), wetdep_dblchek(pcols), wetdep_trac_qqcw(pcols)
    real(r8), target :: wetdep_tracer_incu(pcols), wetdep_tracer_mean(pcols)
    real(r8), pointer :: fldcw(:,:)
    character(len=192) :: wrap_proof_line

    interface
       subroutine aero_model_wetdep_mode_phase_stage_dispatch_codon(m_c, lphase_c, ncol_c, pcols_c, pver_c, pcnst_c, ntot_amode_c, &
            nslot_max_c, nimptblgrow_mind_c, nimptblgrow_maxd_c, dt_c, gravit_c, omsm_c, dgnum_mode_c, &
            dlndg_nimptblgrow_c, sol_factb_c, sol_facti_c, sol_factic_scalar_c, base_f_act_scalar_c, &
            is_coarse_interstitial_c, f_act_conv_coarse_dust_c, f_act_conv_coarse_nacl_c, pmid_p, q1_p, pdel_p, &
            cldt_p, cldcu_p, cmfdqr_p, evapc_p, conicw_p, prain_p, qme_p, evapr_p, totcond_p, cldvcu_p, cldvst_p, &
            dlf_p, isprx_mask_p, dgnumwet_p, scavimptblnum_mode_p, scavimptblvol_mode_p, state_q_p, ptend_q_p, &
            qaerwat_mode_p, fracis_full_p, f_act_conv_coarse_p, qqcw_mode_phase_p, q_tmp_work_p, hygro_sum_old_p, &
            hygro_sum_del_p, scavcoefnum_p, scavcoefvol_p, scavcoef_work_p, iscavt_work_p, f_act_conv_work_p, &
            sol_factic_work_p, slot_active_p, slot_mm_p, slot_jnv_p, slot_mass_kind_p, slot_hygro_scale_p, &
            diag_dqdt_p, diag_icscavt_p, diag_isscavt_p, diag_bcscavt_p, diag_bsscavt_p, diag_sflx_p, &
            diag_sflx_ics_p, diag_sflx_iss_p, diag_sflx_bcs_p, diag_sflx_bss_p, clds_p, fracev_p, fracev_cu_p, &
            fracp_p, pdog_p, rpdog_p, precabc_p, precabs_p, rat_p, scavab_p, scavabc_p, srcc_p, srcs_p, srct_p, &
            fins_p, finc_p, conv_scav_ic_p, conv_scav_bc_p, st_scav_ic_p, st_scav_bc_p, odds_p, dblchek_p, &
            trac_qqcw_p, tracer_incu_p, tracer_mean_p, fracis_dummy_p, dblchek_hist_p, srct_hist_p, rat_hist_p, &
            fracev_hist_p) bind(c, name="aero_model_wetdep_mode_phase_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: m_c, lphase_c, ncol_c, pcols_c, pver_c, pcnst_c, ntot_amode_c, nslot_max_c
         integer(c_int64_t), value :: nimptblgrow_mind_c, nimptblgrow_maxd_c, is_coarse_interstitial_c
         real(c_double), value :: dt_c, gravit_c, omsm_c, dgnum_mode_c, dlndg_nimptblgrow_c, sol_factb_c, sol_facti_c
         real(c_double), value :: sol_factic_scalar_c, base_f_act_scalar_c, f_act_conv_coarse_dust_c
         real(c_double), value :: f_act_conv_coarse_nacl_c
         type(c_ptr), value :: pmid_p, q1_p, pdel_p, cldt_p, cldcu_p, cmfdqr_p, evapc_p, conicw_p, prain_p, qme_p
         type(c_ptr), value :: evapr_p, totcond_p, cldvcu_p, cldvst_p, dlf_p, isprx_mask_p, dgnumwet_p
         type(c_ptr), value :: scavimptblnum_mode_p, scavimptblvol_mode_p, state_q_p, ptend_q_p, qaerwat_mode_p
         type(c_ptr), value :: fracis_full_p, f_act_conv_coarse_p, qqcw_mode_phase_p, q_tmp_work_p, hygro_sum_old_p
         type(c_ptr), value :: hygro_sum_del_p, scavcoefnum_p, scavcoefvol_p, scavcoef_work_p, iscavt_work_p
         type(c_ptr), value :: f_act_conv_work_p, sol_factic_work_p, slot_active_p, slot_mm_p, slot_jnv_p
         type(c_ptr), value :: slot_mass_kind_p, slot_hygro_scale_p, diag_dqdt_p, diag_icscavt_p, diag_isscavt_p
         type(c_ptr), value :: diag_bcscavt_p, diag_bsscavt_p, diag_sflx_p, diag_sflx_ics_p, diag_sflx_iss_p
         type(c_ptr), value :: diag_sflx_bcs_p, diag_sflx_bss_p, clds_p, fracev_p, fracev_cu_p, fracp_p, pdog_p
         type(c_ptr), value :: rpdog_p, precabc_p, precabs_p, rat_p, scavab_p, scavabc_p, srcc_p, srcs_p, srct_p
         type(c_ptr), value :: fins_p, finc_p, conv_scav_ic_p, conv_scav_bc_p, st_scav_ic_p, st_scav_bc_p, odds_p
         type(c_ptr), value :: dblchek_p, trac_qqcw_p, tracer_incu_p, tracer_mean_p, fracis_dummy_p
         type(c_ptr), value :: dblchek_hist_p, srct_hist_p, rat_hist_p, fracev_hist_p
       end subroutine aero_model_wetdep_mode_phase_stage_dispatch_codon
    end interface

    slot_active(:) = 0_c_int64_t
    slot_mm(:) = 0_c_int64_t
    slot_jnv(:) = 0_c_int64_t
    slot_mass_kind(:) = 0_c_int64_t
    slot_hygro_scale(:) = 0.0_r8

    if (lphase == 1) then
       sol_factb_local = sol_factb_interstitial
       sol_facti_local = 0.0_r8
       sol_factic_scalar = sol_factic_interstitial
       if (m == modeptr_pcarbon) then
          base_f_act_scalar = 0.0_r8
       else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
          base_f_act_scalar = 0.4_r8
       else
          base_f_act_scalar = 0.8_r8
       end if
       is_coarse_interstitial = merge(1_c_int64_t, 0_c_int64_t, m == modeptr_coarse)
    else
       sol_factb_local = 0.0_r8
       sol_facti_local = sol_facti_cloud_borne
       sol_factic_scalar = 0.0_r8
       base_f_act_scalar = 0.0_r8
       is_coarse_interstitial = 0_c_int64_t
    end if

    do lspec = 0, nspec_amode(m)+1
       slot = lspec + 1
       if (lspec == 0) then
          if (lphase == 1) then
             mm_local = numptr_amode(m)
             jnv_local = 1
          else
             mm_local = numptrcw_amode(m)
             jnv_local = 0
          end if
       else if (lspec <= nspec_amode(m)) then
          if (lphase == 1) then
             mm_local = lmassptr_amode(lspec,m)
             jnv_local = 2
          else
             mm_local = lmassptrcw_amode(lspec,m)
             jnv_local = 0
          end if
       else
          cycle
       end if

       if (mm_local <= 0) cycle

       slot_active(slot) = 1_c_int64_t
       slot_mm(slot) = int(mm_local, c_int64_t)
       slot_jnv(slot) = int(jnv_local, c_int64_t)

       if ((lphase == 1) .and. (lspec > 0)) then
          slot_hygro_scale(slot) = spechygro(lspectype_amode(lspec,m)) / specdens_amode(lspectype_amode(lspec,m))
          if (m == modeptr_coarse) then
             if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                slot_mass_kind(slot) = 1_c_int64_t
             else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                slot_mass_kind(slot) = 2_c_int64_t
             end if
          end if
       end if
    end do

    isprx_mask(:,:) = 0_c_int64_t
    do k = 1, pver
       do i = 1, ncol
          isprx_mask(i,k) = merge(1_c_int64_t, 0_c_int64_t, isprx(i,k))
       end do
    end do

    dgnum_mode = dgnum_amode(m)
    scavimptblnum_mode(:) = scavimptblnum(:,m)
    scavimptblvol_mode(:) = scavimptblvol(:,m)
    omsm = 1._r8 - 2*epsilon(1._r8)

    if (masterproc .and. .not. aero_model_wetdep_mode_phase_wrap_proof_written) then
       wrap_proof_line = 'aero_model_wetdep_mode_phase_codon_wrap entered (unified wetdep inputs/cldst/clddiag/qqcw/scavcoef/mode phase stage dispatch = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_wetdep_append_impl_proof('AERO_MODEL_WETDEP_PROOF_FILE', trim(wrap_proof_line))
       aero_model_wetdep_mode_phase_wrap_proof_written = .true.
       call flush(iulog)
    end if

    qqcw_mode_phase(:,:,:) = 0.0_r8
    do slot = 1, wetdep_mode_phase_nslot
       if (slot_active(slot) == 0_c_int64_t) cycle
       call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, int(slot_mm(slot)), fldcw)
       qqcw_mode_phase(:ncol,:,slot) = fldcw(:ncol,:)
    end do

    call aero_model_wetdep_mode_phase_stage_dispatch_codon( &
         int(m, c_int64_t), int(lphase, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(pcnst, c_int64_t), int(ntot_amode, c_int64_t), int(wetdep_mode_phase_nslot, c_int64_t), &
         int(nimptblgrow_mind, c_int64_t), int(nimptblgrow_maxd, c_int64_t), real(dt, c_double), real(gravit, c_double), &
         real(omsm, c_double), real(dgnum_mode, c_double), real(dlndg_nimptblgrow, c_double), real(sol_factb_local, c_double), &
         real(sol_facti_local, c_double), real(sol_factic_scalar, c_double), real(base_f_act_scalar, c_double), &
         is_coarse_interstitial, real(0.40_r8, c_double), real(0.80_r8, c_double), c_loc(pmid), c_loc(q1), c_loc(pdel), &
         c_loc(cldt), c_loc(cldcu), c_loc(cmfdqr), c_loc(evapc), c_loc(conicw), c_loc(prain), c_loc(qme), c_loc(evapr), &
         c_loc(totcond), c_loc(cldvcu), c_loc(cldvst), c_loc(dlf), c_loc(isprx_mask), c_loc(dgnumwet), &
         c_loc(scavimptblnum_mode(nimptblgrow_mind)), c_loc(scavimptblvol_mode(nimptblgrow_mind)), c_loc(state_q), &
         c_loc(ptend_q), c_loc(qaerwat_mode), c_loc(fracis_full), c_loc(f_act_conv_coarse), c_loc(qqcw_mode_phase), &
         c_loc(q_tmp_work), c_loc(hygro_sum_old), c_loc(hygro_sum_del), c_loc(scavcoefnum), c_loc(scavcoefvol), &
         c_loc(scavcoef_work), c_loc(iscavt_work), c_loc(f_act_conv_work), c_loc(sol_factic_work), c_loc(slot_active), &
         c_loc(slot_mm), c_loc(slot_jnv), c_loc(slot_mass_kind), c_loc(slot_hygro_scale), c_loc(diag_dqdt), &
         c_loc(diag_icscavt), c_loc(diag_isscavt), c_loc(diag_bcscavt), c_loc(diag_bsscavt), c_loc(diag_sflx), &
         c_loc(diag_sflx_ics), c_loc(diag_sflx_iss), c_loc(diag_sflx_bcs), c_loc(diag_sflx_bss), c_loc(wetdep_clds), &
         c_loc(wetdep_fracev), c_loc(wetdep_fracev_cu), c_loc(wetdep_fracp), c_loc(wetdep_pdog), c_loc(wetdep_rpdog), &
         c_loc(wetdep_precabc), c_loc(wetdep_precabs), c_loc(wetdep_rat), c_loc(wetdep_scavab), c_loc(wetdep_scavabc), &
         c_loc(wetdep_srcc), c_loc(wetdep_srcs), c_loc(wetdep_srct), c_loc(wetdep_fins), c_loc(wetdep_finc), &
         c_loc(wetdep_conv_scav_ic), c_loc(wetdep_conv_scav_bc), c_loc(wetdep_st_scav_ic), c_loc(wetdep_st_scav_bc), &
         c_loc(wetdep_odds), c_loc(wetdep_dblchek), c_loc(wetdep_trac_qqcw), c_loc(wetdep_tracer_incu), &
         c_loc(wetdep_tracer_mean), c_loc(fracis_dummy_work), c_loc(wetdep_dblchek_hist), c_loc(wetdep_srct_hist), &
         c_loc(wetdep_rat_hist), c_loc(wetdep_fracev_hist) &
    )

    ! The active interstitial path cycles before the aerosol-water slot, so only
    ! cloud-borne qqcw is written back here.
    if (lphase /= 2) return

    do slot = 1, wetdep_mode_phase_nslot
       if (slot_active(slot) == 0_c_int64_t) cycle
       call aero_model_wetdep_resolve_qqcw_ptr(qqcw_ptrs, int(slot_mm(slot)), fldcw)
       fldcw(:ncol,:) = qqcw_mode_phase(:ncol,:,slot)
    end do

  end subroutine aero_model_wetdep_mode_phase_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_bcscavcoef_codon_wrap(m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol)

    use modal_aero_data, only: dgnum_amode
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: m, ncol
    logical, intent(in) :: isprx(pcols,pver)
    real(r8), target, intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), target, intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

    integer :: i, k
    real(r8) :: dgnum_mode
    integer(c_int64_t), target :: isprx_mask(pcols,pver)
    real(r8), target :: scavimptblnum_mode(nimptblgrow_mind:nimptblgrow_maxd)
    real(r8), target :: scavimptblvol_mode(nimptblgrow_mind:nimptblgrow_maxd)

    interface
       subroutine modal_aero_bcscavcoef_get_codon(m_c, ncol_c, pcols_c, pver_c, ntot_amode_c, nimptblgrow_mind_c, &
            nimptblgrow_maxd_c, dlndg_nimptblgrow_c, dgnum_mode_c, isprx_mask_p, dgn_awet_p, scavimptblnum_mode_p, &
            scavimptblvol_mode_p, scavcoefnum_p, scavcoefvol_p) bind(c, name="modal_aero_bcscavcoef_get_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: m_c, ncol_c, pcols_c, pver_c, ntot_amode_c, nimptblgrow_mind_c
         integer(c_int64_t), value :: nimptblgrow_maxd_c
         real(c_double), value :: dlndg_nimptblgrow_c, dgnum_mode_c
         type(c_ptr), value :: isprx_mask_p, dgn_awet_p, scavimptblnum_mode_p, scavimptblvol_mode_p
         type(c_ptr), value :: scavcoefnum_p, scavcoefvol_p
       end subroutine modal_aero_bcscavcoef_get_codon
    end interface

    do k = 1, pver
       do i = 1, ncol
          isprx_mask(i,k) = merge(1_c_int64_t, 0_c_int64_t, isprx(i,k))
       end do
    end do

    dgnum_mode = dgnum_amode(m)
    scavimptblnum_mode(:) = scavimptblnum(:,m)
    scavimptblvol_mode(:) = scavimptblvol(:,m)

    call modal_aero_bcscavcoef_get_codon( &
         int(m, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(ntot_amode, c_int64_t), int(nimptblgrow_mind, c_int64_t), int(nimptblgrow_maxd, c_int64_t), &
         real(dlndg_nimptblgrow, c_double), real(dgnum_mode, c_double), c_loc(isprx_mask), c_loc(dgn_awet), &
         c_loc(scavimptblnum_mode(nimptblgrow_mind)), c_loc(scavimptblvol_mode(nimptblgrow_mind)), &
         c_loc(scavcoefnum), c_loc(scavcoefvol) &
    )

  end subroutine aero_model_wetdep_bcscavcoef_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_codon_wrap(stage, ncol, dt, tmpa, pdel, state_tracer, ptend_tracer, q_tmp, dqdt_tmp, &
                                          sflx, sflx_ics, sflx_iss, sflx_bcs, sflx_bss, hygro_sum_old, hygro_sum_del, &
                                          qaerwat, fldcw, icscavt, isscavt, bcscavt, bsscavt, aerdep)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: stage, ncol
    real(r8), intent(in) :: dt, tmpa
    real(r8), target, intent(in) :: pdel(pcols,pver), state_tracer(pcols,pver), dqdt_tmp(pcols,pver)
    real(r8), target, intent(in) :: icscavt(pcols,pver), isscavt(pcols,pver), bcscavt(pcols,pver), bsscavt(pcols,pver)
    real(r8), target, intent(inout) :: ptend_tracer(pcols,pver), q_tmp(pcols,pver), hygro_sum_old(pcols,pver), &
                                       hygro_sum_del(pcols,pver), qaerwat(pcols,pver), fldcw(pcols,pver)
    real(r8), target, intent(inout) :: sflx(pcols), sflx_ics(pcols), sflx_iss(pcols), sflx_bcs(pcols), sflx_bss(pcols)
    real(r8), target, intent(inout) :: aerdep(pcols)
    character(len=192) :: wrap_proof_line

    interface
       subroutine aero_model_wetdep_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, dt_c, tmpa_c, gravit_c, pdel_p, &
                                          state_tracer_p, ptend_tracer_p, q_tmp_p, dqdt_p, sflx_p, sflx_ics_p, &
                                          sflx_iss_p, sflx_bcs_p, sflx_bss_p, hygro_sum_old_p, hygro_sum_del_p, &
                                          qaerwat_p, fldcw_p, icscavt_p, isscavt_p, bcscavt_p, bsscavt_p, aerdep_p) &
            bind(c, name="aero_model_wetdep_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: dt_c, tmpa_c, gravit_c
         type(c_ptr), value :: pdel_p, state_tracer_p, ptend_tracer_p, q_tmp_p, dqdt_p
         type(c_ptr), value :: sflx_p, sflx_ics_p, sflx_iss_p, sflx_bcs_p, sflx_bss_p
         type(c_ptr), value :: hygro_sum_old_p, hygro_sum_del_p, qaerwat_p, fldcw_p
         type(c_ptr), value :: icscavt_p, isscavt_p, bcscavt_p, bsscavt_p, aerdep_p
       end subroutine aero_model_wetdep_stage_dispatch_codon
    end interface

    if (masterproc .and. .not. aero_model_wetdep_wrap_proof_written) then
       wrap_proof_line = 'aero_model_wetdep_codon_wrap entered (unified wetdep shell stage dispatch = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_wetdep_append_impl_proof('AERO_MODEL_WETDEP_PROOF_FILE', trim(wrap_proof_line))
       aero_model_wetdep_wrap_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_wetdep_stage_dispatch_codon( &
         int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         real(dt, c_double), real(tmpa, c_double), real(gravit, c_double), c_loc(pdel), c_loc(state_tracer), &
         c_loc(ptend_tracer), c_loc(q_tmp), c_loc(dqdt_tmp), c_loc(sflx), c_loc(sflx_ics), c_loc(sflx_iss), &
         c_loc(sflx_bcs), c_loc(sflx_bss), c_loc(hygro_sum_old), c_loc(hygro_sum_del), c_loc(qaerwat), c_loc(fldcw), &
         c_loc(icscavt), c_loc(isscavt), c_loc(bcscavt), c_loc(bsscavt), c_loc(aerdep) &
    )

  end subroutine aero_model_wetdep_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_scavenging_codon_wrap(branch_mode, ncol, deltat, p, q, pdel, cldt, cldc, cmfdqr, &
                                                     evapc, conicw, precs, conds, evaps, cwat, tracer, scavt, iscavt, &
                                                     cldvcu, cldvst, dlf, fracis, sol_factb, scavcoef, qqcw, f_act_conv, &
                                                     icscavt, isscavt, bcscavt, bsscavt, sol_facti, sol_factic)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: branch_mode, ncol
    real(r8), intent(in) :: deltat, sol_factb, sol_facti
    real(r8), target, intent(in) :: p(pcols,pver), q(pcols,pver), pdel(pcols,pver), cldt(pcols,pver), cldc(pcols,pver)
    real(r8), target, intent(in) :: cmfdqr(pcols,pver), evapc(pcols,pver), conicw(pcols,pver), precs(pcols,pver)
    real(r8), target, intent(in) :: conds(pcols,pver), evaps(pcols,pver), cwat(pcols,pver), cldvcu(pcols,pver)
    real(r8), target, intent(in) :: cldvst(pcols,pver), dlf(pcols,pver), scavcoef(pcols,pver), qqcw(pcols,pver)
    real(r8), target, intent(in) :: f_act_conv(pcols,pver), sol_factic(pcols,pver)
    real(r8), target, intent(inout) :: tracer(pcols,pver)
    real(r8), target, intent(out) :: scavt(pcols,pver), iscavt(pcols,pver), fracis(pcols,pver)
    real(r8), target, intent(out) :: icscavt(pcols,pver), isscavt(pcols,pver), bcscavt(pcols,pver), bsscavt(pcols,pver)

    integer :: i, k
    logical :: found
    real(r8) :: omsm
    real(r8), target :: clds(pcols), fracev(pcols), fracev_cu(pcols), fracp(pcols), pdog(pcols), rpdog(pcols)
    real(r8), target :: precabc(pcols), precabs(pcols), rat(pcols), scavab(pcols), scavabc(pcols), srcc(pcols)
    real(r8), target :: srcs(pcols), srct(pcols), fins(pcols), finc(pcols), conv_scav_ic(pcols), conv_scav_bc(pcols)
    real(r8), target :: st_scav_ic(pcols), st_scav_bc(pcols), odds(pcols), dblchek(pcols), trac_qqcw(pcols)
    real(r8), target :: tracer_incu(pcols), tracer_mean(pcols), dblchek_hist(pcols,pver), srct_hist(pcols,pver)
    real(r8), target :: rat_hist(pcols,pver), fracev_hist(pcols,pver)

    interface
       subroutine wetdepa_v2_codon(pcols_c, pver_c, ncol_c, branch_mode_c, gravit_c, deltat_c, omsm_c, &
            sol_facti_c, sol_factb_c, p_p, q_p, pdel_p, cldt_p, cldc_p, cmfdqr_p, evapc_p, conicw_p, &
            precs_p, conds_p, evaps_p, cwat_p, tracer_p, scavt_p, iscavt_p, cldvcu_p, cldvst_p, dlf_p, &
            fracis_p, scavcoef_p, sol_factic_p, qqcw_p, f_act_conv_p, icscavt_p, isscavt_p, bcscavt_p, &
            bsscavt_p, clds_p, fracev_p, fracev_cu_p, fracp_p, pdog_p, rpdog_p, precabc_p, precabs_p, rat_p, &
            scavab_p, scavabc_p, srcc_p, srcs_p, srct_p, fins_p, finc_p, conv_scav_ic_p, conv_scav_bc_p, &
            st_scav_ic_p, st_scav_bc_p, odds_p, dblchek_p, trac_qqcw_p, tracer_incu_p, tracer_mean_p, &
            dblchek_hist_p, srct_hist_p, rat_hist_p, fracev_hist_p) bind(c, name="wetdepa_v2_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, pver_c, ncol_c, branch_mode_c
         real(c_double), value :: gravit_c, deltat_c, omsm_c, sol_facti_c, sol_factb_c
         type(c_ptr), value :: p_p, q_p, pdel_p, cldt_p, cldc_p, cmfdqr_p, evapc_p, conicw_p, precs_p
         type(c_ptr), value :: conds_p, evaps_p, cwat_p, tracer_p, scavt_p, iscavt_p, cldvcu_p, cldvst_p
         type(c_ptr), value :: dlf_p, fracis_p, scavcoef_p, sol_factic_p, qqcw_p, f_act_conv_p
         type(c_ptr), value :: icscavt_p, isscavt_p, bcscavt_p, bsscavt_p
         type(c_ptr), value :: clds_p, fracev_p, fracev_cu_p, fracp_p, pdog_p, rpdog_p, precabc_p, precabs_p
         type(c_ptr), value :: rat_p, scavab_p, scavabc_p, srcc_p, srcs_p, srct_p, fins_p, finc_p
         type(c_ptr), value :: conv_scav_ic_p, conv_scav_bc_p, st_scav_ic_p, st_scav_bc_p, odds_p, dblchek_p
         type(c_ptr), value :: trac_qqcw_p, tracer_incu_p, tracer_mean_p
         type(c_ptr), value :: dblchek_hist_p, srct_hist_p, rat_hist_p, fracev_hist_p
       end subroutine wetdepa_v2_codon
    end interface

    omsm = 1._r8 - 2*epsilon(1._r8)

    call wetdepa_v2_codon( &
         int(pcols, c_int64_t), int(pver, c_int64_t), int(ncol, c_int64_t), int(branch_mode, c_int64_t), &
         real(gravit, c_double), real(deltat, c_double), real(omsm, c_double), real(sol_facti, c_double), &
         real(sol_factb, c_double), c_loc(p), c_loc(q), c_loc(pdel), c_loc(cldt), c_loc(cldc), c_loc(cmfdqr), &
         c_loc(evapc), c_loc(conicw), c_loc(precs), c_loc(conds), c_loc(evaps), c_loc(cwat), c_loc(tracer), &
         c_loc(scavt), c_loc(iscavt), c_loc(cldvcu), c_loc(cldvst), c_loc(dlf), c_loc(fracis), c_loc(scavcoef), &
         c_loc(sol_factic), c_loc(qqcw), c_loc(f_act_conv), c_loc(icscavt), c_loc(isscavt), c_loc(bcscavt), &
         c_loc(bsscavt), c_loc(clds), c_loc(fracev), c_loc(fracev_cu), c_loc(fracp), c_loc(pdog), c_loc(rpdog), &
         c_loc(precabc), c_loc(precabs), c_loc(rat), c_loc(scavab), c_loc(scavabc), c_loc(srcc), c_loc(srcs), &
         c_loc(srct), c_loc(fins), c_loc(finc), c_loc(conv_scav_ic), c_loc(conv_scav_bc), c_loc(st_scav_ic), &
         c_loc(st_scav_bc), c_loc(odds), c_loc(dblchek), c_loc(trac_qqcw), c_loc(tracer_incu), c_loc(tracer_mean), &
         c_loc(dblchek_hist), c_loc(srct_hist), c_loc(rat_hist), c_loc(fracev_hist) &
    )

    do k = 1, pver
       found = .false.
       do i = 1, ncol
          if (dblchek_hist(i,k) < 0._r8) then
             found = .true.
             exit
          end if
       end do

       if (found) then
          do i = 1, ncol
             if (dblchek_hist(i,k) < 0._r8) then
                write(iulog,*) ' wetdapa: negative value ', i, k, tracer(i,k), &
                     dblchek_hist(i,k), scavt(i,k), srct_hist(i,k), rat_hist(i,k), fracev_hist(i,k)
             end if
          end do
       end if
    end do

  end subroutine aero_model_wetdep_scavenging_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_fill_f_act_conv_coarse(ncol, dt, lcoardust, lcoarnacl, state_q, ptend_q, f_act_conv_coarse)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol, lcoardust, lcoarnacl
    real(r8), intent(in) :: dt
    real(r8), target, intent(in) :: state_q(pcols,pver,pcnst), ptend_q(pcols,pver,pcnst)
    real(r8), target, intent(out) :: f_act_conv_coarse(pcols,pver)

    integer :: i, k
    real(r8) :: tmpdust, tmpnacl

    interface
       subroutine aero_model_wetdep_f_act_conv_coarse_codon(ncol_c, pcols_c, pver_c, pcnst_c, dt_c, &
            lcoardust_c, lcoarnacl_c, state_q_p, ptend_q_p, f_act_conv_coarse_p) &
            bind(c, name="aero_model_wetdep_f_act_conv_coarse_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, lcoardust_c, lcoarnacl_c
         real(c_double), value :: dt_c
         type(c_ptr), value :: state_q_p, ptend_q_p, f_act_conv_coarse_p
       end subroutine aero_model_wetdep_f_act_conv_coarse_codon
    end interface

    if (aero_model_wetdep_use_native_impl) then
       f_act_conv_coarse(:,:) = 0.60_r8
       do k = 1, pver
          do i = 1, ncol
             tmpdust = max( 0.0_r8, state_q(i,k,lcoardust) + ptend_q(i,k,lcoardust)*dt )
             tmpnacl = max( 0.0_r8, state_q(i,k,lcoarnacl) + ptend_q(i,k,lcoarnacl)*dt )
             if ((tmpdust+tmpnacl) > 1.0e-30_r8) then
                f_act_conv_coarse(i,k) = (0.40_r8*tmpdust + 0.80_r8*tmpnacl)/(tmpdust+tmpnacl)
             end if
          end do
       end do
       return
    end if

    call aero_model_wetdep_f_act_conv_coarse_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), real(dt, c_double), &
         int(lcoardust, c_int64_t), int(lcoarnacl, c_int64_t), c_loc(state_q), c_loc(ptend_q), c_loc(f_act_conv_coarse) &
    )

  end subroutine aero_model_wetdep_fill_f_act_conv_coarse

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep_column_flux(ncol, field, pdel, sflx)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: field(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(out) :: sflx(pcols)

    integer :: i, k

    interface
       subroutine aero_model_wetdep_column_flux_codon(ncol_c, pcols_c, pver_c, gravit_c, field_p, pdel_p, sflx_p) &
            bind(c, name="aero_model_wetdep_column_flux_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: gravit_c
         type(c_ptr), value :: field_p, pdel_p, sflx_p
       end subroutine aero_model_wetdep_column_flux_codon
    end interface

    if (aero_model_wetdep_use_native_impl) then
       sflx(:) = 0._r8
       do k = 1, pver
          do i = 1, ncol
             sflx(i) = sflx(i) + field(i,k)*pdel(i,k)/gravit
          end do
       end do
       return
    end if

    call aero_model_wetdep_column_flux_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(gravit, c_double), &
         c_loc(field), c_loc(pdel), c_loc(sflx) &
    )

  end subroutine aero_model_wetdep_column_flux

  !-------------------------------------------------------------------------
  ! provides wet tropospheric aerosol surface area info for modal aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  mmr, radmean, relhum, pmid, temp, strato_sad, &
                  sulfate, rho, ltrop, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total )

    ! dummy args
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: rho(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)

    ! local vars
    real(r8), pointer, dimension(:,:,:) :: dgnumwet
    integer :: beglev(ncol)
    integer :: endlev(ncol)
    integer :: i,k

    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet )

    beglev(:ncol)=ltrop(:ncol)
    endlev(:ncol)=pver
    call surf_area_dens( ncol, mmr, pmid, temp, dgnumwet, beglev, endlev, sad_total, sfc=sfc )

    do i = 1,ncol
       do k = ltrop(i),pver
          dm_aer(i,k,:) = dgnumwet(i,k,:) * 1.e2_r8 ! convert m to cm
       enddo
    enddo

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! provides dry stratospheric aerosol surface area info for modal aerosols
  ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad )
    use iso_c_binding, only : c_int64_t

    ! dummy args
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)

    ! local vars
    real(r8), pointer, dimension(:,:,:) :: dgnum
    integer :: beglev(ncol)
    integer :: endlev(ncol)
    integer(c_int64_t) :: active_c

    interface
       function aero_model_strat_surfarea_codon(active_c) result(out_c) bind(c, name="aero_model_strat_surfarea_codon")
         import :: c_int64_t
         integer(c_int64_t), value :: active_c
         integer(c_int64_t) :: out_c
       end function aero_model_strat_surfarea_codon
    end interface

    strato_sad = 0._r8

    active_c = merge(1_c_int64_t, 0_c_int64_t, modal_strat_sulfate)
    if (.not. aero_model_strat_surfarea_codon_logged) then
       aero_model_strat_surfarea_codon_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'aero_model_strat_surfarea direct = codon; inactive/default branch selected'
          call flush(iulog)
       end if
    end if

    if (active_c == 0_c_int64_t) return

    call pbuf_get_field(pbuf, dgnum_idx, dgnum )

    beglev(:ncol)=top_lev
    endlev(:ncol)=ltrop(:ncol)
    call surf_area_dens( ncol, mmr, pmid, temp, dgnum, beglev, endlev, strato_sad )

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use time_manager,          only : get_nstep
    use modal_aero_coag,       only : modal_aero_coag_sub, modal_aero_coag_sub_direct_codon
    use modal_aero_gasaerexch, only : modal_aero_gasaerexch_sub, modal_aero_gasaerexch_sub_direct_codon
    use modal_aero_newnuc,     only : modal_aero_newnuc_sub, modal_aero_newnuc_sub_direct_codon
    use mo_setsox,             only : setsox, has_sox, setsox_shell_codon_wrap
    use modal_aero_data,       only : cnst_name_cw, qqcw_get_field

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(pcols)
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:) 
    real(r8), intent(in) :: zm(:,:) 
    real(r8), intent(in) :: qh2o(:,:) 
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:) 
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    ! local vars 
    
    integer :: n, m
    integer :: i,k
    integer :: nstep
    integer, parameter :: gasaerexch_stage3_save = 0
    integer, parameter :: gasaerexch_stage3_delta = 1
    integer, parameter :: gasaerexch_stage_aq_save = 4

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8) :: wrk(ncol,gas_pcnst)
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: sulfeq(:,:,:)

    call pbuf_get_field(pbuf, dgnum_idx,      dgnum )
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    call pbuf_get_field(pbuf, pblh_idx,       pblh)

    do n=1,ntot_amode
       call outfld(dgnum_name(n), dgnum(1:ncol,1:pver,n), ncol, lchnk )
       call outfld(dgnumwet_name(n), dgnumwet(1:ncol,1:pver,n), ncol, lchnk )
    end do

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()
    call aero_model_gasaerexch_select_impl()

    ! calculate tendency due to gas phase chemistry and processes
    if (aero_model_gasaerexch_use_native_impl) then
       call aero_model_gasaerexch_gas_tend(ncol, delt, vmr0, vmr, dvmrdt)
    else
       call aero_model_gasaerexch_presetsox_codon(vmr0, vmr, dvmrdt, mbar, pdel, ncol, delt, wrk)
    end if
    do m = 1, gas_pcnst
      if (aero_model_gasaerexch_use_native_impl) then
         call aero_model_gasaerexch_column_flux(ncol, dvmrdt(:,:,m), mbar, pdel, adv_mass(m), wrk(:,m))
      end if
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol,m), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    if (aero_model_gasaerexch_use_native_impl) then
       call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )
    else
       call aero_model_gasaerexch_load_vmrcw_codon( lchnk, vmr, vmrcw, mbar, ncol, loffset, pbuf, dvmrdt, dvmrcwdt )
    end if

    if (aero_model_gasaerexch_use_native_impl) then
       dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
       dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)
    end if

  ! aqueous chemistry ...

    if( has_sox ) then
       if (aero_model_gasaerexch_use_native_impl) then
          call setsox(   &
               ncol,     &
               lchnk,    &
               loffset,  &
               delt,     &
               pmid,     &
               pdel,     &
               tfld,     &
               mbar,     &
               cwat,     &
               cldfr,    &
               cldnum,   &
               airdens,  &
               invariants, &
               vmrcw,    &
               vmr       &
               )
       else
          call setsox_shell_codon_wrap( &
               ncol,     &
               lchnk,    &
               loffset,  &
               delt,     &
               pmid,     &
               pdel,     &
               tfld,     &
               mbar,     &
               cwat,     &
               cldfr,    &
               cldnum,   &
               airdens,  &
               invariants, &
               vmrcw,    &
               vmr       &
               )
       end if
    endif

!   Tendency due to aqueous chemistry 
    if (aero_model_gasaerexch_use_native_impl) then
       call aero_model_gasaerexch_aq_tend(ncol, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)
    else
       call aero_model_gasaerexch_codon_wrap( &
            gasaerexch_stage_aq_save, gasaerexch_stage3_save, ncol, delt, ndx_h2so4, vmr0, vmr, vmrcw, dvmrdt, dvmrcwdt, &
            mbar, pdel, adv_mass, wrk, del_h2so4_aeruptk )
    end if
    do m = 1, gas_pcnst
      if (aero_model_gasaerexch_use_native_impl) then
         call aero_model_gasaerexch_column_flux(ncol, dvmrdt(:,:,m), mbar, pdel, adv_mass(m), wrk(:,m))
      end if
      name = 'AQ_'//trim(solsym(m))
      call outfld( name, wrk(:ncol,m), ncol, lchnk )
    enddo

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    if (aero_model_gasaerexch_use_native_impl) then
       call aero_model_gasaerexch_h2so4_save(ncol, ndx_h2so4, vmr, del_h2so4_aeruptk)
    end if

    call t_startf('modal_gas-aer_exchng')
    
    if ( sulfeq_idx>0 ) then
       call pbuf_get_field( pbuf, sulfeq_idx, sulfeq )
    else
       nullify( sulfeq )
    endif

    if (aero_model_gasaerexch_use_native_impl) then
       call modal_aero_gasaerexch_sub(         &
            lchnk,    ncol,     nstep,         &
            loffset,            delt,          &
            tfld,     pmid,     pdel,          &
            qh2o,               troplev,       &
            vmr,                vmrcw,         &
            dvmrdt,             dvmrcwdt,      &
            dgnum,              dgnumwet,      &
            sulfeq     )
    else
       call modal_aero_gasaerexch_sub_direct_codon( &
            lchnk,    ncol,     nstep,              &
            loffset,            delt,               &
            tfld,     pmid,     pdel,               &
            qh2o,               troplev,            &
            vmr,                vmrcw,              &
            dvmrdt,             dvmrcwdt,           &
            dgnum,              dgnumwet,           &
            sulfeq     )
    end if

    if (aero_model_gasaerexch_use_native_impl) then
       call aero_model_gasaerexch_h2so4_delta(ncol, ndx_h2so4, vmr, del_h2so4_aeruptk)
    else
       call aero_model_gasaerexch_codon_wrap( &
            3, gasaerexch_stage3_delta, ncol, delt, ndx_h2so4, vmr0, vmr, vmrcw, dvmrdt, dvmrcwdt, mbar, pdel, adv_mass, wrk, del_h2so4_aeruptk )
    end if

    call t_stopf('modal_gas-aer_exchng')

    call t_startf('modal_nucl')

    ! do aerosol nucleation (new particle formation)
    if (aero_model_gasaerexch_use_native_impl) then
       call modal_aero_newnuc_sub(                          &
            lchnk,    ncol,     nstep,         &
            loffset,            delt,          &
            tfld,     pmid,     pdel,          &
            zm,       pblh,                    &
            qh2o,     cldfr,                   &
            vmr,                               &
            del_h2so4_gasprod,  del_h2so4_aeruptk )
    else
       call modal_aero_newnuc_sub_direct_codon(            &
            lchnk,    ncol,     nstep,         &
            loffset,            delt,          &
            tfld,     pmid,     pdel,          &
            zm,       pblh,                    &
            qh2o,     cldfr,                   &
            vmr,                               &
            del_h2so4_gasprod,  del_h2so4_aeruptk )
    end if

    call t_stopf('modal_nucl')

    call t_startf('modal_coag')

    ! do aerosol coagulation
    if (aero_model_gasaerexch_use_native_impl) then
       call modal_aero_coag_sub(                            &
            lchnk,    ncol,     nstep,         &
            loffset,            delt,          &
            tfld,     pmid,     pdel,          &
            vmr,                               &
            dgnum,              dgnumwet,      &
            wetdens                          )
    else
       call modal_aero_coag_sub_direct_codon(              &
            lchnk,    ncol,     nstep,         &
            loffset,            delt,          &
            tfld,     pmid,     pdel,          &
            vmr,                               &
            dgnum,              dgnumwet,      &
            wetdens                          )
    end if

    call t_stopf('modal_coag')

    if (aero_model_gasaerexch_use_native_impl) then
       call vmr2qqcw( lchnk, vmrcw, mbar, ncol, loffset, pbuf )
    else
       call aero_model_gasaerexch_store_vmrcw_codon( lchnk, vmrcw, mbar, ncol, loffset, pbuf )
    end if

    ! diagnostics for cloud-borne aerosols... 
    do n = 1,pcnst
       fldcw => qqcw_get_field(pbuf,n,lchnk,errorhandle=.true.)
       if(associated(fldcw)) then
          call outfld( cnst_name_cw(n), fldcw(:,:), pcols, lchnk )
       endif
    end do

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_load_vmrcw_codon(lchnk, vmr, vmrcw, mbar, ncol, im, pbuf, dvmrdt, dvmrcwdt)

    use modal_aero_data, only : qqcw_fill_cptrs
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr, c_associated, c_null_ptr

    integer, intent(in) :: lchnk, ncol, im
    real(r8), target, intent(in) :: mbar(ncol,pver)
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: vmrcw(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: dvmrdt(ncol,pver,gas_pcnst), dvmrcwdt(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: m
    type(c_ptr), target :: qqcw_ptrs(pcnst)
    integer(c_int64_t), target :: qqcw_present(pcnst)
    character(len=160) :: proof_line

    interface
       subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, &
            qqcw_offset_c, mbar_ld1_c, delt_c, gravit_c, qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, &
            dvmrdt_p, dvmrcwdt_p, mbar_p, pdel_p, adv_mass_p, wrk_p) bind(c, name="aero_model_gasaerexch_preset_load_stage_dispatch_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, qqcw_offset_c, mbar_ld1_c
         real(c_double), value :: delt_c, gravit_c
         type(c_ptr), value :: qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p
         type(c_ptr), value :: mbar_p, pdel_p, adv_mass_p, wrk_p
       end subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon
    end interface

    call qqcw_fill_cptrs(pbuf, qqcw_ptrs)
    do m = 1, pcnst
       if (c_associated(qqcw_ptrs(m))) then
          qqcw_present(m) = 1_c_int64_t
       else
          qqcw_present(m) = 0_c_int64_t
       end if
    end do

    if (masterproc .and. .not. aero_model_gasaerexch_load_snapshot_proof_written) then
       proof_line = 'aero_model_gasaerexch preset/load stage shell entered (unified load snapshot stage dispatch = codon)'
       write(iulog,'(A)') trim(proof_line)
       call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', trim(proof_line))
       aero_model_gasaerexch_load_snapshot_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_gasaerexch_preset_load_stage_dispatch_codon( &
         2_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(im, c_int64_t), int(ncol, c_int64_t), 0.0_c_double, real(gravit, c_double), &
         c_loc(qqcw_ptrs(1)), c_loc(qqcw_present(1)), c_null_ptr, c_loc(vmr(1,1,1)), c_loc(vmrcw(1,1,1)), &
         c_loc(dvmrdt(1,1,1)), c_loc(dvmrcwdt(1,1,1)), c_loc(mbar(1,1)), c_null_ptr, c_loc(adv_mass(1)), c_null_ptr &
    )

  end subroutine aero_model_gasaerexch_load_vmrcw_codon

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_presetsox_codon(vmr0, vmr, dvmrdt, mbar, pdel, ncol, delt, wrk)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr, c_null_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: delt
    real(r8), target, intent(in) :: vmr0(ncol,pver,gas_pcnst), vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: mbar(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(inout) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: wrk(ncol,gas_pcnst)

    character(len=160) :: proof_line

    interface
       subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, &
            qqcw_offset_c, mbar_ld1_c, delt_c, gravit_c, qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, &
            dvmrdt_p, dvmrcwdt_p, mbar_p, pdel_p, adv_mass_p, wrk_p) bind(c, name="aero_model_gasaerexch_preset_load_stage_dispatch_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, qqcw_offset_c, mbar_ld1_c
         real(c_double), value :: delt_c, gravit_c
         type(c_ptr), value :: qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p
         type(c_ptr), value :: mbar_p, pdel_p, adv_mass_p, wrk_p
       end subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon
    end interface

    if (masterproc .and. .not. aero_model_gasaerexch_presetsox_proof_written) then
       proof_line = 'aero_model_gasaerexch preset/load stage shell entered (unified presetsox stage dispatch = codon)'
       write(iulog,'(A)') trim(proof_line)
       call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', trim(proof_line))
       aero_model_gasaerexch_presetsox_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_gasaerexch_preset_load_stage_dispatch_codon( &
         1_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), 0_c_int64_t, int(pcols, c_int64_t), real(delt, c_double), real(gravit, c_double), &
         c_null_ptr, c_null_ptr, c_loc(vmr0(1,1,1)), c_loc(vmr(1,1,1)), c_null_ptr, c_loc(dvmrdt(1,1,1)), &
         c_null_ptr, c_loc(mbar(1,1)), c_loc(pdel(1,1)), c_loc(adv_mass(1)), c_loc(wrk(1,1)) &
    )

  end subroutine aero_model_gasaerexch_presetsox_codon

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_store_vmrcw_codon(lchnk, vmr, mbar, ncol, im, pbuf)

    use modal_aero_data, only : qqcw_fill_cptrs
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr, c_associated, c_null_ptr

    integer, intent(in) :: lchnk, ncol, im
    real(r8), target, intent(in) :: mbar(ncol,pver)
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: m
    type(c_ptr), target :: qqcw_ptrs(pcnst)
    integer(c_int64_t), target :: qqcw_present(pcnst)
    character(len=160) :: proof_line

    interface
       subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, &
            qqcw_offset_c, mbar_ld1_c, delt_c, gravit_c, qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, &
            dvmrdt_p, dvmrcwdt_p, mbar_p, pdel_p, adv_mass_p, wrk_p) bind(c, name="aero_model_gasaerexch_preset_load_stage_dispatch_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, qqcw_offset_c, mbar_ld1_c
         real(c_double), value :: delt_c, gravit_c
         type(c_ptr), value :: qqcw_ptrs_p, qqcw_present_p, vmr0_p, vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p
         type(c_ptr), value :: mbar_p, pdel_p, adv_mass_p, wrk_p
       end subroutine aero_model_gasaerexch_preset_load_stage_dispatch_codon
    end interface

    call qqcw_fill_cptrs(pbuf, qqcw_ptrs)
    do m = 1, pcnst
       if (c_associated(qqcw_ptrs(m))) then
          qqcw_present(m) = 1_c_int64_t
       else
          qqcw_present(m) = 0_c_int64_t
       end if
    end do

    if (masterproc .and. .not. aero_model_gasaerexch_store_snapshot_proof_written) then
       proof_line = 'aero_model_gasaerexch preset/load/store stage shell entered (unified store snapshot stage dispatch = codon)'
       write(iulog,'(A)') trim(proof_line)
       call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', trim(proof_line))
       aero_model_gasaerexch_store_snapshot_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_gasaerexch_preset_load_stage_dispatch_codon( &
         3_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(im, c_int64_t), int(ncol, c_int64_t), 0.0_c_double, real(gravit, c_double), &
         c_loc(qqcw_ptrs(1)), c_loc(qqcw_present(1)), c_null_ptr, c_loc(vmr(1,1,1)), c_null_ptr, &
         c_null_ptr, c_null_ptr, c_loc(mbar(1,1)), c_null_ptr, c_loc(adv_mass(1)), c_null_ptr &
    )

  end subroutine aero_model_gasaerexch_store_vmrcw_codon

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_append_impl_proof(env_name, proof_line)

    character(len=*), intent(in) :: env_name, proof_line

    character(len=512) :: proof_path
    integer :: status, n, unit_id

    call get_environment_variable(env_name, value=proof_path, length=n, status=status)
    if (status /= 0 .or. n <= 0) return

    open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
         position='append', iostat=status)
    if (status /= 0) return

    write(unit_id,'(A)') trim(proof_line)
    close(unit_id)

  end subroutine aero_model_gasaerexch_append_impl_proof

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_select_impl()

    character(len=48) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_use_native_impl = .false.
    end if

    aero_model_gasaerexch_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch implementation = codon'
          if (.not. aero_model_gasaerexch_proof_written) then
             call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', &
                  'aero_model_gasaerexch selector entered implementation = codon')
             aero_model_gasaerexch_proof_written = .true.
          end if
       end if
       call flush(iulog)
    end if

  end subroutine aero_model_gasaerexch_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_codon_wrap(stage, stage3_mode, ncol, delt, ndx_h2so4_in, &
                                              vmr0, vmr, vmrcw, dvmrdt, dvmrcwdt, mbar, pdel, &
                                              adv_mass_in, wrk, del_h2so4_aeruptk)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr

    integer, intent(in) :: stage, stage3_mode, ncol, ndx_h2so4_in
    real(r8), intent(in) :: delt
    real(r8), target, intent(in) :: vmr0(ncol,pver,gas_pcnst), vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: vmrcw(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: mbar(pcols,pver), pdel(pcols,pver), adv_mass_in(gas_pcnst)
    real(r8), target, intent(inout) :: dvmrdt(ncol,pver,gas_pcnst), dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: wrk(ncol,gas_pcnst), del_h2so4_aeruptk(ncol,pver)
    type(c_ptr) :: vmrcw_p, dvmrcwdt_p, wrk_p, del_h2so4_aeruptk_p
    character(len=160) :: wrap_proof_line

    interface
       subroutine aero_model_gasaerexch_stage_dispatch_codon(stage_c, stage3_mode_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, &
                                              ndx_h2so4_c, delt_c, gravit_c, vmr0_p, vmr_p, vmrcw_p, dvmrdt_p, &
                                              dvmrcwdt_p, mbar_p, pdel_p, adv_mass_p, wrk_p, del_h2so4_aeruptk_p) &
            bind(c, name="aero_model_gasaerexch_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, stage3_mode_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, ndx_h2so4_c
         real(c_double), value :: delt_c, gravit_c
         type(c_ptr), value :: vmr0_p, vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p
         type(c_ptr), value :: mbar_p, pdel_p, adv_mass_p, wrk_p, del_h2so4_aeruptk_p
       end subroutine aero_model_gasaerexch_stage_dispatch_codon
    end interface

    if (stage == 2 .or. stage == 4) then
       vmrcw_p = c_loc(vmrcw(1,1,1))
       dvmrcwdt_p = c_loc(dvmrcwdt(1,1,1))
    else
       vmrcw_p = c_null_ptr
       dvmrcwdt_p = c_null_ptr
    end if

    if (stage <= 2 .or. stage == 4) then
       wrk_p = c_loc(wrk(1,1))
    else
       wrk_p = c_null_ptr
    end if

    if (stage == 3 .or. stage == 4) then
       del_h2so4_aeruptk_p = c_loc(del_h2so4_aeruptk(1,1))
    else
       del_h2so4_aeruptk_p = c_null_ptr
    end if

    if (masterproc .and. .not. aero_model_gasaerexch_wrap_proof_written) then
       wrap_proof_line = 'aero_model_gasaerexch_codon_wrap entered (unified setsox/gaexch/newnuc/coag/qqcw stage dispatch = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', trim(wrap_proof_line))
       aero_model_gasaerexch_wrap_proof_written = .true.
       call flush(iulog)
    end if

    if (masterproc .and. stage == 4 .and. .not. aero_model_gasaerexch_aq_save_proof_written) then
       wrap_proof_line = 'aero_model_gasaerexch aq/save shell entered (unified aq/save stage dispatch = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_gasaerexch_append_impl_proof('AERO_MODEL_GASAEREXCH_PROOF_FILE', trim(wrap_proof_line))
       aero_model_gasaerexch_aq_save_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_gasaerexch_stage_dispatch_codon( &
         int(stage, c_int64_t), int(stage3_mode, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndx_h2so4_in, c_int64_t), real(delt, c_double), &
         real(gravit, c_double), c_loc(vmr0(1,1,1)), c_loc(vmr(1,1,1)), vmrcw_p, c_loc(dvmrdt(1,1,1)), &
         dvmrcwdt_p, c_loc(mbar(1,1)), c_loc(pdel(1,1)), c_loc(adv_mass_in(1)), wrk_p, del_h2so4_aeruptk_p &
    )

  end subroutine aero_model_gasaerexch_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_h2so4_save_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_h2so4_save_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_H2SO4_SAVE_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_h2so4_save_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_h2so4_save_use_native_impl = .false.
    end if

    aero_model_gasaerexch_h2so4_save_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_h2so4_save_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch_h2so4_save implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch_h2so4_save implementation = codon'
       end if
    end if

  end subroutine aero_model_gasaerexch_h2so4_save_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_h2so4_delta_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_h2so4_delta_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_H2SO4_DELTA_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_h2so4_delta_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_h2so4_delta_use_native_impl = .false.
    end if

    aero_model_gasaerexch_h2so4_delta_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_h2so4_delta_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch_h2so4_delta implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch_h2so4_delta implementation = codon'
       end if
    end if

  end subroutine aero_model_gasaerexch_h2so4_delta_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_h2so4_save(ncol, ndx_h2so4_in, vmr, del_h2so4_aeruptk)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol, ndx_h2so4_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: del_h2so4_aeruptk(ncol,pver)

    interface
       subroutine aero_model_gasaerexch_h2so4_save_codon(ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c, vmr_p, del_h2so4_aeruptk_p) &
            bind(c, name="aero_model_gasaerexch_h2so4_save_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c
         type(c_ptr), value :: vmr_p, del_h2so4_aeruptk_p
       end subroutine aero_model_gasaerexch_h2so4_save_codon
    end interface

    call aero_model_gasaerexch_h2so4_save_select_impl()

    if (aero_model_gasaerexch_h2so4_save_use_native_impl) then
       if (ndx_h2so4_in > 0) then
          del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4_in)
       else
          del_h2so4_aeruptk(:,:) = 0.0_r8
       end if
       return
    end if

    call aero_model_gasaerexch_h2so4_save_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndx_h2so4_in, c_int64_t), &
         c_loc(vmr(1,1,1)), c_loc(del_h2so4_aeruptk(1,1)) &
    )

  end subroutine aero_model_gasaerexch_h2so4_save

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_h2so4_delta(ncol, ndx_h2so4_in, vmr, del_h2so4_aeruptk)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol, ndx_h2so4_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: del_h2so4_aeruptk(ncol,pver)

    interface
       subroutine aero_model_gasaerexch_h2so4_delta_codon(ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c, vmr_p, del_h2so4_aeruptk_p) &
            bind(c, name="aero_model_gasaerexch_h2so4_delta_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c
         type(c_ptr), value :: vmr_p, del_h2so4_aeruptk_p
       end subroutine aero_model_gasaerexch_h2so4_delta_codon
    end interface

    call aero_model_gasaerexch_h2so4_delta_select_impl()

    if (aero_model_gasaerexch_h2so4_delta_use_native_impl) then
       if (ndx_h2so4_in > 0) then
          del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4_in) - del_h2so4_aeruptk(1:ncol,:)
       end if
       return
    end if

    call aero_model_gasaerexch_h2so4_delta_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndx_h2so4_in, c_int64_t), &
         c_loc(vmr(1,1,1)), c_loc(del_h2so4_aeruptk(1,1)) &
    )

  end subroutine aero_model_gasaerexch_h2so4_delta

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_gas_tend_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_gas_tend_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_GAS_TEND_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_gas_tend_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_gas_tend_use_native_impl = .false.
    end if

    aero_model_gasaerexch_gas_tend_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_gas_tend_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch_gas_tend implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch_gas_tend implementation = codon'
       end if
    end if

  end subroutine aero_model_gasaerexch_gas_tend_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_gas_tend(ncol, delt, vmr0, vmr, dvmrdt)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: delt
    real(r8), target, intent(in) :: vmr0(ncol,pver,gas_pcnst), vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: dvmrdt(ncol,pver,gas_pcnst)

    interface
       subroutine aero_model_gasaerexch_gas_tend_codon(ncol_c, pver_c, gas_pcnst_c, delt_c, vmr0_p, vmr_p, dvmrdt_p) &
            bind(c, name="aero_model_gasaerexch_gas_tend_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c
         real(c_double), value :: delt_c
         type(c_ptr), value :: vmr0_p, vmr_p, dvmrdt_p
       end subroutine aero_model_gasaerexch_gas_tend_codon
    end interface

    call aero_model_gasaerexch_gas_tend_select_impl()

    if (aero_model_gasaerexch_gas_tend_use_native_impl) then
       dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
       return
    end if

    call aero_model_gasaerexch_gas_tend_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), real(delt, c_double), &
         c_loc(vmr0(1,1,1)), c_loc(vmr(1,1,1)), c_loc(dvmrdt(1,1,1)) &
    )

  end subroutine aero_model_gasaerexch_gas_tend

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_aq_tend_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_aq_tend_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_AQ_TEND_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_aq_tend_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_aq_tend_use_native_impl = .false.
    end if

    aero_model_gasaerexch_aq_tend_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_aq_tend_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch_aq_tend implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch_aq_tend implementation = codon'
       end if
    end if

  end subroutine aero_model_gasaerexch_aq_tend_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_aq_tend(ncol, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: delt
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst), vmrcw(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: dvmrdt(ncol,pver,gas_pcnst), dvmrcwdt(ncol,pver,gas_pcnst)

    interface
       subroutine aero_model_gasaerexch_aq_tend_codon(ncol_c, pver_c, gas_pcnst_c, delt_c, vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p) &
            bind(c, name="aero_model_gasaerexch_aq_tend_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c
         real(c_double), value :: delt_c
         type(c_ptr), value :: vmr_p, vmrcw_p, dvmrdt_p, dvmrcwdt_p
       end subroutine aero_model_gasaerexch_aq_tend_codon
    end interface

    call aero_model_gasaerexch_aq_tend_select_impl()

    if (aero_model_gasaerexch_aq_tend_use_native_impl) then
       dvmrdt = (vmr - dvmrdt) / delt
       dvmrcwdt = (vmrcw - dvmrcwdt) / delt
       return
    end if

    call aero_model_gasaerexch_aq_tend_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), real(delt, c_double), &
         c_loc(vmr(1,1,1)), c_loc(vmrcw(1,1,1)), c_loc(dvmrdt(1,1,1)), c_loc(dvmrcwdt(1,1,1)) &
    )

  end subroutine aero_model_gasaerexch_aq_tend

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_column_flux_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_gasaerexch_column_flux_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_GASAEREXCH_COLUMN_FLUX_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_gasaerexch_column_flux_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_gasaerexch_column_flux_use_native_impl = .false.
    end if

    aero_model_gasaerexch_column_flux_impl_selected = .true.

    if (masterproc) then
       if (aero_model_gasaerexch_column_flux_use_native_impl) then
          write(iulog,*) 'aero_model_gasaerexch_column_flux implementation = native'
       else
          write(iulog,*) 'aero_model_gasaerexch_column_flux implementation = codon'
       end if
    end if

  end subroutine aero_model_gasaerexch_column_flux_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch_column_flux(ncol, field, mbar, pdel, adv_mass_in, wrk)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: field(ncol,pver)
    real(r8), target, intent(in) :: mbar(pcols,pver), pdel(pcols,pver)
    real(r8), intent(in) :: adv_mass_in
    real(r8), target, intent(out) :: wrk(ncol)

    integer :: k

    interface
       subroutine aero_model_gasaerexch_column_flux_codon(ncol_c, pcols_c, pver_c, adv_mass_c, gravit_c, field_p, mbar_p, pdel_p, wrk_p) &
            bind(c, name="aero_model_gasaerexch_column_flux_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: adv_mass_c, gravit_c
         type(c_ptr), value :: field_p, mbar_p, pdel_p, wrk_p
       end subroutine aero_model_gasaerexch_column_flux_codon
    end interface

    call aero_model_gasaerexch_column_flux_select_impl()

    if (aero_model_gasaerexch_column_flux_use_native_impl) then
       wrk(:) = 0._r8
       do k = 1,pver
          wrk(:ncol) = wrk(:ncol) + field(:ncol,k) * adv_mass_in/mbar(:ncol,k)*pdel(:ncol,k)/gravit
       end do
       return
    end if

    call aero_model_gasaerexch_column_flux_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(adv_mass_in, c_double), real(gravit, c_double), &
         c_loc(field), c_loc(mbar), c_loc(pdel), c_loc(wrk) &
    )

  end subroutine aero_model_gasaerexch_column_flux

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_names, seasalt_indices, seasalt_active, seasalt_nbin, seasalt_nnum, &
         seasalt_emis_scale, seasalt_sz_range_lo, seasalt_sz_range_hi
    use dust_model,    only: dust_names, dust_indices, dust_active, dust_nbin, dust_nnum, dust_emis_sclfctr, dust_dmt_vwr
    use physics_types, only: physics_state
    use soil_erod_mod, only: soil_erodibility, soil_erod_fact
    use sslt_sections, only: nsections, Dg, rdry, consta, constb
    use mo_constants,  only: dust_density, dns_aer_sst=>seasalt_density
    use physconst,     only: pi
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr, c_null_ptr

    ! Arguments:

    type(physics_state), target, intent(in)      :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol, ndstflx
    integer :: i, m, mm, isec
    integer(c_int64_t), target :: dust_indices_c(dust_nbin+dust_nnum)
    integer(c_int64_t), target :: seasalt_indices_c(seasalt_nbin+seasalt_nnum)
    real(c_double), target :: dust_emis_sclfctr_c(dust_nbin)
    real(c_double), target :: seasalt_sz_range_lo_c(seasalt_nbin)
    real(c_double), target :: seasalt_sz_range_hi_c(seasalt_nbin)
    real(c_double), target :: dg_c(nsections)
    real(c_double), target :: rdry_c(nsections)
    real(r8), target :: soil_erod_tmp(pcols)
    real(r8), target :: dust_flux_sum(pcols)
    real(r8), target :: dust_sflx(pcols)   ! accumulate over dust bins for output
    real(r8), target :: seasalt_sflx(pcols)   ! accumulate over seasalt bins for output
    real(r8), target :: u10cubed(pcols)
    real(r8), target :: whitecap(pcols)
    real(r8), target :: fi(pcols,nsections)
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model
    real(r8), parameter :: soil_erod_threshold = 0.1_r8

    interface
       function aero_model_emissions_codon(stage_c) result(stage_out) bind(c, name="aero_model_emissions_codon")
         import :: c_int64_t
         integer(c_int64_t), value :: stage_c
         integer(c_int64_t) :: stage_out
       end function aero_model_emissions_codon
    end interface

    call aero_model_emissions_select_impl()

    if (aero_model_emissions_use_native_impl) then
       call aero_model_emissions_native(state, cam_in)
       return
    end if

    if (aero_model_emissions_codon(1_c_int64_t) /= 1_c_int64_t) then
       call endrun('aero_model_emissions :: Codon entry token failed')
    end if

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then
       ndstflx = size(cam_in%dstflx, 2)
       dust_indices_c(:) = int(dust_indices(:), c_int64_t)
       dust_emis_sclfctr_c(:) = real(dust_emis_sclfctr(:), c_double)
    end if

    if (seasalt_active) then
       seasalt_indices_c(:) = int(seasalt_indices(:), c_int64_t)
       seasalt_sz_range_lo_c(:) = real(seasalt_sz_range_lo(:), c_double)
       seasalt_sz_range_hi_c(:) = real(seasalt_sz_range_hi(:), c_double)
       do isec = 1, nsections
          dg_c(isec) = real(Dg(isec), c_double)
          rdry_c(isec) = real(rdry(isec), c_double)
       end do

       call aero_model_emissions_seasalt_wind_select_impl()
       if (aero_model_emissions_seasalt_wind_use_native_impl) then
          call aero_model_emissions_seasalt_wind(ncol, state%u, state%v, state%zm, z0, u10cubed)
       end if
    end if

    if (dust_active .and. seasalt_active) then
       call aero_model_emissions_codon_wrap( &
            3, ncol, ndstflx, real(soil_erod_fact, c_double), real(seasalt_emis_scale, c_double), &
            real(dust_density, c_double), real(dns_aer_sst, c_double), c_loc(cam_in%dstflx(1,1)), c_loc(soil_erod_tmp(1)), &
            c_loc(dust_flux_sum(1)), c_loc(fi(1,1)), c_loc(cam_in%ocnfrac(1)), c_loc(cam_in%cflx(1,1)), &
            c_loc(dust_sflx(1)), c_loc(seasalt_sflx(1)), c_loc(dust_indices_c(1)), c_loc(dust_emis_sclfctr_c(1)), &
            c_loc(dust_dmt_vwr(1)), c_loc(seasalt_indices_c(1)), c_loc(seasalt_sz_range_lo_c(1)), &
            c_loc(seasalt_sz_range_hi_c(1)), c_loc(dg_c(1)), c_loc(rdry_c(1)), c_loc(soil_erodibility(1,lchnk)), &
            real(soil_erod_threshold, c_double), c_loc(cam_in%sst(1)), c_loc(u10cubed(1)), c_loc(whitecap(1)), &
            c_loc(consta(1,1)), c_loc(constb(1,1)), merge(1_c_int64_t, 0_c_int64_t, .not. aero_model_emissions_seasalt_wind_use_native_impl), &
            real(z0, c_double), c_loc(state%u), c_loc(state%v), c_loc(state%zm) &
       )
    else if (dust_active) then
       call aero_model_emissions_codon_wrap( &
            1, ncol, ndstflx, real(soil_erod_fact, c_double), 0._c_double, real(dust_density, c_double), 0._c_double, &
            c_loc(cam_in%dstflx(1,1)), c_loc(soil_erod_tmp(1)), c_loc(dust_flux_sum(1)), c_null_ptr, c_null_ptr, &
            c_loc(cam_in%cflx(1,1)), c_loc(dust_sflx(1)), c_null_ptr, c_loc(dust_indices_c(1)), c_loc(dust_emis_sclfctr_c(1)), &
            c_loc(dust_dmt_vwr(1)), &
            c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_loc(soil_erodibility(1,lchnk)), &
            real(soil_erod_threshold, c_double), c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
            0_c_int64_t, 0._c_double, c_null_ptr, c_null_ptr, c_null_ptr &
       )
    else if (seasalt_active) then
       call aero_model_emissions_codon_wrap( &
            2, ncol, 0, 0._c_double, real(seasalt_emis_scale, c_double), 0._c_double, real(dns_aer_sst, c_double), &
            c_null_ptr, c_null_ptr, c_null_ptr, c_loc(fi(1,1)), c_loc(cam_in%ocnfrac(1)), c_loc(cam_in%cflx(1,1)), &
            c_null_ptr, c_loc(seasalt_sflx(1)), c_null_ptr, c_null_ptr, c_null_ptr, c_loc(seasalt_indices_c(1)), &
            c_loc(seasalt_sz_range_lo_c(1)), c_loc(seasalt_sz_range_hi_c(1)), c_loc(dg_c(1)), c_loc(rdry_c(1)), &
            c_null_ptr, 0._c_double, c_loc(cam_in%sst(1)), c_loc(u10cubed(1)), c_loc(whitecap(1)), c_loc(consta(1,1)), &
            c_loc(constb(1,1)), merge(1_c_int64_t, 0_c_int64_t, .not. aero_model_emissions_seasalt_wind_use_native_impl), &
            real(z0, c_double), c_loc(state%u), c_loc(state%v), c_loc(state%zm) &
       )
    end if

    if (dust_active) then
       ! some dust emis diagnostics ...
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',dust_sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       do m=1,seasalt_nbin
          mm = seasalt_indices(m)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       call outfld('SSTSFMBL',seasalt_sflx(:),pcols,lchnk)
    endif

  end subroutine aero_model_emissions

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_native( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_names, seasalt_indices, seasalt_active,seasalt_nbin
    use dust_model,    only: dust_emis, dust_names, dust_indices, dust_active,dust_nbin, dust_nnum
    use physics_types, only: physics_state

    type(physics_state),    intent(in)    :: state
    type(cam_in_t),         intent(inout) :: cam_in

    integer :: lchnk, ncol
    integer :: m, mm
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)
    real(r8) :: u10cubed(pcols)
    real (r8), parameter :: z0=0.0001_r8

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then
       call dust_emis( ncol, lchnk, cam_in%dstflx, cam_in%cflx, soil_erod_tmp )
       call aero_model_emissions_accumulate_sflx(ncol, dust_nbin, dust_indices(1:dust_nbin), cam_in%cflx, sflx)
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       call aero_model_emissions_seasalt_wind(ncol, state%u, state%v, state%zm, z0, u10cubed)
       call seasalt_emis( u10cubed, cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx )
       call aero_model_emissions_accumulate_sflx(ncol, seasalt_nbin, seasalt_indices(1:seasalt_nbin), cam_in%cflx, sflx)
       do m=1,seasalt_nbin
          mm = seasalt_indices(m)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       call outfld('SSTSFMBL',sflx(:),pcols,lchnk)
    endif

  end subroutine aero_model_emissions_native

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_append_impl_proof(env_name, proof_line)

    character(len=*), intent(in) :: env_name, proof_line

    character(len=512) :: proof_path
    integer :: status, n, unit_id

    call get_environment_variable(env_name, value=proof_path, length=n, status=status)
    if (status /= 0 .or. n <= 0) return

    open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
         position='append', iostat=status)
    if (status /= 0) return

    write(unit_id,'(A)') trim(proof_line)
    close(unit_id)

  end subroutine aero_model_emissions_append_impl_proof

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_emissions_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_EMISSIONS_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_emissions_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_emissions_use_native_impl = .false.
    end if

    aero_model_emissions_impl_selected = .true.

    if (masterproc) then
       if (aero_model_emissions_use_native_impl) then
          write(iulog,*) 'aero_model_emissions implementation = native'
       else
          write(iulog,*) 'aero_model_emissions implementation = codon'
          if (.not. aero_model_emissions_proof_written) then
             call aero_model_emissions_append_impl_proof('AERO_MODEL_EMISSIONS_PROOF_FILE', &
                  'aero_model_emissions selector entered implementation = codon')
             aero_model_emissions_proof_written = .true.
          end if
       end if
       call flush(iulog)
    end if

  end subroutine aero_model_emissions_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_codon_wrap(stage, ncol, ndstflx, soil_erod_fact_in, seasalt_emis_scale_in, dust_density_in, &
                                             seasalt_density_in, dstflx_p, soil_erod_p, dust_flux_sum_p, fi_p, ocnfrac_p, cflx_p, &
                                             dust_sflx_p, seasalt_sflx_p, dust_indices_p, &
                                             dust_emis_sclfctr_p, dust_x_mton_p, seasalt_indices_p, seasalt_sz_range_lo_p, &
                                             seasalt_sz_range_hi_p, dg_p, rdry_p, soil_erodibility_p, soil_erod_threshold_in, &
                                             sst_p, u10cubed_p, whitecap_p, consta_p, constb_p, compute_wind_c, z0_in, &
                                             state_u_p, state_v_p, state_zm_p)

    use dust_model,       only: dust_nbin
    use seasalt_model,    only: seasalt_nbin
    use sslt_sections,    only: nsections
    use physconst,        only: pi
    use iso_c_binding,    only: c_double, c_int64_t, c_ptr

    integer, intent(in) :: stage, ncol, ndstflx
    real(c_double), intent(in) :: soil_erod_fact_in, seasalt_emis_scale_in, dust_density_in, seasalt_density_in
    real(c_double), intent(in) :: soil_erod_threshold_in, z0_in
    integer(c_int64_t), intent(in) :: compute_wind_c
    type(c_ptr), value :: dstflx_p, soil_erod_p, dust_flux_sum_p, fi_p, ocnfrac_p, cflx_p
    type(c_ptr), value :: dust_sflx_p, seasalt_sflx_p, dust_indices_p, dust_emis_sclfctr_p
    type(c_ptr), value :: dust_x_mton_p, seasalt_indices_p, seasalt_sz_range_lo_p, seasalt_sz_range_hi_p
    type(c_ptr), value :: dg_p, rdry_p, soil_erodibility_p
    type(c_ptr), value :: sst_p, u10cubed_p, whitecap_p, consta_p, constb_p
    type(c_ptr), value :: state_u_p, state_v_p, state_zm_p

    character(len=256) :: wrap_proof_line

    interface
       subroutine aero_model_emissions_shell_wind_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, ndstflx_c, dust_nbin_c, seasalt_nbin_c, &
            nsections_c, soil_erod_fact_c, seasalt_emis_scale_c, pi_c, dust_density_c, seasalt_density_c, dstflx_p, soil_erod_p, &
            dust_flux_sum_p, fi_p, ocnfrac_p, cflx_p, dust_sflx_p, seasalt_sflx_p, dust_indices_p, dust_emis_sclfctr_p, &
            dust_x_mton_p, seasalt_indices_p, &
            seasalt_sz_range_lo_p, seasalt_sz_range_hi_p, dg_p, rdry_p, soil_erodibility_p, soil_erod_threshold_c, &
            sst_p, u10cubed_p, whitecap_p, consta_p, constb_p, compute_wind_c, z0_c, state_u_p, state_v_p, state_zm_p) &
            bind(c, name="aero_model_emissions_shell_wind_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, ndstflx_c, dust_nbin_c, seasalt_nbin_c, nsections_c
         real(c_double), value :: soil_erod_fact_c, seasalt_emis_scale_c, pi_c, dust_density_c, seasalt_density_c
         real(c_double), value :: soil_erod_threshold_c, z0_c
         integer(c_int64_t), value :: compute_wind_c
         type(c_ptr), value :: dstflx_p, soil_erod_p, dust_flux_sum_p, fi_p, ocnfrac_p, cflx_p
         type(c_ptr), value :: dust_sflx_p, seasalt_sflx_p, dust_indices_p, dust_emis_sclfctr_p
         type(c_ptr), value :: dust_x_mton_p, seasalt_indices_p, seasalt_sz_range_lo_p, seasalt_sz_range_hi_p
         type(c_ptr), value :: dg_p, rdry_p, soil_erodibility_p
         type(c_ptr), value :: sst_p, u10cubed_p, whitecap_p, consta_p, constb_p
         type(c_ptr), value :: state_u_p, state_v_p, state_zm_p
       end subroutine aero_model_emissions_shell_wind_stage_dispatch_codon
    end interface

    if (masterproc .and. .not. aero_model_emissions_wrap_proof_written) then
       wrap_proof_line = 'aero_model_emissions_codon_wrap entered (parent shell = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_emissions_append_impl_proof('AERO_MODEL_EMISSIONS_PROOF_FILE', trim(wrap_proof_line))
       aero_model_emissions_wrap_proof_written = .true.
       call flush(iulog)
    end if

    if (masterproc .and. stage == 1 .and. .not. aero_model_emissions_dust_stage_proof_written) then
       wrap_proof_line = 'aero_model_emissions dust shell entered (unified dust stage dispatch = codon)'
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_emissions_append_impl_proof('AERO_MODEL_EMISSIONS_PROOF_FILE', trim(wrap_proof_line))
       aero_model_emissions_dust_stage_proof_written = .true.
       call flush(iulog)
    else if (masterproc .and. stage == 2 .and. .not. aero_model_emissions_seasalt_stage_proof_written) then
       if (aero_model_emissions_seasalt_wind_use_native_impl) then
          wrap_proof_line = 'aero_model_emissions seasalt shell entered (unified seasalt stage dispatch = codon; wind = native)'
       else
          wrap_proof_line = 'aero_model_emissions seasalt shell entered (unified wind/seasalt stage dispatch = codon)'
       end if
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_emissions_append_impl_proof('AERO_MODEL_EMISSIONS_PROOF_FILE', trim(wrap_proof_line))
       aero_model_emissions_seasalt_stage_proof_written = .true.
       call flush(iulog)
    else if (masterproc .and. stage == 3 .and. .not. aero_model_emissions_all_stage_proof_written) then
       if (compute_wind_c /= 0_c_int64_t) then
          wrap_proof_line = 'aero_model_emissions all shell entered (unified wind/dust/seasalt stage dispatch = codon; outfld = native)'
       else
          wrap_proof_line = 'aero_model_emissions all shell entered (unified dust/seasalt stage dispatch = codon; wind = native; outfld = native)'
       end if
       write(iulog,'(A)') trim(wrap_proof_line)
       call aero_model_emissions_append_impl_proof('AERO_MODEL_EMISSIONS_PROOF_FILE', trim(wrap_proof_line))
       aero_model_emissions_all_stage_proof_written = .true.
       call flush(iulog)
    end if

    call aero_model_emissions_shell_wind_stage_dispatch_codon( &
         int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(ndstflx, c_int64_t), int(dust_nbin, c_int64_t), &
         int(seasalt_nbin, c_int64_t), int(nsections, c_int64_t), soil_erod_fact_in, seasalt_emis_scale_in, real(pi, c_double), &
         dust_density_in, seasalt_density_in, dstflx_p, soil_erod_p, dust_flux_sum_p, fi_p, ocnfrac_p, cflx_p, dust_sflx_p, &
         seasalt_sflx_p, dust_indices_p, &
         dust_emis_sclfctr_p, dust_x_mton_p, seasalt_indices_p, seasalt_sz_range_lo_p, seasalt_sz_range_hi_p, dg_p, rdry_p, &
         soil_erodibility_p, soil_erod_threshold_in, sst_p, u10cubed_p, whitecap_p, consta_p, constb_p, &
         compute_wind_c, z0_in, state_u_p, state_v_p, state_zm_p &
    )

  end subroutine aero_model_emissions_codon_wrap

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_accumulate_sflx_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_emissions_accumulate_sflx_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_EMISSIONS_ACCUMULATE_SFLX_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_emissions_accumulate_sflx_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_emissions_accumulate_sflx_use_native_impl = .false.
    end if

    aero_model_emissions_accumulate_sflx_impl_selected = .true.

    if (masterproc) then
       if (aero_model_emissions_accumulate_sflx_use_native_impl) then
          write(iulog,*) 'aero_model_emissions_accumulate_sflx implementation = native'
       else
          write(iulog,*) 'aero_model_emissions_accumulate_sflx implementation = codon'
       end if
    end if

  end subroutine aero_model_emissions_accumulate_sflx_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_accumulate_sflx(ncol, nindices, indices, cflx, sflx)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol, nindices
    integer, intent(in) :: indices(nindices)
    real(r8), target, intent(in) :: cflx(pcols,pcnst)
    real(r8), target, intent(out) :: sflx(pcols)

    integer :: m
    integer(c_int64_t), target :: indices_c(nindices)

    interface
       subroutine aero_model_emissions_accumulate_sflx_codon(ncol_c, pcols_c, nindices_c, indices_p, cflx_p, sflx_p) &
            bind(c, name="aero_model_emissions_accumulate_sflx_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, nindices_c
         type(c_ptr), value :: indices_p, cflx_p, sflx_p
       end subroutine aero_model_emissions_accumulate_sflx_codon
    end interface

    call aero_model_emissions_accumulate_sflx_select_impl()

    if (aero_model_emissions_accumulate_sflx_use_native_impl) then
       sflx(:) = 0._r8
       do m = 1, nindices
          sflx(:ncol) = sflx(:ncol) + cflx(:ncol,indices(m))
       end do
       return
    end if

    indices_c(:) = int(indices(:), c_int64_t)
    call aero_model_emissions_accumulate_sflx_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(nindices, c_int64_t), &
         c_loc(indices_c), c_loc(cflx), c_loc(sflx) &
    )

  end subroutine aero_model_emissions_accumulate_sflx

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_seasalt_wind_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (aero_model_emissions_seasalt_wind_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('AERO_MODEL_EMISSIONS_SEASALT_WIND_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       aero_model_emissions_seasalt_wind_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       aero_model_emissions_seasalt_wind_use_native_impl = .false.
    end if

    aero_model_emissions_seasalt_wind_impl_selected = .true.

    if (masterproc) then
       if (aero_model_emissions_seasalt_wind_use_native_impl) then
          write(iulog,*) 'aero_model_emissions_seasalt_wind implementation = native'
       else
          write(iulog,*) 'aero_model_emissions_seasalt_wind implementation = codon'
       end if
    end if

  end subroutine aero_model_emissions_seasalt_wind_select_impl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions_seasalt_wind(ncol, state_u, state_v, state_zm, z0, u10cubed)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_u(pcols,pver), state_v(pcols,pver), state_zm(pcols,pver)
    real(r8), intent(in) :: z0
    real(r8), target, intent(out) :: u10cubed(pcols)

    interface
       subroutine aero_model_emissions_seasalt_wind_codon(ncol_c, pcols_c, pver_c, z0_c, state_u_p, state_v_p, state_zm_p, u10cubed_p) &
            bind(c, name="aero_model_emissions_seasalt_wind_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: z0_c
         type(c_ptr), value :: state_u_p, state_v_p, state_zm_p, u10cubed_p
       end subroutine aero_model_emissions_seasalt_wind_codon
    end interface

    call aero_model_emissions_seasalt_wind_select_impl()

    if (aero_model_emissions_seasalt_wind_use_native_impl) then
       u10cubed(:ncol)=sqrt(state_u(:ncol,pver)**2+state_v(:ncol,pver)**2)
       ! move the winds to 10m high from the midpoint of the gridbox:
       ! follows Tie and Seinfeld and Pandis, p.859 with math.
       u10cubed(:ncol)=u10cubed(:ncol)*log(10._r8/z0)/log(state_zm(:ncol,pver)/z0)
       ! we need them to the 3.41 power, according to Gong et al., 1997:
       u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8
       return
    end if

    call aero_model_emissions_seasalt_wind_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(z0, c_double), &
         c_loc(state_u), c_loc(state_v), c_loc(state_zm), c_loc(u10cubed) &
    )

  end subroutine aero_model_emissions_seasalt_wind

  !===============================================================================
  ! private methods


  !=============================================================================
  !=============================================================================
  subroutine surf_area_dens( ncol, mmr, pmid, temp, diam, beglev, endlev, sad, sfc )
    use mo_constants,    only : pi
    use modal_aero_data, only : nspec_amode, alnsg_amode

    ! dummy args
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: mmr(:,:,:)
    real(r8), intent(in)  :: pmid(:,:)
    real(r8), intent(in)  :: temp(:,:)
    real(r8), intent(in)  :: diam(:,:,:)
    integer,  intent(in)  :: beglev(:)
    integer,  intent(in)  :: endlev(:)
    real(r8), intent(out) :: sad(:,:)
    real(r8),optional, intent(out) :: sfc(:,:,:)

    ! local vars
    real(r8) :: sad_mode(pcols,pver,ntot_amode)
    real(r8) :: rho_air
    integer  :: i,k,l,m 
    real(r8) :: chm_mass, tot_mass

    !
    ! Compute surface aero for each mode.
    ! Total over all modes as the surface area for chemical reactions.
    !

    sad = 0._r8
    sad_mode = 0._r8

    do i = 1,ncol
       do k = beglev(i),endlev(i)
          rho_air = pmid(i,k)/(temp(i,k)*287.04_r8)
          do l=1,ntot_amode
             !
             ! compute a mass weighting of the number
             !
             tot_mass = 0._r8
             chm_mass = 0._r8
             do m=1,nspec_amode(l)
               if ( index_tot_mass(l,m) > 0 ) &
                    tot_mass = tot_mass + mmr(i,k,index_tot_mass(l,m))
               if ( index_chm_mass(l,m) > 0 ) &
                    chm_mass = chm_mass + mmr(i,k,index_chm_mass(l,m))
             end do
             if ( tot_mass > 0._r8 ) then
               sad_mode(i,k,l) = chm_mass/tot_mass * &
                    mmr(i,k,num_idx(l))*rho_air*pi*diam(i,k,l)**2*&
                    exp(2*alnsg_amode(l)**2)  ! m^2/m^3
               sad_mode(i,k,l) = 1.e-2_r8 * sad_mode(i,k,l) ! cm^2/cm^3
             else
               sad_mode(i,k,l) = 0._r8
             end if
          end do
          sad(i,k) = sum(sad_mode(i,k,:))

       enddo
    enddo

    if (present(sfc)) then
       sfc(:,:,:) = sad_mode(:,:,:) 
    endif

  end subroutine surf_area_dens

  !===============================================================================
  !===============================================================================
  subroutine modal_aero_bcscavcoef_init
    use iso_c_binding, only : c_int64_t
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    !
    !-----------------------------------------------------------------------
    
    use shr_kind_mod,    only: r8 => shr_kind_r8
    use modal_aero_data
    use cam_abortutils,  only: endrun

    implicit none

    interface
       function modal_aero_bcscavcoef_init_codon(tag) result(tag_out) bind(c, name='modal_aero_bcscavcoef_init_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function modal_aero_bcscavcoef_init_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.



    !   local variables
    integer nnfit_maxd
    parameter (nnfit_maxd=27)

    integer i, jgrow, jdens, jpress, jtemp, ll, mode, nnfit
    integer lunerr

    real(r8) dg0, dg0_cgs, press, &
         rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
         scavratenum, scavratevol, sigmag,                &
         temp, wetdiaratio, wetvolratio
    real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
    real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('MODAL_AERO_BCSCAVCOEF_INIT_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = modal_aero_bcscavcoef_init_codon(int(168, c_int64_t))
       if (rt_codon_tag_out /= int(168, c_int64_t)) then
          write(iulog,*) 'modal_aero_bcscavcoef_init_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'modal_aero_bcscavcoef_init implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

    lunerr = 6
    dlndg_nimptblgrow = log( 1.25_r8 )

    modeloop: do mode = 1, ntot_amode

       sigmag = sigmag_amode(mode)

       ll = lspectype_amode(1,mode)
       rhodryaero = specdens_amode(ll)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          dg0 = dgnum_amode(mode)*wetdiaratio

          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

          !
          !   compute impaction scavenging rates at 1 temp-press pair and save
          !
          nnfit = 0

          temp = 273.16_r8
          press = 0.75e6_r8   ! dynes/cm2
          rhowetaero = rhodryaero

          dg0_cgs = dg0*1.0e2_r8   ! m to cm
          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3
          call calc_1_impact_rate( &
               dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
               scavratenum, scavratevol, lunerr )

          nnfit = nnfit + 1
          if (nnfit .gt. nnfit_maxd) then
             write(lunerr,9110)
             call endrun()
          end if
9110      format( '*** subr. modal_aero_bcscavcoef_init -- nnfit too big' )

          xxfitnum(1,nnfit) = 1._r8
          yyfitnum(nnfit) = log( scavratenum )

          xxfitvol(1,nnfit) = 1._r8
          yyfitvol(nnfit) = log( scavratevol )

5900      continue

          !
          ! skip mlinfit stuff because scav table no longer has dependencies on
          !    air temp, air press, and particle wet density
          ! just load the log( scavrate--- ) values
          !
          !!
          !!   do linear regression
          !!	log(scavrate) = a1 + a2*log(wetdens)
          !!
          !	call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 1, 1, rmserr )
          !	call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 1, 1, rmserr )
          !
          !	scavimptblnum(jgrow,mode) = aafitnum(1)
          !	scavimptblvol(jgrow,mode) = aafitvol(1)

          scavimptblnum(jgrow,mode) = yyfitnum(1)
          scavimptblvol(jgrow,mode) = yyfitvol(1)

       enddo growloop
    enddo modeloop
    return
  end subroutine modal_aero_bcscavcoef_init

  !===============================================================================
  !===============================================================================
  subroutine modal_aero_depvel_part( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
                                     radius_part, density_part, sig_part, moment, fraction_landuse_lcl )

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    !
    ! !USES
    !
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst,     only: pi,boltz, gravit, rair
    use mo_drydep,     only: n_land_type

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8), target, intent(in) :: t(pcols,pver)       !atm temperature (K)
    real(r8), target, intent(in) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8), target, intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8), target, intent(in) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8), target, intent(in) :: radius_part(pcols,pver)    ! mean (volume/number) particle radius (m)
    real(r8), target, intent(in) :: density_part(pcols,pver)   ! density of particle material (kg/m3)
    real(r8), target, intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particles
    integer,  intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol
    real(r8), target, intent(in) :: fraction_landuse_lcl(pcols,n_land_type)

    real(r8), target, intent(out) :: vlc_trb(pcols)       !Turbulent deposn velocity (m/s)
    real(r8), target, intent(out) :: vlc_grv(pcols,pver)  !grav deposn velocity (m/s)
    real(r8), target, intent(out) :: vlc_dry(pcols,pver)  !dry deposn velocity (m/s)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: i,k,ix                !indices
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8), target :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8), target :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8), target :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8), target :: slp_crc(pcols,pver) ![frc] Slip correction factor
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface
    real(r8), target :: radius_moment(pcols,pver) ! median radius (m) for moment
    real(r8) :: lnsig         ! ln(sig_part)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt
    real(r8) :: lnd_frc
    real(r8) :: wrk1

    interface
       subroutine modal_aero_depvel_part_codon(ncol_c, pcols_c, pver_c, n_land_type_c, moment_c, pi_c, boltz_c, &
            gravit_c, rair_c, t_p, pmid_p, ram1_p, fv_p, vlc_dry_p, vlc_trb_p, vlc_grv_p, radius_part_p, &
            density_part_p, sig_part_p, fraction_landuse_p, vsc_dyn_atm_p, vsc_knm_atm_p, mfp_atm_p, slp_crc_p, &
            radius_moment_p) bind(c, name="modal_aero_depvel_part_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, n_land_type_c, moment_c
         real(c_double), value :: pi_c, boltz_c, gravit_c, rair_c
         type(c_ptr), value :: t_p, pmid_p, ram1_p, fv_p, vlc_dry_p, vlc_trb_p, vlc_grv_p, radius_part_p
         type(c_ptr), value :: density_part_p, sig_part_p, fraction_landuse_p, vsc_dyn_atm_p, vsc_knm_atm_p
         type(c_ptr), value :: mfp_atm_p, slp_crc_p, radius_moment_p
       end subroutine modal_aero_depvel_part_codon
    end interface

    if (.not. aero_model_drydep_use_native_impl) then
       call modal_aero_depvel_part_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(n_land_type, c_int64_t), &
            int(moment, c_int64_t), real(pi, c_double), real(boltz, c_double), real(gravit, c_double), &
            real(rair, c_double), c_loc(t), c_loc(pmid), c_loc(ram1), c_loc(fv), c_loc(vlc_dry), c_loc(vlc_trb), &
            c_loc(vlc_grv), c_loc(radius_part), c_loc(density_part), c_loc(sig_part), c_loc(fraction_landuse_lcl), &
            c_loc(vsc_dyn_atm), c_loc(vsc_knm_atm), c_loc(mfp_atm), c_loc(slp_crc), c_loc(radius_moment) &
       )
       return
    end if

    call modal_aero_depvel_part_native( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
         radius_part, density_part, sig_part, moment, fraction_landuse_lcl )

    return
  end subroutine modal_aero_depvel_part

  !===============================================================================
  !===============================================================================
  subroutine modal_aero_depvel_part_native( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
                                            radius_part, density_part, sig_part, moment, fraction_landuse_lcl )

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    !
    ! !USES
    !
    use physconst,     only: pi,boltz, gravit, rair
    use mo_drydep,     only: n_land_type

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8), intent(in) :: t(pcols,pver)       !atm temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8), intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume/number) particle radius (m)
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material (kg/m3)
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particles
    integer,  intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: fraction_landuse_lcl(pcols,n_land_type)

    real(r8), intent(out) :: vlc_trb(pcols)       !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(pcols,pver)       !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(pcols,pver)       !dry deposn velocity (m/s)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: i,k,ix                !indices
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: slp_crc(pcols,pver) ![frc] Slip correction factor
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface
    real(r8) :: radius_moment(pcols,pver) ! median radius (m) for moment
    real(r8) :: lnsig         ! ln(sig_part)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt
    real(r8) :: lnd_frc
    real(r8) :: wrk1, wrk2, wrk3

    ! constants
    real(r8) gamma(11)      ! exponent of schmidt number
!   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
!              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
!              0.50d+00/
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
!   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
!               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
!             100.00d+00/
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
!   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
!                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
!                         -1.00d+00/
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
!   data iwet/1,   -1,   -1,   -1,   -1,  &
!            -1,   -1,   -1,    1,   -1,  &
!             1/
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet


    !------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol

          lnsig = log(sig_part(i,k))
! use a maximum radius of 50 microns when calculating deposition velocity
          radius_moment(i,k) = min(50.0e-6_r8,radius_part(i,k))*   &
                          exp((float(moment)-1.5_r8)*lnsig*lnsig)
          dispersion = exp(2._r8*lnsig*lnsig)

          rho=pmid(i,k)/rair/t(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(i,k) = 1.0_r8 + mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment(i,k)/(mfp_atm(i,k)))) / &
                  radius_moment(i,k)   ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(i,k) = (4.0_r8/18.0_r8) * radius_moment(i,k)*radius_moment(i,k)*density_part(i,k)* &
                  gravit*slp_crc(i,k) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(i,k) = vlc_grv(i,k) * dispersion

          vlc_dry(i,k)=vlc_grv(i,k)
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       dff_aer = boltz * t(i,k) * slp_crc(i,k) / &    ![m2 s-1]
                 (6.0_r8*pi*vsc_dyn_atm(i,k)*radius_moment(i,k)) !SeP97 p.474
       shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972

       wrk2 = 0._r8
       wrk3 = 0._r8
       do lt = 1,n_land_type
          lnd_frc = fraction_landuse_lcl(i,lt)
          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(lt))
             if (radius_collector(lt) > 0.0_r8) then
!       vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) / (gravit*radius_collector(lt))
                interception = 2.0_r8*(radius_moment(i,k)/radius_collector(lt))**2.0_r8
             else
!       non-vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**2.0_r8   

             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(i) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i,k)

             wrk1 = 1.0_r8 / rss_trb
             wrk2 = wrk2 + lnd_frc*( wrk1 )
             wrk3 = wrk3 + lnd_frc*( wrk1 + vlc_grv(i,k) )
          endif
       enddo  ! n_land_type
       vlc_trb(i) = wrk2
       vlc_dry(i,k) = wrk3
    enddo !ncol

    return
  end subroutine modal_aero_depvel_part_native

  !===============================================================================
  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol )

    use modal_aero_data
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    !-----------------------------------------------------------------------
    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    real(r8), target, intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), target, intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

    integer i, k, jgrow
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum
    integer(c_int64_t), target :: isprx_mask(pcols,pver)
    real(r8), target :: scavimptblnum_mode(nimptblgrow_mind:nimptblgrow_maxd)
    real(r8), target :: scavimptblvol_mode(nimptblgrow_mind:nimptblgrow_maxd)
    real(r8) :: dgnum_mode

    interface
       subroutine modal_aero_bcscavcoef_get_codon(m_c, ncol_c, pcols_c, pver_c, ntot_amode_c, nimptblgrow_mind_c, &
            nimptblgrow_maxd_c, dlndg_nimptblgrow_c, dgnum_mode_c, isprx_mask_p, dgn_awet_p, scavimptblnum_mode_p, &
            scavimptblvol_mode_p, scavcoefnum_p, scavcoefvol_p) bind(c, name="modal_aero_bcscavcoef_get_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: m_c, ncol_c, pcols_c, pver_c, ntot_amode_c, nimptblgrow_mind_c
         integer(c_int64_t), value :: nimptblgrow_maxd_c
         real(c_double), value :: dlndg_nimptblgrow_c, dgnum_mode_c
         type(c_ptr), value :: isprx_mask_p, dgn_awet_p, scavimptblnum_mode_p, scavimptblvol_mode_p
         type(c_ptr), value :: scavcoefnum_p, scavcoefvol_p
       end subroutine modal_aero_bcscavcoef_get_codon
    end interface

    if (.not. aero_model_wetdep_use_native_impl) then
       do k = 1, pver
          do i = 1, ncol
             isprx_mask(i,k) = merge(1_c_int64_t, 0_c_int64_t, isprx(i,k))
          end do
       end do
       dgnum_mode = dgnum_amode(m)
       scavimptblnum_mode(:) = scavimptblnum(:,m)
       scavimptblvol_mode(:) = scavimptblvol(:,m)

       call modal_aero_bcscavcoef_get_codon( &
            int(m, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
            int(ntot_amode, c_int64_t), int(nimptblgrow_mind, c_int64_t), int(nimptblgrow_maxd, c_int64_t), &
            real(dlndg_nimptblgrow, c_double), real(dgnum_mode, c_double), c_loc(isprx_mask), c_loc(dgn_awet), &
            c_loc(scavimptblnum_mode(nimptblgrow_mind)), c_loc(scavimptblvol_mode(nimptblgrow_mind)), &
            c_loc(scavcoefnum), c_loc(scavcoefvol) &
       )
       return
    end if


    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) ) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

             if ((dumdgratio .ge. 0.99_r8) .and. (dumdgratio .le. 1.01_r8)) then
                scavimpvol = scavimptblvol(0,m)
                scavimpnum = scavimptblnum(0,m)
             else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow .lt. 0._r8) jgrow = jgrow - 1
                if (jgrow .lt. nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                end if

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                     dumfhi*scavimptblvol(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                     dumfhi*scavimptblnum(jgrow+1,m)

             end if

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          end if

       end do
    end do

    return
  end subroutine modal_aero_bcscavcoef_get

  !===============================================================================
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )
   !
   !   routine computes a single impaction scavenging rate
   !	for precipitation rate of 1 mm/h
   !
   !   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
   !   sigmag = geometric standard deviation of size distrib.
   !   rhoaero = density of aerosol particles (g/cm^3)
   !   temp = temperature (K)
   !   press = pressure (dyne/cm^2)
   !   scavratenum = number scavenging rate (1/h)
   !   scavratevol = volume or mass scavenging rate (1/h)
   !   lunerr = logical unit for error message
   !
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, &
                           gravity => gravity_cgs, rgas => rgas_cgs

   implicit none

   !   subr. parameters
   integer lunerr
   real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

   !   local variables
   integer nrainsvmax
   parameter (nrainsvmax=50)
   real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
        vfallrainsv(nrainsvmax)

   integer naerosvmax
   parameter (naerosvmax=51)
   real(r8) aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

   integer i, ja, jr, na, nr
   real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
   real(r8) anumsum, avolsum, cair, chi
   real(r8) d, dr, dum, dumfuchs, dx
   real(r8) ebrown, eimpact, eintercept, etotal, freepath
   real(r8) precip, precipmmhr, precipsum
   real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
   real(r8) scavsumnum, scavsumnumbb
   real(r8) scavsumvol, scavsumvolbb
   real(r8) schmidt, sqrtreynolds, sstar, stokes, sx              
   real(r8) taurelax, vfall, vfallstp
   real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair                     

   
   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr .gt. nrainsvmax) then
      write(lunerr,9110)
      call endrun()
   end if

9110 format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

   precipmmhr = 1.0_r8
   precip = precipmmhr/36000._r8

   ag0 = dg0/2._r8
   sx = log( sigmag )
   xg0 = log( ag0 )
   xg3 = xg0 + 3._r8*sx*sx

   xlo = xg3 - 4._r8*sx
   xhi = xg3 + 4._r8*sx
   dx = 0.2_r8*sx

   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

   na = 1 + nint( (xhi-xlo)/dx )
   if (na .gt. naerosvmax) then
      write(lunerr,9120)
      call endrun()
   end if

9120 format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

   !   air molar density
   cair = press/(rgas*temp)
   !   air mass density
   rhoair = 28.966_r8*cair
   !   molecular freepath
   freepath = 2.8052e-10_r8/cair
   !   air dynamic viscosity
   airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
        ((temp/296.16_r8)**1.5_r8)
   !   air kinemaic viscosity
   airkinvisc = airdynvisc/rhoair
   !   ratio of water viscosity to air viscosity (from Slinn)
   xmuwaterair = 60.0_r8

   !
   !   compute rain drop number concentrations
   !	rrainsv = raindrop radius (cm)
   !	xnumrainsv = raindrop number concentration (#/cm^3)
   !		(number in the bin, not number density)
   !	vfallrainsv = fall velocity (cm/s)
   !
   precipsum = 0._r8
   do i = 1, nr
      r = rlo + (i-1)*dr
      rrainsv(i) = r
      xnumrainsv(i) = exp( -r/2.7e-2_r8 )

      d = 2._r8*r
      if (d .le. 0.007_r8) then
         vfallstp = 2.88e5_r8 * d**2._r8
      else if (d .le. 0.025_r8) then
         vfallstp = 2.8008e4_r8 * d**1.528_r8
      else if (d .le. 0.1_r8) then
         vfallstp = 4104.9_r8 * d**1.008_r8
      else if (d .le. 0.25_r8) then
         vfallstp = 1812.1_r8 * d**0.638_r8
      else
         vfallstp = 1069.8_r8 * d**0.235_r8
      end if

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(i) = vfall
      precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
   end do
   precipsum = precipsum*pi*1.333333_r8

   rnumsum = 0._r8
   do i = 1, nr
      xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
      rnumsum = rnumsum + xnumrainsv(i)
   end do

   !
   !   compute aerosol concentrations
   !	aaerosv = particle radius (cm)
   !	fnumaerosv = fraction of total number in the bin (--)
   !	fvolaerosv = fraction of total volume in the bin (--)
   !
   anumsum = 0._r8
   avolsum = 0._r8
   do i = 1, na
      x = xlo + (i-1)*dx
      a = exp( x )
      aaerosv(i) = a
      dum = (x - xg0)/sx
      ynumaerosv(i) = exp( -0.5_r8*dum*dum )
      yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
      anumsum = anumsum + ynumaerosv(i)
      avolsum = avolsum + yvolaerosv(i)
   end do

   do i = 1, na
      ynumaerosv(i) = ynumaerosv(i)/anumsum
      yvolaerosv(i) = yvolaerosv(i)/avolsum
   end do


   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !
   !   outer loop for rain drop radius
   !
   jr_loop: do jr = 1, nr

      r = rrainsv(jr)
      vfall = vfallrainsv(jr)

      reynolds = r * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      !
      !   inner loop for aerosol particle radius
      !
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8

      ja_loop: do ja = 1, na

         a = aaerosv(ja)

         chi = a/r

         dum = freepath/a
         dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
         taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)

         aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
         aerodiffus = boltz_cgs*temp*taurelax/aeromass

         schmidt = airkinvisc/aerodiffus
         stokes = vfall*taurelax/r

         ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
              (reynolds*schmidt)

         dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
              (1._r8 + xmuwaterair/sqrtreynolds)
         eintercept = 4._r8*chi*(chi + dum)

         dum = log( 1._r8 + reynolds )
         sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
         eimpact = 0._r8
         if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
         end if

         etotal = ebrown + eintercept + eimpact
         etotal = min( etotal, 1.0_r8 )

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
 end subroutine calc_1_impact_rate
  
  !=============================================================================
  !=============================================================================
  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), target, intent(in)    :: mbar(ncol,pver)
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    interface
       subroutine qqcw2vmr_codon(ncol_c, pver_c, fldcw_ld1_c, mbar_p, fldcw_p, adv_mass_c, vmr_p) &
            bind(c, name="qqcw2vmr_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, fldcw_ld1_c
         type(c_ptr), value :: mbar_p, fldcw_p, vmr_p
         real(c_double), value :: adv_mass_c
       end subroutine qqcw2vmr_codon
    end interface

    call qqcw2vmr_select_impl()

    if (qqcw2vmr_use_native_impl) then
       do m=1,gas_pcnst
          if( adv_mass(m) /= 0._r8 ) then
             fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
             if(associated(fldcw)) then
                do k=1,pver
                   vmr(:ncol,k,m) = mbar(:ncol,k) * fldcw(:ncol,k) / adv_mass(m)
                end do
             else
                vmr(:,:,m) = 0.0_r8
             end if
          end if
       end do
       return
    end if

    do m=1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             call qqcw2vmr_codon( &
                  int(ncol, c_int64_t), int(pver, c_int64_t), int(size(fldcw,1), c_int64_t), &
                  c_loc(mbar(1,1)), c_loc(fldcw(1,1)), &
                  real(adv_mass(m), c_double), c_loc(vmr(1,1,m)) &
             )
          else
             vmr(:,:,m) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr


  !=============================================================================
  !=============================================================================
  subroutine vmr2qqcw( lchnk, vmr, mbar, ncol, im, pbuf )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use chem_mods,       only : adv_mass, gas_pcnst
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer,  only : physics_buffer_desc
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), target, intent(in)    :: mbar(ncol,pver)
    real(r8), target, intent(in)    :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    interface
       subroutine vmr2qqcw_codon(ncol_c, pver_c, fldcw_ld1_c, vmr_p, mbar_p, adv_mass_c, fldcw_p) &
            bind(c, name="vmr2qqcw_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, fldcw_ld1_c
         type(c_ptr), value :: vmr_p, mbar_p, fldcw_p
         real(c_double), value :: adv_mass_c
       end subroutine vmr2qqcw_codon
    end interface

    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    call vmr2qqcw_select_impl()

    if (vmr2qqcw_use_native_impl) then
       do m = 1,gas_pcnst
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
             do k = 1,pver
                fldcw(:ncol,k) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
             end do
          end if
       end do
       return
    end if

    do m = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
       if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
          call vmr2qqcw_codon( &
               int(ncol, c_int64_t), int(pver, c_int64_t), int(size(fldcw,1), c_int64_t), &
               c_loc(vmr(1,1,m)), c_loc(mbar(1,1)), real(adv_mass(m), c_double), c_loc(fldcw(1,1)) &
          )
       end if
    end do

  end subroutine vmr2qqcw

end module aero_model

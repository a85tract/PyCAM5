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
  use physconst,      only: gravit, mwdry, rair, rhoh2o, tmelt, pi, boltz
  use spmd_utils,     only: masterproc
  use ap_aero_model_emissions_scheme, only: aero_model_emissions_run
  use ap_aero_model_gasaerexch_scheme, only: aero_model_gasaerexch_run
  use ap_aero_model_drydep_scheme, only: aero_model_drydep_run
  use ap_modal_aerosol_timer_hooks, only: register_modal_aerosol_timer_hooks

  use cam_history,    only: outfld, fieldname_len
  use chem_mods,      only: gas_pcnst, adv_mass
  use mo_tracname,    only: solsym

  use modal_aero_data,only: cnst_name_cw, qqcw_get_field
  use modal_aero_data,only: ntot_amode, modename_amode, maxd_aspectype, &
       alnsg_amode, sigmag_amode, nspec_amode, numptr_amode, &
       numptrcw_amode, lmassptr_amode, lmassptrcw_amode
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

contains
  
  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    character(len=16) :: aer_wetdep_list(pcnst) = ' '
    character(len=16) :: aer_drydep_list(pcnst) = ' '

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

    wetdep_list = aer_wetdep_list
    drydep_list = aer_drydep_list

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register
    use modal_aero_initialize_data, only : modal_aero_register

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

    call register_modal_aerosol_timer_hooks(t_startf, t_stopf)

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
  subroutine aero_model_drydep(state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend)

    use mo_drydep, only : n_land_type, fraction_landuse
    use modal_aero_deposition, only : set_srf_drydep

    type(physics_state), intent(in) :: state
    real(r8), intent(in) :: obklen(:)
    real(r8), intent(in) :: ustar(:)
    type(cam_in_t), target, intent(in) :: cam_in
    real(r8), intent(in) :: dt
    type(cam_out_t), intent(inout) :: cam_out
    type(physics_ptend), intent(out) :: ptend
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: l, m, mm
    integer :: lchnk, ncol
    integer :: errflg
    logical :: qcw_needed(pcnst)
    character(len=512) :: errmsg
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)
    real(r8), pointer :: fldcw(:,:)
    real(r8), allocatable :: qcw(:,:,:)
    real(r8), allocatable :: fv(:), ram1(:)
    real(r8), allocatable :: ddv(:,:,:)
    real(r8), allocatable :: dep_flux_is(:,:), dep_trb_is(:,:), dep_grv_is(:,:)
    real(r8), allocatable :: dep_flux_cw(:,:), dep_trb_cw(:,:), dep_grv_cw(:,:)
    real(r8), allocatable :: aerdepdryis(:,:), aerdepdrycw(:,:)
    logical, allocatable :: active_is(:), active_cw(:)

    call t_startf('ap_aero_model_drydep_run')

    lchnk = state%lchnk
    ncol = state%ncol
    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', &
         lq=drydep_lq)

    call pbuf_get_field(pbuf, dgnumwet_idx, dgncur_awet, &
         start=(/1,1,1/), kount=(/pcols,pver,nmodes/))
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens, &
         start=(/1,1,1/), kount=(/pcols,pver,nmodes/))

    allocate(qcw(pcols,pver,pcnst), fv(pcols), ram1(pcols), &
         ddv(pcols,pver,pcnst), dep_flux_is(pcols,pcnst), &
         dep_trb_is(pcols,pcnst), dep_grv_is(pcols,pcnst), &
         dep_flux_cw(pcols,pcnst), dep_trb_cw(pcols,pcnst), &
         dep_grv_cw(pcols,pcnst), aerdepdryis(pcols,pcnst), &
         aerdepdrycw(pcols,pcnst), active_is(pcnst), active_cw(pcnst))

    qcw(:,:,:) = 0._r8
    qcw_needed(:) = .false.
    do m = 1, ntot_amode
       mm = numptrcw_amode(m)
       if (mm > 0) qcw_needed(mm) = .true.
       do l = 1, nspec_amode(m)
          mm = lmassptrcw_amode(l,m)
          if (mm > 0) qcw_needed(mm) = .true.
       end do
    end do
    do mm = 1, pcnst
       if (qcw_needed(mm)) then
          fldcw => qqcw_get_field(pbuf, mm, lchnk)
          qcw(:,:,mm) = fldcw(:,:)
       end if
    end do

    call aero_model_drydep_run( &
         pcols, pver, pverp, pcnst, ncol, ntot_amode, maxd_aspectype, &
         n_land_type, dt, pi, boltz, gravit, rair, rhoh2o, &
         cam_in%landfrac, cam_in%icefrac, cam_in%ocnfrac, obklen, ustar, &
         cam_in%ram1, cam_in%fv, fraction_landuse(:,:,lchnk), state%t, &
         state%pmid, state%pdel, state%pint, state%q, dgncur_awet, wetdens, &
         alnsg_amode, sigmag_amode, nspec_amode, numptr_amode, &
         numptrcw_amode, lmassptr_amode, lmassptrcw_amode, qcw, ptend%lq, &
         ptend%q, fv, ram1, ddv, dep_flux_is, dep_trb_is, dep_grv_is, &
         dep_flux_cw, dep_trb_cw, dep_grv_cw, aerdepdryis, aerdepdrycw, &
         active_is, active_cw, errmsg, errflg)

    if (errflg /= 0) then
       write(iulog,*) trim(errmsg)
       call endrun('aero_model_drydep: standalone kernel failure')
    end if

    do mm = 1, pcnst
       if (qcw_needed(mm)) then
          fldcw => qqcw_get_field(pbuf, mm, lchnk)
          fldcw(1:ncol,:) = qcw(1:ncol,:,mm)
       end if
    end do

    call outfld('airFV', fv, pcols, lchnk)
    call outfld('RAM1', ram1, pcols, lchnk)
    do mm = 1, pcnst
       if (active_is(mm)) then
          call outfld(trim(cnst_name(mm))//'DDV', ddv(:,:,mm), pcols, lchnk)
          call outfld(trim(cnst_name(mm))//'DDF', dep_flux_is(:,mm), pcols, lchnk)
          call outfld(trim(cnst_name(mm))//'TBF', dep_trb_is(:,mm), pcols, lchnk)
          call outfld(trim(cnst_name(mm))//'GVF', dep_grv_is(:,mm), pcols, lchnk)
          call outfld(trim(cnst_name(mm))//'DTQ', ptend%q(:,:,mm), pcols, lchnk)
       end if
       if (active_cw(mm)) then
          call outfld(trim(cnst_name_cw(mm))//'DDF', dep_flux_cw(:,mm), pcols, lchnk)
          call outfld(trim(cnst_name_cw(mm))//'TBF', dep_trb_cw(:,mm), pcols, lchnk)
          call outfld(trim(cnst_name_cw(mm))//'GVF', dep_grv_cw(:,mm), pcols, lchnk)
       end if
    end do

    if (.not. aerodep_flx_prescribed()) then
       call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    end if

    deallocate(qcw, fv, ram1, ddv, dep_flux_is, dep_trb_is, dep_grv_is, &
         dep_flux_cw, dep_trb_cw, dep_grv_cw, aerdepdryis, aerdepdrycw, &
         active_is, active_cw)

    call t_stopf('ap_aero_model_drydep_run')

  end subroutine aero_model_drydep


  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use ap_aero_model_wetdep_scheme, only : aero_model_wetdep_run


    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)


    call t_startf('ap_aero_model_wetdep_run')
    call aero_model_wetdep_run(state, dt, dlf, cam_out, ptend, pbuf, &
         nmodes, dgnumwet_idx, qaerwat_idx, fracis_idx, &
         sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial, &
         nwetdep, wetdep_lq, dlndg_nimptblgrow, scavimptblnum, scavimptblvol)
    call t_stopf('ap_aero_model_wetdep_run')

  end subroutine aero_model_wetdep

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

    strato_sad = 0._r8

    if (.not.modal_strat_sulfate) return

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

    use modal_aero_data, only : maxd_aspectype, ntot_aspectype, &
         modeptr_accum, modeptr_aitken, modeptr_pcarbon,        &
         numptr_amode, mprognum_amode, nspec_amode,            &
         lmassptr_amode, lspectype_amode,                       &
         lptr_so4_a_amode, lptr_nh4_a_amode,                   &
         lptr_soa_a_amode, lptr_pom_a_amode,                   &
         sigmag_amode, alnsg_amode, specmw_amode,              &
         specdens_amode, specmw_so4_amode,                     &
         specmw_nh4_amode, specmw_soa_amode,                   &
         specdens_so4_amode, specdens_nh4_amode,               &
         specdens_soa_amode, dgnumlo_amode, dgnum_amode,       &
         dgnumhi_amode, qqcw_get_field
    use modal_aero_gasaerexch, only : maxspec_pcage,            &
         modefrm_pcage, modetoo_pcage, nspecfrm_pcage,         &
         lspecfrm_pcage, lspectoo_pcage,                        &
         n_so4_monolayers_pcage, dr_so4_monolayers_pcage,     &
         soa_equivso4_factor
    use modal_aero_coag, only : pair_option_acoag,              &
         maxpair_acoag, maxspec_acoag, npair_acoag,            &
         modefrm_acoag, modetoo_acoag, nspecfrm_acoag,         &
         lspecfrm_acoag, lspectoo_acoag
    use modal_aero_newnuc, only : l_h2so4_sv, l_nh3_sv,        &
         lnumait_sv, lnh4ait_sv, lso4ait_sv
    use mo_setsox, only : setsox, has_sox
    use time_manager, only : get_nstep
    use wv_saturation, only : qsat

    integer,  intent(in) :: loffset
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk
    integer,  intent(in) :: troplev(pcols)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: reaction_rates(:,:,:)
    real(r8), intent(in) :: tfld(:,:)
    real(r8), intent(in) :: pmid(:,:)
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(in) :: mbar(:,:)
    real(r8), intent(in) :: relhum(:,:)
    real(r8), intent(in) :: airdens(:,:)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: vmr0(:,:,:)
    real(r8), intent(inout) :: vmr(:,:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: i, k, l, lb, lmz, m, n, nstep
    integer :: l_so4g, l_nh4g, l_msag, l_soag
    integer :: mait, lptr_nh4_aitken
    integer :: ga_ierr, rename_ierr, co_ierr
    real(r8) :: dgnumlo_aitken, dgnum_aitken, dgnumhi_aitken
    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:), sulfeq_ptr(:,:,:)
    real(r8), pointer :: pblh(:), fldcw(:,:)
    real(r8), target :: sulfeq_zero(pcols,pver,ntot_amode)
    real(r8) :: ev_sat(pcols,pver), qv_sat(pcols,pver)
    real(r8) :: vmr_before_setsox(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw_before_setsox(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)
    real(r8) :: gs_flux(pcols,gas_pcnst)
    real(r8) :: aq_flux(pcols,gas_pcnst)
    real(r8) :: ga_qsrflx(pcols,gas_pcnst,2)
    real(r8) :: ga_qqcwsrflx(pcols,gas_pcnst,2)
    real(r8) :: nn_qsrflx(pcols,pcnst)
    real(r8) :: co_qsrflx(pcols,pcnst)
    logical :: ga_dotend(gas_pcnst), ga_dotendqqcw(gas_pcnst)
    logical :: ga_dotendrn(gas_pcnst), ga_dotendqqcwrn(gas_pcnst)
    logical :: nn_dotend(pcnst), co_dotend(pcnst)
    logical :: has_sulfeq_local
    character(len=fieldname_len+3) :: fieldname
    character(len=32) :: name

    call pbuf_get_field(pbuf, dgnum_idx, dgnum)
    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet)
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
    call pbuf_get_field(pbuf, pblh_idx, pblh)

    do n = 1, ntot_amode
       call outfld(dgnum_name(n), dgnum(1:ncol,1:pver,n), ncol, lchnk)
       call outfld(dgnumwet_name(n), dgnumwet(1:ncol,1:pver,n), ncol, lchnk)
    end do

    nstep = get_nstep()

    do m = 1, gas_pcnst
       if (adv_mass(m) /= 0.0_r8) then
          fldcw => qqcw_get_field(pbuf, m+loffset, lchnk, errorhandle=.true.)
          if (associated(fldcw)) then
             do k = 1, pver
                vmrcw(:,k,m) = mbar(1:ncol,k) * fldcw(1:ncol,k) / adv_mass(m)
             end do
          else
             vmrcw(:,:,m) = 0.0_r8
          end if
       end if
    end do

    vmr_before_setsox = vmr(1:ncol,1:pver,1:gas_pcnst)
    vmrcw_before_setsox = vmrcw
    if (has_sox) then
       call setsox(ncol, lchnk, loffset, delt, pmid, pdel, tfld, &
            mbar, cwat, cldfr, cldnum, airdens, invariants,    &
            vmrcw, vmr)
    end if

    call qsat(tfld(1:ncol,1:pver), pmid(1:ncol,1:pver),        &
         ev_sat(1:ncol,1:pver), qv_sat(1:ncol,1:pver))

    call cnst_get_ind('H2SO4', l_so4g, .false.)
    call cnst_get_ind('NH3', l_nh4g, .false.)
    call cnst_get_ind('MSA', l_msag, .false.)
    call cnst_get_ind('SOAG', l_soag, .false.)
    l_so4g = l_so4g - loffset
    l_nh4g = l_nh4g - loffset
    l_msag = l_msag - loffset
    l_soag = l_soag - loffset

    mait = modeptr_aitken
    lptr_nh4_aitken = 0
    dgnumlo_aitken = 1.0_r8
    dgnum_aitken = 1.0_r8
    dgnumhi_aitken = 1.0_r8
    if (mait > 0 .and. mait <= ntot_amode) then
       lptr_nh4_aitken = lptr_nh4_a_amode(mait)
       dgnumlo_aitken = dgnumlo_amode(mait)
       dgnum_aitken = dgnum_amode(mait)
       dgnumhi_aitken = dgnumhi_amode(mait)
    end if

    has_sulfeq_local = sulfeq_idx > 0
    if (has_sulfeq_local) then
       call pbuf_get_field(pbuf, sulfeq_idx, sulfeq_ptr)
    else
       sulfeq_zero = 0.0_r8
       sulfeq_ptr => sulfeq_zero
    end if

    call t_startf('ap_aero_model_gasaerexch_run')
    call t_startf('modal_gas-aer_exchng')
    call t_startf('modal_nucl')
    call t_startf('modal_coag')
    call aero_model_gasaerexch_run(                              &
         lchnk, ncol, nstep, pcols, pver, gas_pcnst, pcnst,     &
         ntot_amode, top_lev, maxd_aspectype, ntot_aspectype,   &
         maxpair_acoag, loffset, delt,                          &
         l_so4g, l_nh4g, l_msag, l_soag,                       &
         modefrm_pcage, modetoo_pcage, nspecfrm_pcage,         &
         lspecfrm_pcage, lspectoo_pcage,                       &
         modeptr_accum, modeptr_aitken, modeptr_pcarbon,       &
         numptr_amode, mprognum_amode, nspec_amode,            &
         lmassptr_amode, lspectype_amode,                      &
         lptr_so4_a_amode, lptr_nh4_a_amode,                  &
         lptr_soa_a_amode, lptr_pom_a_amode,                  &
         sigmag_amode, alnsg_amode, specmw_amode,             &
         specdens_amode, specmw_so4_amode,                    &
         specmw_nh4_amode, specmw_soa_amode,                  &
         specdens_so4_amode, specdens_nh4_amode,              &
         specdens_soa_amode, dr_so4_monolayers_pcage,         &
         n_so4_monolayers_pcage, soa_equivso4_factor,         &
         pair_option_acoag, npair_acoag,                      &
         modefrm_acoag, modetoo_acoag, nspecfrm_acoag,        &
         lspecfrm_acoag, lspectoo_acoag,                      &
         l_h2so4_sv, l_nh3_sv, lnumait_sv, lnh4ait_sv,       &
         lso4ait_sv, lptr_nh4_aitken,                         &
         dgnumlo_aitken, dgnum_aitken, dgnumhi_aitken,       &
         gravit, mwdry, rair, tmelt, adv_mass,               &
         tfld, pmid, pdel, mbar, zm, pblh, qh2o, qv_sat,     &
         cldfr(1:ncol,1:pver), troplev,                       &
         vmr0(1:ncol,1:pver,1:gas_pcnst), vmr_before_setsox, &
         vmrcw_before_setsox, vmr, vmrcw,                    &
         del_h2so4_gasprod(1:ncol,1:pver), dgnum, dgnumwet,  &
         wetdens, has_sulfeq_local, sulfeq_ptr,              &
         gs_flux, aq_flux, ga_qsrflx, ga_qqcwsrflx,          &
         ga_dotend, ga_dotendqqcw, ga_dotendrn,              &
         ga_dotendqqcwrn, nn_qsrflx, nn_dotend,              &
         co_qsrflx, co_dotend, ga_ierr, rename_ierr, co_ierr)
    call t_stopf('modal_coag')
    call t_stopf('modal_nucl')
    call t_stopf('modal_gas-aer_exchng')
    call t_stopf('ap_aero_model_gasaerexch_run')

    if (ga_ierr /= 0) then
       write(*,'(/a/a,2i7)') &
            '*** aero_model_gasaerexch -- cannot find H2SO4 species', &
            '    l_so4g, loffset =', l_so4g, loffset
       call endrun('aero_model_gasaerexch error')
    end if
    if (rename_ierr /= 0) then
       write(*,'(/a,1x,i7)') &
            '*** aero_model_gasaerexch -- renaming error pair =', rename_ierr
       call endrun('aero_model_gasaerexch error')
    end if
    if (co_ierr /= 0) then
       write(*,*) '*** aero_model_gasaerexch coagulation error'
       write(*,*) '    pair_option_acoag =', pair_option_acoag
       call endrun('aero_model_gasaerexch coagulation error')
    end if

    do m = 1, gas_pcnst
       name = 'GS_' // trim(solsym(m))
       call outfld(name, gs_flux(:,m), ncol, lchnk)
       name = 'AQ_' // trim(solsym(m))
       call outfld(name, aq_flux(:,m), ncol, lchnk)
    end do

    do l = 1, gas_pcnst
       lb = l + loffset
       if (ga_dotend(l)) then
          fieldname = trim(cnst_name(lb)) // '_sfgaex1'
          call outfld(fieldname, ga_qsrflx(:,l,1), pcols, lchnk)
       end if
       if (ga_dotendrn(l)) then
          fieldname = trim(cnst_name(lb)) // '_sfgaex2'
          call outfld(fieldname, ga_qsrflx(:,l,2), pcols, lchnk)
       end if
       if (ga_dotendqqcwrn(l)) then
          fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
          call outfld(fieldname, ga_qqcwsrflx(:,l,2), pcols, lchnk)
       end if
    end do

    do l = loffset+1, pcnst
       lmz = l - loffset
       if (nn_dotend(lmz)) then
          do i = 1, ncol
             nn_qsrflx(i,lmz) = nn_qsrflx(i,lmz) * &
                  (adv_mass(lmz) / mwdry)
          end do
          fieldname = trim(cnst_name(l)) // '_sfnnuc1'
          call outfld(fieldname, nn_qsrflx(:,lmz), pcols, lchnk)
       end if
       if (co_dotend(lmz)) then
          co_qsrflx(:,lmz) = co_qsrflx(:,lmz) * &
               (adv_mass(lmz) / (gravit*mwdry))
          fieldname = trim(cnst_name(l)) // '_sfcoag1'
          call outfld(fieldname, co_qsrflx(:,lmz), pcols, lchnk)
       end if
    end do

    do m = 1, gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+loffset, lchnk, errorhandle=.true.)
       if (adv_mass(m) /= 0.0_r8 .and. associated(fldcw)) then
          do k = 1, pver
             fldcw(1:ncol,k) = adv_mass(m) * vmrcw(:,k,m) / mbar(1:ncol,k)
          end do
       end if
    end do

    do n = 1, pcnst
       fldcw => qqcw_get_field(pbuf, n, lchnk, errorhandle=.true.)
       if (associated(fldcw)) then
          call outfld(cnst_name_cw(n), fldcw(:,:), pcols, lchnk)
       end if
    end do

  end subroutine aero_model_gasaerexch


  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_names, seasalt_indices, seasalt_active,seasalt_nbin
    use dust_model,    only: dust_emis, dust_names, dust_indices, dust_active,dust_nbin, dust_nnum
    use physics_types, only: physics_state

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol
    integer :: m, mm
    integer :: errflg
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    real(r8) :: u10cubed(pcols)
    character(len=512) :: errmsg
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then

       call dust_emis( ncol, lchnk, cam_in%dstflx, cam_in%cflx, soil_erod_tmp )

       ! some dust emis diagnostics ...
       sflx(:)=0._r8
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          if (m<=dust_nbin) sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       call t_startf('ap_aero_model_emissions_run')
       call aero_model_emissions_run(ncol, pcols, state%u(:,pver), &
            state%v(:,pver), state%zm(:,pver), z0, u10cubed, errmsg, errflg)
       call t_stopf('ap_aero_model_emissions_run')

       sflx(:)=0._r8

       call seasalt_emis( u10cubed, cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx )

       do m=1,seasalt_nbin
          mm = seasalt_indices(m)
          sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       call outfld('SSTSFMBL',sflx(:),pcols,lchnk)
    endif

  end subroutine aero_model_emissions

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
  !===============================================================================


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


  !=============================================================================
  !=============================================================================

end module aero_model

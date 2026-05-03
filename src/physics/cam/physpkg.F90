module physpkg
!-----------------------------------------------------------------------
! Purpose:
!
! Provides the interface to CAM physics package
!
! Revision history:
! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
!                            initialization of grid info in phys_state.
! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
!
! Mar 2011,   J. Nusbaumer  Added calls to water tracer/isotope routines.
!
!-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
       physics_ptend, physics_tend_init, physics_update,    &
       physics_type_alloc, physics_ptend_dealloc,&
       physics_state_alloc, physics_state_dealloc, physics_tend_alloc, physics_tend_dealloc
  use phys_grid,        only: get_ncols_p
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use constituents,     only: pcnst, cnst_name, cnst_get_ind
  use camsrfexch,       only: cam_out_t, cam_in_t

  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, phys_getopts, waccmx_is
  use scamMod,          only: single_column, scm_crm_mode
  use flux_avg,         only: flux_avg_init
  use infnan,           only: posinf, assignment(=)
#ifdef SPMD
  use mpishorthand
#endif
  use perf_mod
  use cam_logfile,     only: iulog
  use camsrfexch,      only: cam_export

  use modal_aero_calcsize,    only: modal_aero_calcsize_init, modal_aero_calcsize_diag, modal_aero_calcsize_reg
  use modal_aero_wateruptake, only: modal_aero_wateruptake_init, modal_aero_wateruptake_dr, modal_aero_wateruptake_reg

  implicit none
  private
  save

  ! Public methods
  public phys_register ! was initindx  - register physics methods
  public phys_init   ! Public initialization method
  public phys_run1   ! First phase of the public run method
  public phys_run2   ! Second phase of the public run method
  public phys_final  ! Public finalization method

  ! Private module data

  ! Physics package options
  character(len=16) :: shallow_scheme
  character(len=16) :: macrop_scheme
  character(len=16) :: microp_scheme 
  integer           :: cld_macmic_num_steps    ! Number of macro/micro substeps
  logical           :: do_clubb_sgs
  logical           :: use_subcol_microp   ! if true, use subcolumns in microphysics
  logical           :: state_debug_checks  ! Debug physics_state.
  logical           :: clim_modal_aero     ! climate controled by prognostic or prescribed modal aerosols
  logical           :: prog_modal_aero     ! Prognostic modal aerosols present
  logical           :: micro_do_icesupersat
  logical           :: use_native_phys_tstep_impl = .false.
  logical           :: phys_tstep_impl_selected = .false.
  integer           :: phys_tstep_branch_mask = 0
  logical           :: phys_tstep_branch_selected = .false.
  logical           :: use_native_tphysac_flx_net_update_impl = .false.
  logical           :: tphysac_flx_net_update_impl_selected = .false.
  logical           :: use_native_tphysbc_precip_ops_impl = .false.
  logical           :: tphysbc_precip_ops_impl_selected = .false.
  logical           :: use_native_tphysac_t_update_impl = .false.
  logical           :: tphysac_t_update_impl_selected = .false.
  logical           :: use_native_tphysac_q_snapshot_impl = .false.
  logical           :: tphysac_q_snapshot_impl_selected = .false.
  logical           :: use_native_tphysbc_qini_snapshot_impl = .false.
  logical           :: tphysbc_qini_snapshot_impl_selected = .false.
  logical           :: use_native_tphysbc_dadadj_input_impl = .false.
  logical           :: tphysbc_dadadj_input_impl_selected = .false.
  logical           :: use_native_tphysbc_dadadj_output_impl = .false.
  logical           :: tphysbc_dadadj_output_impl_selected = .false.
  logical           :: use_native_tphysbc_dtcore_update_impl = .false.
  logical           :: tphysbc_dtcore_update_impl_selected = .false.
  logical           :: use_native_tphysbc_tini_copy_impl = .false.
  logical           :: tphysbc_tini_copy_impl_selected = .false.
  logical           :: use_native_tphysbc_flx_cnd_sum_impl = .false.
  logical           :: tphysbc_flx_cnd_sum_impl_selected = .false.
  logical           :: use_native_tphysbc_macrop_fluxes_impl = .false.
  logical           :: tphysbc_macrop_fluxes_impl_selected = .false.
  logical           :: use_native_tphysbc_init_fields_impl = .false.
  logical           :: tphysbc_init_fields_impl_selected = .false.
  logical           :: use_native_tphysbc_radheat_flx_net_impl = .false.
  logical           :: tphysbc_radheat_flx_net_impl_selected = .false.
  logical           :: use_native_tphysbc_zero_buffers_impl = .false.
  logical           :: tphysbc_zero_buffers_impl_selected = .false.
  logical           :: use_native_tphysbc_trace_water_clip_impl = .false.
  logical           :: tphysbc_trace_water_clip_impl_selected = .false.
  logical           :: use_native_tphysbc_dadadj_lq_init_impl = .false.
  logical           :: tphysbc_dadadj_lq_init_impl_selected = .false.
  logical           :: use_native_phys_inidat_qpert_expand_impl = .false.
  logical           :: phys_inidat_qpert_expand_impl_selected = .false.
  logical           :: use_native_phys_inidat_qpert_default_impl = .false.
  logical           :: phys_inidat_qpert_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_pblh_default_impl = .false.
  logical           :: phys_inidat_pblh_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_tpert_default_impl = .false.
  logical           :: phys_inidat_tpert_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_cush_default_impl = .false.
  logical           :: phys_inidat_cush_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_tke_default_impl = .false.
  logical           :: phys_inidat_tke_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_kvm_default_impl = .false.
  logical           :: phys_inidat_kvm_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_kvh_default_impl = .false.
  logical           :: phys_inidat_kvh_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_qcwat_default_impl = .false.
  logical           :: phys_inidat_qcwat_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_iccwat_default_impl = .false.
  logical           :: phys_inidat_iccwat_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_lcwat_default_impl = .false.
  logical           :: phys_inidat_lcwat_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_tcwat_default_impl = .false.
  logical           :: phys_inidat_tcwat_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_cloud_default_impl = .false.
  logical           :: phys_inidat_cloud_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_concld_default_impl = .false.
  logical           :: phys_inidat_concld_default_impl_selected = .false.
  logical           :: use_native_phys_inidat_tbot_init_impl = .false.
  logical           :: phys_inidat_tbot_init_impl_selected = .false.
  logical           :: use_native_phys_inidat_batch_impl = .false.
  logical           :: phys_inidat_batch_impl_selected = .false.
  logical           :: phys_inidat_batch_entered_logged = .false.

  !  Physics buffer index
  integer ::  teout_idx          = 0  

  integer ::  tini_idx           = 0
  integer ::  qini_idx           = 0 
  integer ::  cldliqini_idx      = 0 
  integer ::  cldiceini_idx      = 0 

  integer ::  prec_str_idx       = 0
  integer ::  snow_str_idx       = 0
  integer ::  prec_sed_idx       = 0
  integer ::  snow_sed_idx       = 0
  integer ::  prec_pcw_idx       = 0
  integer ::  snow_pcw_idx       = 0
  integer ::  prec_dp_idx        = 0
  integer ::  snow_dp_idx        = 0
  integer ::  prec_sh_idx        = 0
  integer ::  snow_sh_idx        = 0

!======================================================================= 
contains
!======================================================================= 

subroutine phys_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Register constituents and physics buffer fields.
    ! 
    ! Author:    CSM Contact: M. Vertenstein, Aug. 1997
    !            B.A. Boville, Oct 2001
    !            A. Gettelman, Nov 2010 - put micro/macro physics into separate routines
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer,     only: pbuf_init_time
    use physics_buffer,     only: pbuf_add_field, dtype_r8, pbuf_register_subcol
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use spmd_utils,         only: masterproc
    use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name

    use cam_control_mod,    only: moist_physics
    use chemistry,          only: chem_register
    use cloud_fraction,     only: cldfrc_register
    use stratiform,         only: stratiform_register
    use microp_driver,      only: microp_driver_register
    use microp_aero,        only: microp_aero_register
    use macrop_driver,      only: macrop_driver_register
    use clubb_intr,         only: clubb_register_cam
    use conv_water,         only: conv_water_register
    use physconst,          only: mwdry, cpair, mwh2o, cpwv
    use tracers,            only: tracers_register
    use check_energy,       only: check_energy_register
    use carma_intr,         only: carma_register
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
    use ghg_data,           only: ghg_data_register
    use vertical_diffusion, only: vd_register
    use convect_deep,       only: convect_deep_register
    use convect_shallow,    only: convect_shallow_register
    use radiation,          only: radiation_register
    use co2_cycle,          only: co2_register
    use flux_avg,           only: flux_avg_register
    use iondrag,            only: iondrag_register
    use ionosphere,         only: ionos_register
    use string_utils,       only: to_lower
    use prescribed_ozone,   only: prescribed_ozone_register
    use prescribed_volcaero,only: prescribed_volcaero_register
    use prescribed_strataero,only: prescribed_strataero_register
    use prescribed_aero,    only: prescribed_aero_register
    use prescribed_ghg,     only: prescribed_ghg_register
    use sslt_rebin,         only: sslt_rebin_register
    use aoa_tracers,        only: aoa_tracers_register
    use aircraft_emit,      only: aircraft_emit_register
    use cam_diagnostics,    only: diag_register
    use cloud_diagnostics,  only: cloud_diagnostics_register
    use rad_constituents,   only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not
    use subcol,             only: subcol_register
    use subcol_utils,       only: is_subcol_on


   !water isotopes:
    use water_tracer_vars,  only: trace_water
    use water_tracers,      only: wtrc_register

    implicit none
    !---------------------------Local variables-----------------------------
    !
    integer  :: m        ! loop index
    integer  :: mm       ! constituent index 
    integer  :: nmodes
    !-----------------------------------------------------------------------

    ! Get physics options
    call phys_getopts(shallow_scheme_out       = shallow_scheme, &
                      macrop_scheme_out        = macrop_scheme,   &
                      microp_scheme_out        = microp_scheme,   &
                      cld_macmic_num_steps_out = cld_macmic_num_steps, &
                      do_clubb_sgs_out         = do_clubb_sgs,     &
                      use_subcol_microp_out    = use_subcol_microp, &
                      state_debug_checks_out   = state_debug_checks, &
                      micro_do_icesupersat_out = micro_do_icesupersat)

    ! Initialize dyn_time_lvls
    call pbuf_init_time()

    ! Register the subcol scheme
    call subcol_register()

    ! Register water vapor.
    ! ***** N.B. ***** This must be the first call to cnst_add so that
    !                  water vapor is constituent 1.
    if (moist_physics) then
       call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, &
            longname='Specific humidity', readiv=.true., is_convtran1=.true.)
    else
       call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, &
            longname='Specific humidity', readiv=.false., is_convtran1=.true.)
    end if

    ! Fields for physics package diagnostics
    call pbuf_add_field('TINI',      'physpkg', dtype_r8, (/pcols,pver/), tini_idx)
    call pbuf_add_field('QINI',      'physpkg', dtype_r8, (/pcols,pver/), qini_idx)
    call pbuf_add_field('CLDLIQINI', 'physpkg', dtype_r8, (/pcols,pver/), cldliqini_idx)
    call pbuf_add_field('CLDICEINI', 'physpkg', dtype_r8, (/pcols,pver/), cldiceini_idx)

    ! check energy package
    call check_energy_register

    ! If using an ideal/adiabatic physics option, the CAM physics parameterizations 
    ! aren't called.
    if (moist_physics) then

       ! register fluxes for saving across time
       if (phys_do_flux_avg()) call flux_avg_register()

       call cldfrc_register()

       ! cloud water
       if( microp_scheme == 'RK' ) then
          call stratiform_register()
       elseif( microp_scheme == 'MG' ) then
          if (.not. do_clubb_sgs) call macrop_driver_register()
          call microp_aero_register()
          call microp_driver_register()
       end if
       
       ! Register CLUBB_SGS here
       if (do_clubb_sgs) call clubb_register_cam()
       

       call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
       call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
       call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
       call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
       call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
       call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)
       if (is_subcol_on()) then
         call pbuf_register_subcol('PREC_STR', 'phys_register', prec_str_idx)
         call pbuf_register_subcol('SNOW_STR', 'phys_register', snow_str_idx)
         call pbuf_register_subcol('PREC_PCW', 'phys_register', prec_pcw_idx)
         call pbuf_register_subcol('SNOW_PCW', 'phys_register', snow_pcw_idx)
         call pbuf_register_subcol('PREC_SED', 'phys_register', prec_sed_idx)
         call pbuf_register_subcol('SNOW_SED', 'phys_register', snow_sed_idx)
       end if

    ! Who should add FRACIS? 
    ! -- It does not seem that aero_intr should add it since FRACIS is used in convection
    !     even if there are no prognostic aerosols ... so do it here for now 
       call pbuf_add_field('FRACIS','physpkg',dtype_r8,(/pcols,pver,pcnst/),m)

       call conv_water_register()
       
       ! Determine whether its a 'modal' aerosol simulation  or not
       call rad_cnst_get_info(0, nmodes=nmodes)
       clim_modal_aero = (nmodes > 0)

       if (clim_modal_aero) then
          call modal_aero_calcsize_reg()
          call modal_aero_wateruptake_reg()
       endif


       ! water tracers/isotopes
       if(trace_water) then
         call wtrc_register()
       end if

       ! register chemical constituents including aerosols ...
       call chem_register()

       ! co2 constituents
       call co2_register()

       ! register data model ozone with pbuf
       if (cam3_ozone_data_on) then
          call cam3_ozone_data_register()
       end if
       call prescribed_volcaero_register()
       call prescribed_strataero_register()
       call prescribed_ozone_register()
       call prescribed_aero_register()
       call prescribed_ghg_register()
       call sslt_rebin_register

       ! CAM3 prescribed aerosols
       if (cam3_aero_data_on) then
          call cam3_aero_data_register()
       end if

       ! register various data model gasses with pbuf
       call ghg_data_register()

       ! carma microphysics
       ! 
       call carma_register()

       ! Register iondrag variables with pbuf
       call iondrag_register()

       ! Register ionosphere variables with pbuf if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_register()
       endif

       call aircraft_emit_register()

       ! deep convection
       call convect_deep_register

       !  shallow convection
       call convect_shallow_register

       ! radiation
       call radiation_register
       call cloud_diagnostics_register

       ! vertical diffusion
       if (.not. do_clubb_sgs) call vd_register()
    end if

    ! Register diagnostics PBUF
    call diag_register()

    ! Register age of air tracers
    call aoa_tracers_register()

    ! Register test tracers
    ! ***** N.B. ***** This is the last call to register constituents because
    !                  the test tracers fill the remaining available slots up
    !                  to constituent number PCNST -- regardless of what PCNST is set to.
    call tracers_register()

    ! All tracers registered, check that the dimensions are correct
    call cnst_chk_dim()

    ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

end subroutine phys_register



  !======================================================================= 

subroutine phys_inidat( cam_out, pbuf2d )
    use cam_abortutils,      only: endrun

    use physics_buffer,      only: pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, dyn_time_lvls


    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
    use dycore,              only: dycore_is
    use polar_avg,           only: polar_average
    use short_lived_species, only: initialize_short_lived_species
    use comsrf,              only: landm, sgh, sgh30
    use cam_control_mod,     only: aqua_planet

    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer :: lchnk, m, n, i, k, ncol
    type(file_desc_t), pointer :: fh_ini, fh_topo
    character(len=8) :: fieldname
    real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
    real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
    real(r8), pointer :: qpert(:,:)
    real(r8) :: posinf_r8

    character*11 :: subname='phys_inidat' ! subroutine name
    integer :: tpert_idx, qpert_idx, pblh_idx

    logical :: found=.false., found2=.false., found_primary=.false.
    integer :: ierr, qcwat_source, iccwat_source, lcwat_source, tcwat_source
    character(len=4) :: dim1name
    integer :: ixcldice, ixcldliq
    nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

    fh_ini=>initial_file_get_id()

    !   dynamics variables are handled in dyn_init - here we read variables needed for physics 
    !   but not dynamics

    if(dycore_is('UNSTRUCTURED')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if
    if(aqua_planet) then
       sgh = 0._r8
       sgh30 = 0._r8
       landm = 0._r8
    else
       fh_topo=>topo_file_get_id()
       call infld('SGH', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh, found, grid_map='PHYS')
       if(.not. found) call endrun('ERROR: SGH not found on topo file')

       call infld('SGH30', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh30, found, grid_map='PHYS')
       if(.not. found) then
          if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
          if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
          sgh30 = sgh
       end if

       call infld('LANDM_COSLAT', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            landm, found, grid_map='PHYS')

       if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')
    end if

    allocate(tptr(1:pcols,begchunk:endchunk))

    call infld('PBLH', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    call phys_inidat_pblh_default(pcols, endchunk-begchunk+1, found, tptr)
    if(.not. found) then
       if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
    end if
    pblh_idx = pbuf_get_index('pblh')

    call pbuf_set_field(pbuf2d, pblh_idx, tptr)

    call infld('TPERT', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    call phys_inidat_tpert_default(pcols, endchunk-begchunk+1, found, tptr)
    if(.not. found) then
       if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
    end if
    tpert_idx = pbuf_get_index( 'tpert')
    call pbuf_set_field(pbuf2d, tpert_idx, tptr)

    fieldname='QPERT'  
    qpert_idx = pbuf_get_index( 'qpert',ierr)
    if (qpert_idx > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            tptr, found, grid_map='PHYS')
       call phys_inidat_qpert_default(pcols, endchunk-begchunk+1, found, tptr)
       if(.not. found) then
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
       end if

       allocate(tptr3d_2(pcols,pcnst,begchunk:endchunk))
call phys_inidat_qpert_expand(pcols, pcnst, endchunk-begchunk+1, tptr, tptr3d_2)

       call pbuf_set_field(pbuf2d, qpert_idx, tptr3d_2)
       deallocate(tptr3d_2)
    end if

    fieldname='CUSH'
    m = pbuf_get_index('cush')
    call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr, found, grid_map='PHYS')
    call phys_inidat_cush_default(pcols, endchunk-begchunk+1, found, tptr, 1000._r8)
    if(.not.found) then
       if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
    end if
    do n=1,dyn_time_lvls
       call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
    end do
    deallocate(tptr)

    do lchnk=begchunk,endchunk
       posinf_r8 = posinf
       call phys_inidat_tbot_init(pcols, cam_out(lchnk)%tbot, posinf_r8)
    end do

    !
    ! 3-D fields
    !

    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname='CLOUD'
    m = pbuf_get_index('CLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='PHYS')
    call phys_inidat_cloud_default(pcols, pver, endchunk-begchunk+1, found, tptr3d, 0._r8)
    do n = 1, dyn_time_lvls
       call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
    if(.not. found) then
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname='QCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='PHYS')
       found_primary = found
       found2 = .false.
       if(.not. found_primary) then
          call infld('Q',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found2, grid_map='PHYS')
       end if
       call phys_inidat_qcwat_default(found_primary, found2, qcwat_source)
       if (qcwat_source == 2) then
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
          if(dycore_is('LR')) call polar_average(pver, tptr3d)
       else if (qcwat_source == 0) then
          call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    fieldname = 'ICCWAT'
    m = pbuf_get_index(fieldname, ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, grid_map='phys')
       found2 = .false.
       if(.not. found) then
          call cnst_get_ind('CLDICE', ixcldice)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
             tptr3d, found2, grid_map='PHYS')
       end if
       call phys_inidat_iccwat_default(found, found2, iccwat_source)
       if(iccwat_source /= 0) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call pbuf_set_field(pbuf2d, m, 0._r8)
       end if
       if (masterproc .and. iccwat_source /= 1) then
          if (iccwat_source == 2) then
             write(iulog,*) trim(fieldname), ' initialized with CLDICE'
          else
             write(iulog,*) trim(fieldname), ' initialized to 0.0'
          end if
       end if
    end if

    fieldname = 'LCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       found_primary = found
       found2 = .false.
       if(.not. found_primary) then
          allocate(tptr3d_2(pcols,pver,begchunk:endchunk))
          call cnst_get_ind('CLDICE', ixcldice)
          call cnst_get_ind('CLDLIQ', ixcldliq)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='PHYS')
          call infld('CLDLIQ',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d_2, found2, grid_map='PHYS')
       end if

       call phys_inidat_lcwat_default(found_primary, found, found2, lcwat_source)

       if(lcwat_source == 1) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          if(lcwat_source == 2) then
             tptr3d(:,:,:)=tptr3d(:,:,:)+tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
          else if (lcwat_source == 3) then
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE only'
          else if (lcwat_source == 4) then
             tptr3d(:,:,:)=tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDLIQ only'
          end if

          if (lcwat_source /= 0) then
             do n = 1, dyn_time_lvls
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
             if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
          end if
       end if

       if (.not. found_primary) then
          deallocate(tptr3d_2)
       end if
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'TCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       found_primary = found
       if(.not.found_primary) then
          call infld('T', fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='phys')
       end if
       call phys_inidat_tcwat_default(found_primary, tcwat_source)
       if(tcwat_source == 2) then
          if(dycore_is('LR')) call polar_average(pver, tptr3d)
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pverp,begchunk:endchunk))

    fieldname = 'TKE'
    m = pbuf_get_index( 'tke')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    call phys_inidat_tke_default(pcols, pverp, endchunk-begchunk+1, found, tptr3d, 0.01_r8)
    call pbuf_set_field(pbuf2d, m, tptr3d)
    if (.not. found) then
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
    end if


    fieldname = 'KVM'
    m = pbuf_get_index('kvm')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    call phys_inidat_kvm_default(pcols, pverp, endchunk-begchunk+1, found, tptr3d, 0._r8)
    call pbuf_set_field(pbuf2d, m, tptr3d)
    if (.not. found) then
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if


    fieldname = 'KVH'
    m = pbuf_get_index('kvh')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    call phys_inidat_kvh_default(pcols, pverp, endchunk-begchunk+1, found, tptr3d, 0._r8)
    call pbuf_set_field(pbuf2d, m, tptr3d)
    if (.not. found) then
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'CONCLD'
    m = pbuf_get_index('CONCLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    call phys_inidat_concld_default(pcols, pver, endchunk-begchunk+1, found, tptr3d, 0._r8)
    do n = 1, dyn_time_lvls
       call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
    if(.not. found) then
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate (tptr3d)

    call initialize_short_lived_species(fh_ini, pbuf2d)
end subroutine phys_inidat


subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )

    !----------------------------------------------------------------------- 
    ! 
    ! Initialization of physics package.
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
    use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                  karman, rhodair, physconst_init 
    use ref_pres,           only: pref_edge, pref_mid

    use carma_intr,         only: carma_init
    use cloud_rad_props,    only: cloud_rad_props_init
    use cam_control_mod,    only: nsrest  ! restart flag
    use check_energy,       only: check_energy_init
    use chemistry,          only: chem_init
    use prescribed_ozone,   only: prescribed_ozone_init
    use prescribed_ghg,     only: prescribed_ghg_init
    use prescribed_aero,    only: prescribed_aero_init
    use aerodep_flx,        only: aerodep_flx_init
    use aircraft_emit,      only: aircraft_emit_init
    use prescribed_volcaero,only: prescribed_volcaero_init
    use prescribed_strataero,only: prescribed_strataero_init
    use cloud_fraction,     only: cldfrc_init
    use cldfrc2m,           only: cldfrc2m_init
    use co2_cycle,          only: co2_init, co2_transport
    use convect_deep,       only: convect_deep_init
    use convect_shallow,    only: convect_shallow_init
    use cam_diagnostics,    only: diag_init
    use gw_drag,            only: gw_init
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_init
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_init
    use radheat,            only: radheat_init
    use radiation,          only: radiation_init
    use cloud_diagnostics,  only: cloud_diagnostics_init
    use stratiform,         only: stratiform_init
    use wv_saturation,      only: wv_sat_init
    use microp_driver,      only: microp_driver_init
    use microp_aero,        only: microp_aero_init
    use macrop_driver,      only: macrop_driver_init
    use conv_water,         only: conv_water_init
    use tracers,            only: tracers_init
    use aoa_tracers,        only: aoa_tracers_init
    use rayleigh_friction,  only: rayleigh_friction_init
    use pbl_utils,          only: pbl_utils_init
    use vertical_diffusion, only: vertical_diffusion_init
    use dycore,             only: dycore_is
    use phys_debug_util,    only: phys_debug_init
    use phys_debug,         only: phys_debug_state_init
    use rad_constituents,   only: rad_cnst_init
    use aer_rad_props,      only: aer_rad_props_init
    use subcol,             only: subcol_init
    use qbo,                only: qbo_init
    use iondrag,            only: iondrag_init
#if ( defined OFFLINE_DYN )
    use metdata,            only: metdata_phys_init
#endif
    use ionosphere,         only: ionos_init  ! Initialization of ionosphere module (WACCM-X)
    use majorsp_diffusion,  only: mspd_init   ! Initialization of major species diffusion module (WACCM-X)
    use clubb_intr,         only: clubb_ini_cam
    use sslt_rebin,         only: sslt_rebin_init
    use tropopause,         only: tropopause_init
    use solar_data,         only: solar_data_init
    use rad_solar_var,      only: rad_solar_var_init

    !water isotopes:   
     use water_tracer_vars, only: trace_water
     use water_tracers,     only: wtrc_init

    ! Input/output arguments
    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

    ! local variables
    integer :: lchnk
    !-----------------------------------------------------------------------

    call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

    do lchnk = begchunk, endchunk
       call physics_state_set_grid(lchnk, phys_state(lchnk))
    end do

    !-------------------------------------------------------------------------------------------
    ! Initialize any variables in physconst which are not temporally and/or spatially constant
    !------------------------------------------------------------------------------------------- 
    call physconst_init()

    ! Initialize debugging a physics column
    call phys_debug_init()

    call pbuf_initialize(pbuf2d)

    ! Initialize subcol scheme
    call subcol_init(pbuf2d)

    ! diag_init makes addfld calls for dynamics fields that are output from
    ! the physics decomposition
    call diag_init()

    call check_energy_init()

    call tracers_init()

    ! age of air tracers
    call aoa_tracers_init()

    teout_idx = pbuf_get_index( 'TEOUT')

    ! For adiabatic or ideal physics don't need to initialize any of the
    ! parameterizations below:
    if (adiabatic .or. ideal_phys) return

    if (nsrest .eq. 0) then
       call phys_inidat(cam_out, pbuf2d) 
    end if
    
    ! wv_saturation is relatively independent of everything else and
    ! low level, so init it early. Must at least do this before radiation.
    call wv_sat_init

    ! CAM3 prescribed aerosols
    if (cam3_aero_data_on) call cam3_aero_data_init(phys_state)

    ! Initialize rad constituents and their properties
    call rad_cnst_init()
    call aer_rad_props_init()
    call cloud_rad_props_init()

    ! initialize carma
    call carma_init()

    ! solar irradiance data modules
    call solar_data_init()

    ! Prognostic chemistry.
    call chem_init(phys_state,pbuf2d)

    ! Prescribed tracers
    call prescribed_ozone_init()
    call prescribed_ghg_init()
    call prescribed_aero_init()
    call aerodep_flx_init()
    call aircraft_emit_init()
    call prescribed_volcaero_init()
    call prescribed_strataero_init()

    ! co2 cycle            
    if (co2_transport()) then
       call co2_init()
    end if

    ! CAM3 prescribed ozone
    if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

    call gw_init()

    call rayleigh_friction_init()

    call pbl_utils_init(gravit, karman, cpair, rair, zvir)
    if (.not. do_clubb_sgs) call vertical_diffusion_init(pbuf2d)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_init ()
       ! Initialization of ionosphere module if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_init()
       endif
    endif

    call tsinti(tmelt, latvap, rair, stebol, latice)

    call radiation_init

    call rad_solar_var_init()

    call cloud_diagnostics_init()

    call radheat_init(pref_mid)

    call convect_shallow_init(pref_edge, pbuf2d)

    call cldfrc_init()
    call cldfrc2m_init()

    call convect_deep_init(pref_edge)

    if( microp_scheme == 'RK' ) then
       call stratiform_init()
    elseif( microp_scheme == 'MG' ) then 
       if (.not. do_clubb_sgs) call macrop_driver_init(pbuf2d)
       call microp_aero_init()
       call microp_driver_init(pbuf2d)
       call conv_water_init
    end if


    ! initiate CLUBB within CAM
    if (do_clubb_sgs) call clubb_ini_cam(pbuf2d)

    call qbo_init

    call iondrag_init(pref_mid)

#if ( defined OFFLINE_DYN )
    call metdata_phys_init()
#endif
    call sslt_rebin_init()
    call tropopause_init()

   !Water isotopes:
    call wtrc_init

    prec_dp_idx  = pbuf_get_index('PREC_DP')
    snow_dp_idx  = pbuf_get_index('SNOW_DP')
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    if (clim_modal_aero) then

       ! If climate calculations are affected by prescribed modal aerosols, the
       ! the initialization routine for the dry mode radius calculation is called
       ! here.  For prognostic MAM the initialization is called from
       ! modal_aero_initialize
       if (.not. prog_modal_aero) then
          call modal_aero_calcsize_init(pbuf2d)
       endif

       call modal_aero_wateruptake_init(pbuf2d)

    end if

end subroutine phys_init

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! First part of atmospheric physics package before updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use time_manager,   only: get_nstep
    use cam_diagnostics,only: diag_allocate, diag_physvar_ic
    use check_energy,   only: check_energy_gmean

    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
#if (defined BFB_CAM_SCAM_IOP )
    use cam_history,    only: outfld
#endif
    use comsrf,         only: fsns, fsnt, flns, sgh30, flnt, landm, fsds
    use cam_abortutils, only: endrun
#if ( defined OFFLINE_DYN )
     use metdata,       only: get_met_srf1
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    !-----------------------------------------------------------------------
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! indices
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    call t_startf ('physpkg_st1')
    nstep = get_nstep()

#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SNOWH and TS for micro-phys
    !
    call get_met_srf1( cam_in )
#endif

    ! The following initialization depends on the import state (cam_in)
    ! being initialized.  This isn't true when cam_init is called, so need
    ! to postpone this initialization to here.
    if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

    ! Compute total energy of input state and previous output state
    call t_startf ('chk_en_gmean')
    call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
    call t_stopf ('chk_en_gmean')

    call t_stopf ('physpkg_st1')

    if ( adiabatic .or. ideal_phys )then
       call t_startf ('bc_physics')
       call phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
       call t_stopf ('bc_physics')
    else
       call t_startf ('physpkg_st1')

       call pbuf_allocate(pbuf2d, 'physpkg')
       call diag_allocate()

       !-----------------------------------------------------------------------
       ! Advance time information
       !-----------------------------------------------------------------------

       call phys_timestep_init( phys_state, cam_out, pbuf2d)

       call t_stopf ('physpkg_st1')

#ifdef TRACER_CHECK
       call gmean_mass ('before tphysbc DRY', phys_state)
#endif


       !-----------------------------------------------------------------------
       ! Tendency physics before flux coupler invocation
       !-----------------------------------------------------------------------
       !

#if (defined BFB_CAM_SCAM_IOP )
       do c=begchunk, endchunk
          call outfld('Tg',cam_in(c)%ts,pcols   ,c     )
       end do
#endif

       call t_barrierf('sync_bc_physics', mpicom)
       call t_startf ('bc_physics')
       call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, phys_buffer_chunk)
       do c=begchunk, endchunk
          !
          ! Output physics terms to IC file
          !
          phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

          call t_startf ('diag_physvar_ic')
          call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
          call t_stopf ('diag_physvar_ic')

          call tphysbc (ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                       phys_tend(c), phys_buffer_chunk,  fsds(1,c), landm(1,c),          &
                       sgh30(1,c), cam_out(c), cam_in(c) )

       end do

       call t_adj_detailf(-1)
       call t_stopf ('bc_physics')

       ! Don't call the rest in CRM mode
       if(single_column.and.scm_crm_mode) return

#ifdef TRACER_CHECK
       call gmean_mass ('between DRY', phys_state)
#endif
    end if

end subroutine phys_run1

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Physics for adiabatic or idealized physics case.
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_old_tim_idx
    use time_manager,     only: get_nstep
    use cam_diagnostics,  only: diag_phys_writeout
    use check_energy,     only: check_energy_fix, check_energy_chng
    use dycore,           only: dycore_is

    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer             :: c               ! indices
    integer             :: nstep           ! current timestep number
    type(physics_ptend) :: ptend(begchunk:endchunk) ! indivdual parameterization tendencies
    real(r8)            :: flx_heat(pcols) ! effective sensible heat flux
    real(r8)            :: zero(pcols)     ! array of zeros

    ! physics buffer field for total energy
    real(r8), pointer, dimension(:) :: teout
    logical, SAVE :: first_exec_of_phys_run1_adiabatic_or_ideal  = .TRUE.
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    zero  = 0._r8

    ! Associate pointers with physics buffer fields
    if (first_exec_of_phys_run1_adiabatic_or_ideal) then
       first_exec_of_phys_run1_adiabatic_or_ideal  = .FALSE.
    endif

!$OMP PARALLEL DO PRIVATE (C, FLX_HEAT)
    do c=begchunk, endchunk

       ! Initialize the physics tendencies to zero.
       call physics_tend_init(phys_tend(c))

       ! Dump dynamics variables to history buffers
       call diag_phys_writeout(phys_state(c))

       if (dycore_is('LR') .or. dycore_is('SE') ) then
          call check_energy_fix(phys_state(c), ptend(c), nstep, flx_heat)
          call physics_update(phys_state(c), ptend(c), ztodt, phys_tend(c))
          call check_energy_chng(phys_state(c), phys_tend(c), "chkengyfix", nstep, ztodt, &
               zero, zero, zero, flx_heat)
          call physics_ptend_dealloc(ptend(c))
       end if

       if ( ideal_phys )then
          call t_startf('tphysidl')
          call tphysidl(ztodt, phys_state(c), phys_tend(c))
          call t_stopf('tphysidl')
       end if

       ! Save total enery after physics for energy conservation checks
       call pbuf_set_field(pbuf_get_chunk(pbuf2d, c), teout_idx, phys_state(c)%te_cur)

    end do

end subroutine phys_run1_adiabatic_or_ideal

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d,  cam_out, &
       cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Second part of atmospheric physics package after updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
    use mo_lightning,   only: lightning_no_prod


    use cam_diagnostics,only: diag_deallocate, diag_surf
    use comsrf,         only: trefmxav, trefmnav, sgh, sgh30, fsds 
    use physconst,      only: stebol, latvap
    use carma_intr,     only: carma_accumulate_stats
#if ( defined OFFLINE_DYN )
    use metdata,        only: get_met_srf2
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),pointer, dimension(:,:)     :: pbuf2d

    type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
    type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
    !
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! chunk index
    integer :: ncol                              ! number of columns
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk
    !
    ! If exit condition just return
    !

    if(single_column.and.scm_crm_mode) return

    if ( adiabatic .or. ideal_phys ) return
    !-----------------------------------------------------------------------
    ! Tendency physics after coupler 
    ! Not necessary at terminal timestep.
    !-----------------------------------------------------------------------
    !
#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SHFLX QFLX TAUX TAUY for vert diffusion
    !
    call get_met_srf2( cam_in )
#endif
    ! Set lightning production of NO
    call t_startf ('lightning_no_prod')
    call lightning_no_prod( phys_state, pbuf2d,  cam_in )
    call t_stopf ('lightning_no_prod')

    call t_barrierf('sync_ac_physics', mpicom)
    call t_startf ('ac_physics')
    call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, NCOL, phys_buffer_chunk)

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       !
       ! surface diagnostics for history files
       !
       call t_startf('diag_surf')
       call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
       call t_stopf('diag_surf')

       call tphysac(ztodt, cam_in(c),  &
            sgh(1,c), sgh30(1,c), cam_out(c),                              &
            phys_state(c), phys_tend(c), phys_buffer_chunk,&
            fsds(1,c))
    end do                    ! Chunk loop

    call t_adj_detailf(-1)
    call t_stopf('ac_physics')

#ifdef TRACER_CHECK
    call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

    call t_startf ('carma_accumulate_stats')
    call carma_accumulate_stats()
    call t_stopf ('carma_accumulate_stats')

    call t_startf ('physpkg_st2')
    call pbuf_deallocate(pbuf2d, 'physpkg')

    call pbuf_update_tim_idx()
    call diag_deallocate()
    call t_stopf ('physpkg_st2')

end subroutine phys_run2

  !
  !----------------------------------------------------------------------- 
  !

subroutine phys_final( phys_state, phys_tend, pbuf2d )
    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
    use chemistry, only : chem_final
    use carma_intr, only : carma_final
    use wv_saturation, only : wv_sat_final
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Finalization of physics package
    ! 
    !-----------------------------------------------------------------------
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend ), pointer :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if(associated(pbuf2d)) then
       call pbuf_deallocate(pbuf2d,'global')
       deallocate(pbuf2d)
    end if
    deallocate(phys_state)
    deallocate(phys_tend)
    call chem_final
    call carma_final
    call wv_sat_final

end subroutine phys_final


subroutine tphysac (ztodt,   cam_in,  &
       sgh,     sgh30,                                     &
       cam_out,  state,   tend,    pbuf,            &
       fsds    )
    !----------------------------------------------------------------------- 
    ! 
    ! Tendency physics after coupling to land, sea, and ice models.
    !
    ! Computes the following:
    !
    !   o Aerosol Emission at Surface
    !   o Source-Sink for Advected Tracers
    !   o Symmetric Turbulence Scheme - Vertical Diffusion
    !   o Rayleigh Friction
    !   o Dry Deposition of Aerosol
    !   o Enforce Charge Neutrality ( Only for WACCM )
    !   o Gravity Wave Drag
    !   o QBO Relaxation ( Only for WACCM )
    !   o Ion Drag ( Only for WACCM )
    !   o Scale Dry Mass Energy
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use chemistry,          only: chem_is_active, chem_timestep_tend, chem_emissions
    use cam_diagnostics,    only: diag_phys_tend_writeout
    use gw_drag,            only: gw_tend
    use vertical_diffusion, only: vertical_diffusion_tend
    use rayleigh_friction,  only: rayleigh_friction_tend
    use constituents,       only: cnst_get_ind
    use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
         physics_dme_adjust, set_dry_to_wet, physics_state_check
    use majorsp_diffusion,  only: mspd_intr  ! WACCM-X major diffusion
    use ionosphere,         only: ionos_tend ! WACCM-X ionosphere
    use tracers,            only: tracers_timestep_tend
    use aoa_tracers,        only: aoa_tracers_timestep_tend
    use physconst,          only: rhoh2o, latvap,latice
    use aero_model,         only: aero_model_drydep
    use carma_intr,         only: carma_emission_tend, carma_timestep_tend
    use carma_flags_mod,    only: carma_do_aerosol, carma_do_emission
    use check_energy,       only: check_energy_chng
    use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
    use time_manager,       only: get_nstep
    use cam_abortutils,     only: endrun
    use dycore,             only: dycore_is
    use cam_control_mod,    only: aqua_planet 
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_set
    use charge_neutrality,  only: charge_fix
    use qbo,                only: qbo_relax
    use iondrag,            only: iondrag_calc, do_waccm_ions
    use clubb_intr,         only: clubb_surface
    use perf_mod
    use flux_avg,           only: flux_avg_run
    use unicon_cam,         only: unicon_cam_org_diags

    !Water Tracers:
    use water_tracer_vars,  only: trace_water, wtrc_iatype, wtrc_bulk_indices
    use water_tracers,      only: wtrc_check_h2o  

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
    real(r8), intent(in) :: fsds(pcols)            ! down solar flux
    real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
    real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)


    type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

    !
    !---------------------------Local workspace-----------------------------
    !
    type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

    integer  :: nstep                              ! current timestep number
    real(r8) :: zero(pcols)                        ! array of zeros

    integer :: lchnk                                ! chunk identifier
    integer :: ncol                                 ! number of atmospheric columns
    integer i,k,m                 ! Longitude, level indices
    integer :: yr, mon, day, tod       ! components of a date
    integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.

    logical :: labort                            ! abort flag

    real(r8) tvm(pcols,pver)           ! virtual temperature
    real(r8) prect(pcols)              ! total precipitation
    real(r8) surfric(pcols)            ! surface friction velocity
    real(r8) obklen(pcols)             ! Obukhov length
    real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry
    real(r8) :: flx_heat(pcols)        ! Heat flux for check_energy_chng.
    real(r8) :: tmp_q     (pcols,pver) ! tmp space
    real(r8) :: tmp_cldliq(pcols,pver) ! tmp space
    real(r8) :: tmp_cldice(pcols,pver) ! tmp space
    real(r8) :: tmp_t     (pcols,pver) ! tmp space

    !Water Tracers
    logical :: isOK                    ! Used to check that water tracer mass is being conserved.

    ! physics buffer fields for total energy and mass adjustment
    integer itim_old, ifld

    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction 

    !-----------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()
    
    ! Adjust the surface fluxes to reduce instabilities in near sfc layer
    if (phys_do_flux_avg()) then 
       call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)
    endif

    ! Validate the physics state.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysac")

    call t_startf('tphysac_init')
    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()


    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

    ifld = pbuf_get_index('AST')
    call pbuf_get_field(pbuf, ifld, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    !
    ! accumulate fluxes into net flux array for spectral dycores
    ! jrm Include latent heat of fusion for snow
    !
    call tphysac_flx_net_update(ncol, pcols, tend%flx_net, cam_in%shf, cam_out%precc, &
         cam_out%precl, cam_out%precsc, cam_out%precsl, latvap, latice, rhoh2o)

    ! emissions of aerosols and gas-phase chemistry constituents at surface
    call chem_emissions( state, cam_in )

    if (carma_do_emission) then
       ! carma emissions
       call carma_emission_tend (state, ptend, cam_in, ztodt)
       call physics_update(state, ptend, ztodt, tend)
    end if

    ! get nstep and zero array for energy checker
    zero = 0._r8
    nstep = get_nstep()
    call check_tracers_init(state, tracerint)

    ! Check if latent heat flux exceeds the total moisture content of the
    ! lowest model layer, thereby creating negative moisture.

    call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,               &
         state%q(1,pver,1),state%rpdel(1,pver) ,cam_in%shf ,         &
         cam_in%lhf , cam_in%cflx )

    call t_stopf('tphysac_init')
    !===================================================
    ! Source/sink terms for advected tracers.
    !===================================================
    call t_startf('adv_tracer_src_snk')
    ! Test tracers

    call tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    !===================================================
    ! Chemistry and MAM calculation
    ! MAM core aerosol conversion process is performed in the below 'chem_timestep_tend'.
    ! In addition, surface flux of aerosol species other than 'dust' and 'sea salt', and
    ! elevated emission of aerosol species are treated in 'chem_timestep_tend' before
    ! Gas chemistry and MAM core aerosol conversion. 
    ! Note that surface flux is not added into the atmosphere, but elevated emission is
    ! added into the atmosphere as tendency.
    !===================================================
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, &
            pbuf,  fh2o, fsds)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, &
            cam_in%cflx)
    end if
    call t_stopf('adv_tracer_src_snk')

     if(trace_water) then
       isOK = wtrc_check_h2o("after-tracer source/sink", state, state%q, ztodt)
     end if

    !===================================================
    ! Vertical diffusion/pbl calculation
    ! Call vertical diffusion code (pbl, free atmosphere and molecular)
    !===================================================

    ! If CLUBB is called, do not call vertical diffusion, but obukov length and
    !   surface friction velocity still need to be computed.  In addition, 
    !   surface fluxes need to be updated here for constituents 
    if (do_clubb_sgs) then

       call clubb_surface ( state, ptend, ztodt, cam_in, surfric, obklen)
       
       ! Update surface flux constituents 
       call physics_update(state, ptend, ztodt, tend)

    else

       call t_startf('vertical_diffusion_tend')
       call vertical_diffusion_tend (ztodt ,state ,cam_in%wsx, cam_in%wsy,   &
            cam_in%shf     ,cam_in%cflx     ,surfric  ,obklen   ,ptend    ,ast    ,&
            cam_in%ocnfrac  , cam_in%landfrac ,        &
            sgh30    ,pbuf )

    !------------------------------------------
    ! Call major diffusion for extended model
    !------------------------------------------
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
         call mspd_intr (ztodt    ,state    ,ptend)
       endif

       call physics_update(state, ptend, ztodt, tend)
       call t_stopf ('vertical_diffusion_tend')
 
    endif

   !-----------------------
   !check water tracer mass
   !-----------------------
   if(trace_water) then
     isOk = wtrc_check_h2o("after-pbl", state, state%q, ztodt) !<-Will always lose mass until CLM is up and running.
   end if

    !===================================================
    ! Rayleigh friction calculation
    !===================================================
    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    if (do_clubb_sgs) then
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    else
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, &
           zero, cam_in%shf)
    endif
    
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

    !  aerosol dry deposition processes
    call t_startf('aero_drydep')
    call aero_model_drydep( state, pbuf, obklen, surfric, cam_in, ztodt, cam_out, ptend )
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

   ! CARMA microphysics
   !
   ! NOTE: This does both the timestep_tend for CARMA aerosols as well as doing the dry
   ! deposition for CARMA aerosols. It needs to follow vertical_diffusion_tend, so that
   ! obklen and surfric have been calculated. It needs to follow aero_model_drydep, so
   ! that cam_out%xxxdryxxx fields have already been set for CAM aerosols and cam_out
   ! can be added to for CARMA aerosols.
   if (carma_do_aerosol) then
     call t_startf('carma_timestep_tend')
     call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, obklen=obklen, ustar=surfric)
     call physics_update(state, ptend, ztodt, tend)
   
     call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, zero, zero, zero)
     call t_stopf('carma_timestep_tend')
   end if


    !---------------------------------------------------------------------------------
    !   ... enforce charge neutrality
    !---------------------------------------------------------------------------------
    call charge_fix(state, pbuf)
     
    !===================================================
    ! Gravity wave drag
    !===================================================
    call t_startf('gw_tend')

    call gw_tend(state, sgh, pbuf, ztodt, ptend, cam_in, flx_heat)

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, &
         zero, zero, flx_heat)
    call t_stopf('gw_tend')

    ! QBO relaxation
    call qbo_relax(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "qborelax", nstep, ztodt, zero, zero, zero, zero)

    ! Ion drag calculation
    call t_startf ( 'iondrag' )

    if ( do_waccm_ions ) then
       call iondrag_calc( lchnk, ncol, state, ptend, pbuf,  ztodt )
    else
       call iondrag_calc( lchnk, ncol, state, ptend)
    endif
    !----------------------------------------------------------------------------
    ! Call ionosphere routines for extended model if mode is set to ionosphere
    !----------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call ionos_tend(state, ptend, pbuf, ztodt)
    endif

    call physics_update(state, ptend, ztodt, tend)

    !---------------------------------------------------------------------------------
    ! Enforce charge neutrality after O+ change from ionos_tend
    !---------------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call charge_fix(state, pbuf)
    endif

    ! Check energy integrals
    call check_energy_chng(state, tend, "iondrag", nstep, ztodt, zero, zero, zero, zero)

    call t_stopf  ( 'iondrag' )

    !-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    ! Save total energy for global fixer in next timestep (FV and SE dycores)
    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim_old/),(/pcols,1/))       

    if (shallow_scheme .eq. 'UNICON') then

       ! ------------------------------------------------------------------------
       ! Insert the organization-related heterogeneities computed inside the
       ! UNICON into the tracer arrays here before performing advection.
       ! This is necessary to prevent any modifications of organization-related
       ! heterogeneities by non convection-advection process, such as
       ! dry and wet deposition of aerosols, MAM, etc.
       ! Again, note that only UNICON and advection schemes are allowed to
       ! changes to organization at this stage, although we can include the
       ! effects of other physical processes in future.
       ! ------------------------------------------------------------------------

       call unicon_cam_org_diags(state, pbuf)

    end if

    !*** BAB's FV heating kludge *** apply the heating as temperature tendency.
    !*** BAB's FV heating kludge *** modify the temperature in the state structure
    call tphysac_t_update(ncol, pcols, pver, ztodt, state%t, tini, tend%dtdt, dtcore, tmp_t)

    !
    ! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
    !     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.
    if ( dycore_is('LR') .or. dycore_is('SE')) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist

    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call tphysac_q_snapshot(ncol, pcols, pver, pcnst, ixcldliq, ixcldice, &
         state%q, tmp_q, tmp_cldliq, tmp_cldice)
    call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)

    !-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if (aqua_planet) then
       labort = .false.
       do i=1,ncol
          if (cam_in%ocnfrac(i) /= 1._r8) labort = .true.
       end do
       if (labort) then
          call endrun ('TPHYSAC error:  grid contains non-ocean point')
       endif
    endif

    call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
         tmp_t, qini, cldliqini, cldiceini)

    call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

end subroutine tphysac

subroutine tphysbc (ztodt,               &
       fsns,    fsnt,    flns,    flnt,    state,   &
       tend,    pbuf,     fsds,    landm,            &
       sgh30, cam_out, cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Evaluate and apply physical processes that are calculated BEFORE 
    ! coupling to land, sea, and ice models.  
    !
    ! Processes currently included are: 
    !
    !  o Resetting Negative Tracers to Positive
    !  o Global Mean Total Energy Fixer
    !  o Dry Adjustment
    !  o Asymmetric Turbulence Scheme : Deep Convection & Shallow Convection
    !  o Stratiform Macro-Microphysics
    !  o Wet Scavenging of Aerosol
    !  o Radiation
    !
    ! Method: 
    !
    ! Each parameterization should be implemented with this sequence of calls:
    !  1)  Call physics interface
    !  2)  Check energy
    !  3)  Call physics_update
    ! See Interface to Column Physics and Chemistry Packages 
    !   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,          only : physics_buffer_desc, pbuf_get_field
    use physics_buffer,          only : pbuf_get_index, pbuf_old_tim_idx
    use physics_buffer,          only : col_type_subcol, dyn_time_lvls
    use shr_kind_mod,    only: r8 => shr_kind_r8

    use stratiform,      only: stratiform_tend
    use microp_driver,   only: microp_driver_tend
    use microp_aero,     only: microp_aero_run
    use macrop_driver,   only: macrop_driver_tend
    use physics_types,   only: physics_state, physics_tend, physics_ptend, &
         physics_update, physics_ptend_init, physics_ptend_sum, &
         physics_state_alloc, physics_state_dealloc, &
         physics_tend_alloc, physics_tend_init, physics_tend_dealloc, &
         physics_state_check, physics_ptend_scale
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld
    use physconst,       only: cpair, latvap
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
    use time_manager,    only: is_first_step, get_nstep
    use convect_shallow, only: convect_shallow_tend
    use check_energy,    only: check_energy_chng, check_energy_fix, check_energy_timestep_init
    use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
    use dycore,          only: dycore_is
    use aero_model,      only: aero_model_wetdep
    use carma_intr,      only: carma_wetdep_tend, carma_timestep_tend
    use carma_flags_mod, only: carma_do_detrain, carma_do_cldice, carma_do_cldliq,  carma_do_wetdep
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_adj
    use clubb_intr,      only: clubb_tend_cam
    use sslt_rebin,      only: sslt_rebin_adv
    use tropopause,      only: tropopause_output
    use cam_abortutils,  only: endrun
    use subcol,          only: subcol_gen, subcol_ptend_avg
    use subcol_utils,    only: subcol_ptend_copy, is_subcol_on

    !water tracers
    use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_bulk_indices,&
                                 wtrc_nwset, iwspec, wisotope
    use water_tracers,   only: wtrc_check_h2o, wtrc_ratio, wtrc_get_rstd, &
                               wtrc_is_tagged, wtrc_mass_fixer
    use water_types,     only: pwtype, iwtvap, iwtliq, iwtice


    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
    real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
    real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
    real(r8), intent(in) :: sgh30(pcols)                   ! Std. deviation of 30 s orography for tms

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_sc         ! state for sub-columns
    type(physics_ptend)   :: ptend_sc         ! ptend for sub-columns
    type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
    type(physics_ptend)   :: ptend_aero_sc    ! ptend for microp_aero on sub-columns
    type(physics_tend)    :: tend_sc          ! tend for sub-columns

    integer :: nstep                          ! current timestep number

    real(r8) :: net_flx(pcols)

    real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
    real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

    real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
    real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
    real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from shallow + deep convections
    real(r8) dlf2(pcols,pver)                  ! Detraining cld H20 from shallow convections
    real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
    real(r8) rtdt                              ! 1./ztodt

    integer lchnk                              ! chunk identifier
    integer ncol                               ! number of atmospheric columns
    integer ierr

    integer  i,k,m,n                           ! Longitude, level, constituent indices
    integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
    ! for macro/micro co-substepping
    integer :: macmic_it                       ! iteration variables
    real(r8) :: cld_macmic_ztodt               ! modified timestep
    ! physics buffer fields to compute tendencies for stratiform package
    integer itim_old, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore

    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

    ! convective precipitation variables
    real(r8),pointer :: prec_dp(:)                ! total precipitation from ZM convection
    real(r8),pointer :: snow_dp(:)                ! snow from ZM convection
    real(r8),pointer :: prec_sh(:)                ! total precipitation from Hack convection
    real(r8),pointer :: snow_sh(:)                ! snow from Hack convection

    ! carma precipitation variables
    real(r8) :: prec_sed_carma(pcols)          ! total precip from cloud sedimentation (CARMA)
    real(r8) :: snow_sed_carma(pcols)          ! snow from cloud ice sedimentation (CARMA)

    ! stratiform precipitation variables
    real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
    real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
    real(r8),pointer :: prec_str_sc(:)  ! sfc flux of precip from stratiform (m/s) -- for subcolumns
    real(r8),pointer :: snow_str_sc(:)  ! sfc flux of snow from stratiform   (m/s) -- for subcolumns
    real(r8),pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme
    real(r8),pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme
    real(r8),pointer :: prec_sed(:)     ! total precip from cloud sedimentation
    real(r8),pointer :: snow_sed(:)     ! snow from cloud ice sedimentation

    ! Local copies for substepping
    real(r8) :: prec_pcw_macmic(pcols)
    real(r8) :: snow_pcw_macmic(pcols)
    real(r8) :: prec_sed_macmic(pcols)
    real(r8) :: snow_sed_macmic(pcols)

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: zero_sc(pcols*psubcols)        ! array of zeros
    real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
    real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
    real(r8) :: det_s  (pcols)                 ! vertical integral of detrained static energy from ice
    real(r8) :: det_ice(pcols)                 ! vertical integral of detrained ice
    real(r8) :: flx_cnd(pcols)
    real(r8) :: flx_heat(pcols)
    type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
    real(r8) :: zero_tracers(pcols,pcnst)

    logical   :: lq(pcnst)
    ! For water tracers:
    logical  :: isOK
    integer  :: p    !for H2O mass fixing
    real(r8) :: R    !for H2O mass fixing
    real(r8) :: diff !for H2O mass fixing
    real(r8) :: oval !for H2O mass fixing 
    real(r8) :: wtdlf(pcols,pver,wtrc_nwset)

    !  pass macro to micro
    character(len=16) :: microp_scheme 
    character(len=16) :: macrop_scheme

    ! Debug physics_state.
    logical :: state_debug_checks
    !-----------------------------------------------------------------------

    call t_startf('bc_init')

    call phys_getopts( microp_scheme_out      = microp_scheme, &
                       macrop_scheme_out      = macrop_scheme, &
                       state_debug_checks_out = state_debug_checks)
    
    !-----------------------------------------------------------------------

    zero = 0._r8
    call tphysbc_zero_buffers(pcols, pcnst, pcols*psubcols, zero_tracers, zero_sc)

    lchnk = state%lchnk
    ncol  = state%ncol

    rtdt = 1._r8/ztodt

    nstep = get_nstep()

    !***********************************
    !Correct water tracer/isotope masses
    !***********************************
    if(trace_water) then
      call wtrc_mass_fixer(state)
    end if
   !*************************************************
   !Remove all water tracer values that are too large
   !*************************************************
   !NOTE:  In theory,  water tracers and tags can never
   !be larger than the bulk water. Any tracer quantity
   !that is larger than the bulk water is assumed to be
   !too large due to numerical errors, and thus it is ok
   !to go ahead and destroy that extra mass.  Water isotopes
   !can be larger than bulk water, at least in the condensed
   !phase, but only to a certain point (say, ~50 permil for
   !HDO).  Thus a cutoff of 1.1*bulk, instead of just bulk,
   !is used.  Still, this routine DESTROYS MASS, which is
   !generally not a good thing when modeling a physical system,
   !and it could potentially cover up actual physical errors
   !in the water tracer routines. Still, for now its better
   !to have the model stable via this fixer than for it to
   !be more "accurate" but blow up at high resolutions.
   !**************************************************
   if(trace_water) then
     call tphysbc_trace_water_clip(lchnk, ncol, pcols, pver, pcnst, state%q)
   end if
  !*******************************

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim_old/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )
    call tphysbc_init_fields(ncol, pcols, pver, pcnst, fracis, tend%dTdt, tend%dudt, tend%dvdt)

    !
    ! Make sure that input tracers are all positive (probably unnecessary)
    !
    call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )

    ! Verify state coming from the dynamics
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

    ! Since clybry_fam_adj operates directly on the tracers, and has no
    ! physics_update call, re-run qneg3.

    call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )

    ! Validate output of clybry_fam_adj.
    if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

    !
    ! Dump out "before physics" state
    !
    call diag_state_b4_phys_write (state)

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    call t_stopf('bc_init')

    !===================================================
    ! Global mean total energy fixer
    !===================================================
    call t_startf('energy_fixer')

    !*** BAB's FV heating kludge *** save the initial temperature
    call tphysbc_tini_copy(ncol, pcols, pver, state%t, tini)
    if (dycore_is('LR') .or. dycore_is('SE'))  then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    end if
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call tphysbc_qini_snapshot(ncol, pcols, pver, pcnst, ixcldliq, ixcldice, &
         state%q, qini, cldliqini, cldiceini)

    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! T tendency due to dynamics
    if( nstep > dyn_time_lvls-1 ) then
       call tphysbc_dtcore_update(ncol, pcols, pver, ztodt, tini, dtcore, tend%dTdt)
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')
    !
    !===================================================
    ! Dry adjustment
    ! This code block is not a good example of interfacing a parameterization
    !===================================================
    call t_startf('dry_adjustment')

    ! Copy state info for input to dadadj
    ! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

    call tphysbc_dadadj_lq_init(pcnst, lq)
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
    call tphysbc_dadadj_input(ncol, pcols, pver, pcnst, state%t, state%q, ptend%s, ptend%q)

    call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
         ptend%s, ptend%q(1,1,1))
    call tphysbc_dadadj_output(ncol, pcols, pver, pcnst, ztodt, cpair, state%t, state%q, ptend%s, ptend%q)
    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')
    !
    !===================================================
    ! Moist convection
    !===================================================
    call t_startf('moist_convection')
    !
    ! Since the PBL doesn't pass constituent perturbations, they
    ! are zeroed here for input to the moist convection routine
    !

    !-----------------------
    !Check water tracer mass
    !-----------------------

    if (trace_water) then
      isOK = wtrc_check_h2o("before-convection", state, state%q, ztodt) !<-will always indicates mass loss until Land is up and running
    end if

    call t_startf ('convect_deep_tend')
    call convect_deep_tend(  &
         cmfmc,      cmfcme,             &
         dlf,        pflx,    zdu,       &
         rliq,    &
         ztodt,   &
         state,   ptend, cam_in%landfrac, pbuf, wtdlf ) 
    call t_stopf('convect_deep_tend')

    call physics_update(state, ptend, ztodt, tend)

    !-----------------------
    !Check water tracer mass
    !-----------------------
    if (trace_water) then
      isOK = wtrc_check_h2o("after-deep tphysbc", state, state%q, ztodt)
    end if

    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
    call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
    call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
    call pbuf_get_field(pbuf, prec_str_idx, prec_str )
    call pbuf_get_field(pbuf, snow_str_idx, snow_str )
    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
    call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )
    call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw )
    call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw )

    if (use_subcol_microp) then
      call pbuf_get_field(pbuf, prec_str_idx, prec_str_sc, col_type=col_type_subcol)
      call pbuf_get_field(pbuf, snow_str_idx, snow_str_sc, col_type=col_type_subcol)
    end if

    ! Check energy integrals, including "reserved liquid"
    call tphysbc_flx_cnd_sum(ncol, pcols, prec_dp, rliq, flx_cnd)
    call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_dp, zero)

    !
    ! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
    !
    call t_startf ('convect_shallow_tend')

    call convect_shallow_tend (ztodt   , cmfmc,  cmfmc2  ,&
         dlf        , dlf2   ,  rliq   , rliq2, & 
         state      , ptend  ,  pbuf, sgh30, cam_in, wtdlf)
    call t_stopf ('convect_shallow_tend')

    call physics_update(state, ptend, ztodt, tend)

    !***********************************
    !Correct water tracer/isotope masses
    !***********************************
    if(trace_water) then
      call wtrc_mass_fixer(state)
    end if
    !***********************************

    !-----------------------
    !Check water tracer mass
    !-----------------------
    if (trace_water) then
      isOK = wtrc_check_h2o("after-shallow tphysbc", state, state%q, ztodt)
    end if

    call tphysbc_flx_cnd_sum(ncol, pcols, prec_sh, rliq2, flx_cnd)
    call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_sh, zero)

    call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

    call t_stopf('moist_convection')

    ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
    ! modes that correspond to the available optics data.  This is only necessary
    ! for CAM-RT.  But it's done here so that the microphysics code which is called
    ! from the stratiform interface has access to the same aerosols as the radiation
    ! code.
    call sslt_rebin_adv(pbuf,  state)
    
    !===================================================
    ! Calculate tendencies from CARMA bin microphysics.
    !===================================================
    !
    ! If CARMA is doing detrainment, then on output, rliq no longer represents water reserved
    ! for detrainment, but instead represents potential snow fall. The mass and number of the
    ! snow are stored in the physics buffer and will be incorporated by the MG microphysics.
    !
    ! Currently CARMA cloud microphysics is only supported with the MG microphysics.
    call t_startf('carma_timestep_tend')

    if (carma_do_cldice .or. carma_do_cldliq) then
       call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, dlf=dlf, rliq=rliq, &
            prec_str=prec_str, snow_str=snow_str, prec_sed=prec_sed_carma, snow_sed=snow_sed_carma)
       call physics_update(state, ptend, ztodt, tend)

       ! Before the detrainment, the reserved condensate is all liquid, but if CARMA is doing
       ! detrainment, then the reserved condensate is snow.
       if (carma_do_detrain) then
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str+rliq, snow_str+rliq, zero)
       else
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str, snow_str, zero)
       end if
    end if

    call t_stopf('carma_timestep_tend')

    if( microp_scheme == 'RK' ) then

       !===================================================
       ! Calculate stratiform tendency (sedimentation, detrain, cloud fraction and microphysics )
       !===================================================
       call t_startf('stratiform_tend')

       call stratiform_tend(state, ptend, pbuf, ztodt, &
            cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
            landm, cam_in%snowhland, & ! sediment
            dlf, dlf2, & ! detrain
            rliq  , & ! check energy after detrain
            cmfmc,   cmfmc2, &
            cam_in%ts,      cam_in%sst,        zdu)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call t_stopf('stratiform_tend')

    elseif( microp_scheme == 'MG' ) then
       ! Start co-substepping of macrophysics and microphysics
       cld_macmic_ztodt = ztodt/cld_macmic_num_steps

       ! Clear precip fields that should accumulate.
       call tphysbc_precip_ops(0, ncol, pcols, cld_macmic_num_steps, &
            prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
            prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
            prec_sed_carma, snow_sed_carma)

       do macmic_it = 1, cld_macmic_num_steps

          if (micro_do_icesupersat) then 

            !===================================================
            ! Aerosol Activation
            !===================================================
            call t_startf('microp_aero_run')
            call microp_aero_run(state, ptend, cld_macmic_ztodt, pbuf)
            call t_stopf('microp_aero_run')

            call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

            call physics_update(state, ptend, ztodt, tend)
            call check_energy_chng(state, tend, "mp_aero_tend", nstep, ztodt, zero, zero, zero, zero)      

          endif
          !===================================================
          ! Calculate macrophysical tendency (sedimentation, detrain, cloud fraction)
          !===================================================

          call t_startf('macrop_tend')

          ! don't call Park macrophysics if CLUBB is called
          if (macrop_scheme .ne. 'CLUBB_SGS') then

             call macrop_driver_tend( &
                  state,           ptend,          cld_macmic_ztodt, &
                  cam_in%landfrac, cam_in%ocnfrac, cam_in%snowhland, & ! sediment
                  dlf,             dlf2,           wtdlf,            & ! detrain
                  cmfmc,           cmfmc2,                           &
                  cam_in%ts,       cam_in%sst,     zdu,              &
                  pbuf,            det_s,          det_ice)

             !  Since we "added" the reserved liquid back in this routine, we need 
             !    to account for it in the energy checker
             call tphysbc_macrop_fluxes(1, ncol, pcols, rliq, det_s, flx_cnd, flx_heat)

             ! Unfortunately, physics_update does not know what time period
             ! "tend" is supposed to cover, and therefore can't update it
             ! with substeps correctly. For now, work around this by scaling
             ! ptend down by the number of substeps, then applying it for
             ! the full time (ztodt).
             call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)          
             call physics_update(state, ptend, ztodt, tend)
             call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, &
                  zero, flx_cnd/cld_macmic_num_steps, &
                  det_ice/cld_macmic_num_steps, flx_heat/cld_macmic_num_steps)
 
          else ! Calculate CLUBB macrophysics

             ! =====================================================
             !    CLUBB call (PBL, shallow convection, macrophysics)
             ! =====================================================  
   
             call clubb_tend_cam(state,ptend,pbuf,cld_macmic_ztodt,&
                cmfmc, cmfmc2, cam_in, sgh30, macmic_it, cld_macmic_num_steps, & 
                dlf, det_s, det_ice)

                !  Since we "added" the reserved liquid back in this routine, we need 
                !    to account for it in the energy checker
                call tphysbc_macrop_fluxes(2, ncol, pcols, rliq, det_s, flx_cnd, flx_heat, cam_in%shf)

                ! Unfortunately, physics_update does not know what time period
                ! "tend" is supposed to cover, and therefore can't update it
                ! with substeps correctly. For now, work around this by scaling
                ! ptend down by the number of substeps, then applying it for
                ! the full time (ztodt).
                call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)
                !    Update physics tendencies and copy state to state_eq, because that is 
                !      input for microphysics              
                call physics_update(state, ptend, ztodt, tend)
                call check_energy_chng(state, tend, "clubb_tend", nstep, ztodt, &
                     cam_in%lhf/latvap/cld_macmic_num_steps, flx_cnd/cld_macmic_num_steps, &
                     det_ice/cld_macmic_num_steps, flx_heat/cld_macmic_num_steps)
 
          endif

          call t_stopf('macrop_tend')

          !===================================================
          ! Calculate cloud microphysics 
          !===================================================

          if (is_subcol_on()) then
             ! Allocate sub-column structures. 
             call physics_state_alloc(state_sc, lchnk, psubcols*pcols)
             call physics_tend_alloc(tend_sc, psubcols*pcols)

             ! Generate sub-columns using the requested scheme
             call subcol_gen(state, tend, state_sc, tend_sc, pbuf)

             !Initialize check energy for subcolumns
             call check_energy_timestep_init(state_sc, tend_sc, pbuf, col_type_subcol)
          end if

          if (.not. micro_do_icesupersat) then 

            call t_startf('microp_aero_run')
            call microp_aero_run(state, ptend_aero, cld_macmic_ztodt, pbuf)
            call t_stopf('microp_aero_run')

          endif

          call t_startf('microp_tend')


          if (use_subcol_microp) then
             call microp_driver_tend(state_sc, ptend_sc, cld_macmic_ztodt, pbuf)

             ! Average the sub-column ptend for use in gridded update - will not contain ptend_aero
             call subcol_ptend_avg(ptend_sc, state_sc%ngrdcol, lchnk, ptend)

             ! Copy ptend_aero field to one dimensioned by sub-columns before summing with ptend
             call subcol_ptend_copy(ptend_aero, state_sc, ptend_aero_sc)
             call physics_ptend_sum(ptend_aero_sc, ptend_sc, state_sc%ncol)
             call physics_ptend_dealloc(ptend_aero_sc)

             ! Have to scale and apply for full timestep to get tend right
             ! (see above note for macrophysics).
             call physics_ptend_scale(ptend_sc, 1._r8/cld_macmic_num_steps, ncol)

             call physics_update (state_sc, ptend_sc, ztodt, tend_sc)
             call check_energy_chng(state_sc, tend_sc, "microp_tend_subcol", &
                  nstep, ztodt, zero_sc, prec_str_sc/cld_macmic_num_steps, &
                  snow_str_sc/cld_macmic_num_steps, zero_sc)

             call physics_state_dealloc(state_sc)
             call physics_tend_dealloc(tend_sc)
             call physics_ptend_dealloc(ptend_sc)
          else
             call microp_driver_tend(state, ptend, cld_macmic_ztodt, pbuf)
          end if
          ! combine aero and micro tendencies for the grid
          if (.not. micro_do_icesupersat) then
             call physics_ptend_sum(ptend_aero, ptend, ncol)
             call physics_ptend_dealloc(ptend_aero)
          endif

          ! Have to scale and apply for full timestep to get tend right
          ! (see above note for macrophysics).
          call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

          call physics_update (state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, &
               zero, prec_str/cld_macmic_num_steps, &
               snow_str/cld_macmic_num_steps, zero)

          call t_stopf('microp_tend')
          call tphysbc_precip_ops(1, ncol, pcols, cld_macmic_num_steps, &
               prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
               prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
               prec_sed_carma, snow_sed_carma)

       end do ! end substepping over macrophysics/microphysics

       call tphysbc_precip_ops(2, ncol, pcols, cld_macmic_num_steps, &
            prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
            prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
            prec_sed_carma, snow_sed_carma)

    endif

    ! Add the precipitation from CARMA to the precipitation from stratiform.
    if (carma_do_cldice .or. carma_do_cldliq) then
       call tphysbc_precip_ops(3, ncol, pcols, cld_macmic_num_steps, &
            prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
            prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
            prec_sed_carma, snow_sed_carma)
    end if

   if(trace_water) then
     call wtrc_mass_fixer(state)
   end if

    if ( .not. deep_scheme_does_scav_trans() ) then

       ! -------------------------------------------------------------------------------
       ! 1. Wet Scavenging of Aerosols by Convective and Stratiform Precipitation.
       ! 2. Convective Transport of Non-Water Aerosol Species.
       !
       !  . Aerosol wet chemistry determines scavenging fractions, and transformations
       !  . Then do convective transport of all trace species except qv,ql,qi.
       !  . We needed to do the scavenging first to determine the interstitial fraction.
       !  . When UNICON is used as unified convection, we should still perform
       !    wet scavenging but not 'convect_deep_tend2'.
       ! -------------------------------------------------------------------------------

       call t_startf('bc_aerosols')
       if (clim_modal_aero .and. .not. prog_modal_aero) then
          call modal_aero_calcsize_diag(state, pbuf)
          call modal_aero_wateruptake_dr(state, pbuf)
       endif
       call aero_model_wetdep( state, ztodt, dlf, cam_out, ptend, pbuf)
       call physics_update(state, ptend, ztodt, tend)


       if (carma_do_wetdep) then
          ! CARMA wet deposition
          !
          ! NOTE: It needs to follow aero_model_wetdep, so that cam_out%xxxwetxxx
          ! fields have already been set for CAM aerosols and cam_out can be added
          ! to for CARMA aerosols.
          call t_startf ('carma_wetdep_tend')
          call carma_wetdep_tend(state, ptend, ztodt, pbuf, dlf, cam_out)
          call physics_update(state, ptend, ztodt, tend)
          call t_stopf ('carma_wetdep_tend')
       end if

       call t_startf ('convect_deep_tend2')
       call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf ) 
       call physics_update(state, ptend, ztodt, tend)
       call t_stopf ('convect_deep_tend2')

       ! check tracer integrals
       call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt,  zero_tracers)

       call t_stopf('bc_aerosols')

   endif

   if(trace_water) then
     call wtrc_mass_fixer(state)
   end if

    !===================================================
    ! Moist physical parameteriztions complete: 
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, cam_out%psl)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')


    call radiation_tend(state,ptend, pbuf, &
         cam_out, cam_in, &
         cam_in%landfrac, cam_in%icefrac, cam_in%snowhland, &
         fsns,    fsnt, flns,    flnt,  &
         fsds, net_flx)

    ! Set net flux used by spectral dycores
    call tphysbc_radheat_flx_net(ncol, pcols, tend%flx_net, net_flx)
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

    call t_stopf('radiation')

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

end subroutine tphysbc

subroutine phys_timestep_init(phys_state, cam_out, pbuf2d)
!-----------------------------------------------------------------------------------
!
! Purpose: The place for parameterizations to call per timestep initializations.
!          Generally this is used to update time interpolated fields from boundary
!          datasets.
!
!-----------------------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc
  use carma_intr,          only: carma_timestep_init
  use ghg_data,            only: ghg_data_timestep_init
  use cam3_aero_data,      only: cam3_aero_data_on, cam3_aero_data_timestep_init
  use cam3_ozone_data,     only: cam3_ozone_data_on, cam3_ozone_data_timestep_init
  use radiation,           only: radiation_do
  use tracers,             only: tracers_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use qbo,                 only: qbo_timestep_init
  use efield,              only: get_efield
  use iondrag,             only: do_waccm_ions
  use perf_mod

  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  use prescribed_strataero, only: prescribed_strataero_adv


  implicit none

  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  logical :: do_cam3_aero_data
  logical :: do_cam3_ozone_data
  logical :: do_waccm_ions_local

  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)

  !-----------------------------------------------------------------------------

  call phys_timestep_init_select_impl()
  if (.not. use_native_phys_tstep_impl) then
     call phys_timestep_init_select_branches(cam3_aero_data_on, cam3_ozone_data_on, do_waccm_ions)
     do_cam3_aero_data = iand(phys_tstep_branch_mask, 1) /= 0
     do_cam3_ozone_data = iand(phys_tstep_branch_mask, 2) /= 0
     do_waccm_ions_local = iand(phys_tstep_branch_mask, 4) /= 0
  else
     do_cam3_aero_data = cam3_aero_data_on
     do_cam3_ozone_data = cam3_ozone_data_on
     do_waccm_ions_local = do_waccm_ions
  end if

  ! Chemistry surface values
  call chem_surfvals_set()

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state, pbuf2d)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)
  call prescribed_volcaero_adv(phys_state, pbuf2d)
  call prescribed_strataero_adv(phys_state, pbuf2d)

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, pbuf2d, cam_out)

  ! CAM3 prescribed aerosol masses
  if (do_cam3_aero_data) call cam3_aero_data_timestep_init(pbuf2d,  phys_state)

  ! CAM3 prescribed ozone data
  if (do_cam3_ozone_data) call cam3_ozone_data_timestep_init(pbuf2d,  phys_state)

  ! Time interpolate data models of gasses in pbuf2d
  call ghg_data_timestep_init(pbuf2d,  phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state, pbuf2d)
 
  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(pbuf2d, phys_state)

  !----------------------------------------------------------------------
  ! update QBO data for this time step
  !----------------------------------------------------------------------
  call qbo_timestep_init

  if (do_waccm_ions_local) then
     ! Compute the electric field
     call t_startf ('efield')
     call get_efield
     call t_stopf ('efield')
  endif

  call carma_timestep_init()

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

end subroutine phys_timestep_init

!=======================================================================

subroutine phys_timestep_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_tstep_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_TIMESTEP_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_tstep_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_tstep_impl = .false.
  end if

  phys_tstep_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_tstep_impl) then
        write(iulog,*) 'phys_timestep_init implementation = native'
     else
        write(iulog,*) 'phys_timestep_init implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_timestep_init_select_impl

!=======================================================================

subroutine phys_timestep_init_select_branches(cam3_aero_data_on, cam3_ozone_data_on, do_waccm_ions)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: cam3_aero_data_on
  logical, intent(in) :: cam3_ozone_data_on
  logical, intent(in) :: do_waccm_ions

  integer(c_int64_t), target :: branch_mask

  interface
     subroutine phys_timestep_init_select_branches_codon(cam3_aero_on_c, cam3_ozone_on_c, do_waccm_ions_c, branch_mask_p) &
          bind(c, name="phys_timestep_init_select_branches_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: cam3_aero_on_c, cam3_ozone_on_c, do_waccm_ions_c
       type(c_ptr), value :: branch_mask_p
     end subroutine phys_timestep_init_select_branches_codon
  end interface

  if (phys_tstep_branch_selected) return

  branch_mask = 0_c_int64_t
  call phys_timestep_init_select_branches_codon( &
       merge(1_c_int64_t, 0_c_int64_t, cam3_aero_data_on), &
       merge(1_c_int64_t, 0_c_int64_t, cam3_ozone_data_on), &
       merge(1_c_int64_t, 0_c_int64_t, do_waccm_ions), &
       c_loc(branch_mask) &
  )

  phys_tstep_branch_mask = int(branch_mask)
  phys_tstep_branch_selected = .true.

end subroutine phys_timestep_init_select_branches

!=======================================================================

subroutine tphysbc_precip_ops_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_precip_ops_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_PRECIP_OPS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_precip_ops_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_precip_ops_impl = .false.
  end if

  tphysbc_precip_ops_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_precip_ops_impl) then
        write(iulog,*) 'tphysbc_precip_ops implementation = native'
     else
        write(iulog,*) 'tphysbc_precip_ops implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_precip_ops_select_impl

!=======================================================================

subroutine tphysbc_precip_ops(mode, ncol, pcols_local, cld_macmic_num_steps_local, &
     prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
     prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
     prec_sed_carma, snow_sed_carma)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod,   only: r8 => shr_kind_r8

  integer, intent(in) :: mode
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: cld_macmic_num_steps_local
  real(r8), intent(inout), target :: prec_sed_macmic(pcols_local)
  real(r8), intent(inout), target :: snow_sed_macmic(pcols_local)
  real(r8), intent(inout), target :: prec_pcw_macmic(pcols_local)
  real(r8), intent(inout), target :: snow_pcw_macmic(pcols_local)
  real(r8), intent(inout), target :: prec_sed(pcols_local)
  real(r8), intent(inout), target :: snow_sed(pcols_local)
  real(r8), intent(inout), target :: prec_pcw(pcols_local)
  real(r8), intent(inout), target :: snow_pcw(pcols_local)
  real(r8), intent(inout), target :: prec_str(pcols_local)
  real(r8), intent(inout), target :: snow_str(pcols_local)
  real(r8), intent(in),    target :: prec_sed_carma(pcols_local)
  real(r8), intent(in),    target :: snow_sed_carma(pcols_local)

  integer(c_int64_t), target :: mode_c
  integer(c_int64_t), target :: ncol_c
  integer(c_int64_t), target :: pcols_c
  integer(c_int64_t), target :: cld_macmic_num_steps_c

  interface
     subroutine tphysbc_precip_ops_codon(mode_c, ncol_c, pcols_c, cld_macmic_num_steps_c, &
          prec_sed_macmic_p, snow_sed_macmic_p, prec_pcw_macmic_p, snow_pcw_macmic_p, &
          prec_sed_p, snow_sed_p, prec_pcw_p, snow_pcw_p, prec_str_p, snow_str_p, &
          prec_sed_carma_p, snow_sed_carma_p) bind(c, name="tphysbc_precip_ops_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: mode_c
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: cld_macmic_num_steps_c
       type(c_ptr), value :: prec_sed_macmic_p
       type(c_ptr), value :: snow_sed_macmic_p
       type(c_ptr), value :: prec_pcw_macmic_p
       type(c_ptr), value :: snow_pcw_macmic_p
       type(c_ptr), value :: prec_sed_p
       type(c_ptr), value :: snow_sed_p
       type(c_ptr), value :: prec_pcw_p
       type(c_ptr), value :: snow_pcw_p
       type(c_ptr), value :: prec_str_p
       type(c_ptr), value :: snow_str_p
       type(c_ptr), value :: prec_sed_carma_p
       type(c_ptr), value :: snow_sed_carma_p
     end subroutine tphysbc_precip_ops_codon
  end interface

  call tphysbc_precip_ops_select_impl()

  if (use_native_tphysbc_precip_ops_impl) then
     call tphysbc_precip_ops_native(mode, ncol, pcols_local, cld_macmic_num_steps_local, &
          prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
          prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
          prec_sed_carma, snow_sed_carma)
     return
  end if

  mode_c = int(mode, c_int64_t)
  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  cld_macmic_num_steps_c = int(cld_macmic_num_steps_local, c_int64_t)

  call tphysbc_precip_ops_codon(mode_c, ncol_c, pcols_c, cld_macmic_num_steps_c, &
       c_loc(prec_sed_macmic), c_loc(snow_sed_macmic), c_loc(prec_pcw_macmic), c_loc(snow_pcw_macmic), &
       c_loc(prec_sed), c_loc(snow_sed), c_loc(prec_pcw), c_loc(snow_pcw), c_loc(prec_str), c_loc(snow_str), &
       c_loc(prec_sed_carma), c_loc(snow_sed_carma))

end subroutine tphysbc_precip_ops

!=======================================================================

subroutine tphysbc_precip_ops_native(mode, ncol, pcols_local, cld_macmic_num_steps_local, &
     prec_sed_macmic, snow_sed_macmic, prec_pcw_macmic, snow_pcw_macmic, &
     prec_sed, snow_sed, prec_pcw, snow_pcw, prec_str, snow_str, &
     prec_sed_carma, snow_sed_carma)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: mode
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: cld_macmic_num_steps_local
  real(r8), intent(inout) :: prec_sed_macmic(pcols_local)
  real(r8), intent(inout) :: snow_sed_macmic(pcols_local)
  real(r8), intent(inout) :: prec_pcw_macmic(pcols_local)
  real(r8), intent(inout) :: snow_pcw_macmic(pcols_local)
  real(r8), intent(inout) :: prec_sed(pcols_local)
  real(r8), intent(inout) :: snow_sed(pcols_local)
  real(r8), intent(inout) :: prec_pcw(pcols_local)
  real(r8), intent(inout) :: snow_pcw(pcols_local)
  real(r8), intent(inout) :: prec_str(pcols_local)
  real(r8), intent(inout) :: snow_str(pcols_local)
  real(r8), intent(in)    :: prec_sed_carma(pcols_local)
  real(r8), intent(in)    :: snow_sed_carma(pcols_local)

  select case (mode)
  case (0)
     prec_sed_macmic = 0._r8
     snow_sed_macmic = 0._r8
     prec_pcw_macmic = 0._r8
     snow_pcw_macmic = 0._r8
  case (1)
     prec_sed_macmic(:ncol) = prec_sed_macmic(:ncol) + prec_sed(:ncol)
     snow_sed_macmic(:ncol) = snow_sed_macmic(:ncol) + snow_sed(:ncol)
     prec_pcw_macmic(:ncol) = prec_pcw_macmic(:ncol) + prec_pcw(:ncol)
     snow_pcw_macmic(:ncol) = snow_pcw_macmic(:ncol) + snow_pcw(:ncol)
  case (2)
     prec_sed(:ncol) = prec_sed_macmic(:ncol)/cld_macmic_num_steps_local
     snow_sed(:ncol) = snow_sed_macmic(:ncol)/cld_macmic_num_steps_local
     prec_pcw(:ncol) = prec_pcw_macmic(:ncol)/cld_macmic_num_steps_local
     snow_pcw(:ncol) = snow_pcw_macmic(:ncol)/cld_macmic_num_steps_local
     prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol)
     snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol)
  case (3)
     prec_sed(:ncol) = prec_sed(:ncol) + prec_sed_carma(:ncol)
     snow_sed(:ncol) = snow_sed(:ncol) + snow_sed_carma(:ncol)
  end select

end subroutine tphysbc_precip_ops_native

!=======================================================================

subroutine tphysac_flx_net_update_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysac_flx_net_update_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSAC_FLX_NET_UPDATE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysac_flx_net_update_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysac_flx_net_update_impl = .false.
  end if

  tphysac_flx_net_update_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysac_flx_net_update_impl) then
        write(iulog,*) 'tphysac_flx_net_update implementation = native'
     else
        write(iulog,*) 'tphysac_flx_net_update implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysac_flx_net_update_select_impl

!=======================================================================

subroutine tphysac_flx_net_update(ncol, pcols_local, tend_flx_net, cam_in_shf, cam_out_precc, &
     cam_out_precl, cam_out_precsc, cam_out_precsl, latvap_local, latice_local, rhoh2o_local)

  use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(inout), target :: tend_flx_net(pcols_local)
  real(r8), intent(in), target :: cam_in_shf(pcols_local)
  real(r8), intent(in), target :: cam_out_precc(pcols_local)
  real(r8), intent(in), target :: cam_out_precl(pcols_local)
  real(r8), intent(in), target :: cam_out_precsc(pcols_local)
  real(r8), intent(in), target :: cam_out_precsl(pcols_local)
  real(r8), intent(in) :: latvap_local
  real(r8), intent(in) :: latice_local
  real(r8), intent(in) :: rhoh2o_local

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c

  interface
     subroutine tphysac_flx_net_update_codon(ncol_c, pcols_c, tend_flx_net_p, cam_in_shf_p, cam_out_precc_p, &
          cam_out_precl_p, cam_out_precsc_p, cam_out_precsl_p, latvap_local, latice_local, rhoh2o_local) &
          bind(c, name="tphysac_flx_net_update_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       type(c_ptr), value :: tend_flx_net_p
       type(c_ptr), value :: cam_in_shf_p
       type(c_ptr), value :: cam_out_precc_p
       type(c_ptr), value :: cam_out_precl_p
       type(c_ptr), value :: cam_out_precsc_p
       type(c_ptr), value :: cam_out_precsl_p
       real(c_double), value :: latvap_local
       real(c_double), value :: latice_local
       real(c_double), value :: rhoh2o_local
     end subroutine tphysac_flx_net_update_codon
  end interface

  call tphysac_flx_net_update_select_impl()

  if (use_native_tphysac_flx_net_update_impl) then
     call tphysac_flx_net_update_native(ncol, pcols_local, tend_flx_net, cam_in_shf, cam_out_precc, &
          cam_out_precl, cam_out_precsc, cam_out_precsl, latvap_local, latice_local, rhoh2o_local)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)

  call tphysac_flx_net_update_codon(ncol_c, pcols_c, c_loc(tend_flx_net), c_loc(cam_in_shf), c_loc(cam_out_precc), &
       c_loc(cam_out_precl), c_loc(cam_out_precsc), c_loc(cam_out_precsl), latvap_local, latice_local, rhoh2o_local)

end subroutine tphysac_flx_net_update

!=======================================================================

subroutine tphysac_flx_net_update_native(ncol, pcols_local, tend_flx_net, cam_in_shf, cam_out_precc, &
     cam_out_precl, cam_out_precsc, cam_out_precsl, latvap_local, latice_local, rhoh2o_local)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(inout) :: tend_flx_net(pcols_local)
  real(r8), intent(in) :: cam_in_shf(pcols_local)
  real(r8), intent(in) :: cam_out_precc(pcols_local)
  real(r8), intent(in) :: cam_out_precl(pcols_local)
  real(r8), intent(in) :: cam_out_precsc(pcols_local)
  real(r8), intent(in) :: cam_out_precsl(pcols_local)
  real(r8), intent(in) :: latvap_local
  real(r8), intent(in) :: latice_local
  real(r8), intent(in) :: rhoh2o_local

  integer :: i

  do i = 1, ncol
     tend_flx_net(i) = tend_flx_net(i) + cam_in_shf(i) + (cam_out_precc(i) + cam_out_precl(i))*latvap_local*rhoh2o_local + &
          (cam_out_precsc(i) + cam_out_precsl(i))*latice_local*rhoh2o_local
  end do

end subroutine tphysac_flx_net_update_native

!=======================================================================

subroutine tphysac_t_update_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysac_t_update_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSAC_T_UPDATE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysac_t_update_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysac_t_update_impl = .false.
  end if

  tphysac_t_update_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysac_t_update_impl) then
        write(iulog,*) 'tphysac_t_update implementation = native'
     else
        write(iulog,*) 'tphysac_t_update implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysac_t_update_select_impl

!=======================================================================

subroutine tphysac_t_update(ncol, pcols_local, pver_local, ztodt, state_t, tini, tend_dtdt, dtcore, tmp_t)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod,   only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(inout), target :: state_t(pcols_local, pver_local)
  real(r8), intent(in),    target :: tini(pcols_local, pver_local)
  real(r8), intent(in),    target :: tend_dtdt(pcols_local, pver_local)
  real(r8), intent(inout), target :: dtcore(pcols_local, pver_local)
  real(r8), intent(inout), target :: tmp_t(pcols_local, pver_local)

  integer(c_int64_t), target :: ncol_c
  integer(c_int64_t), target :: pcols_c
  integer(c_int64_t), target :: pver_c

  interface
     subroutine tphysac_t_update_codon(ncol_c, pcols_c, pver_c, ztodt, state_t_p, tini_p, tend_dtdt_p, dtcore_p, tmp_t_p) &
          bind(c, name="tphysac_t_update_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       real(c_double), value :: ztodt
       type(c_ptr), value :: state_t_p
       type(c_ptr), value :: tini_p
       type(c_ptr), value :: tend_dtdt_p
       type(c_ptr), value :: dtcore_p
       type(c_ptr), value :: tmp_t_p
     end subroutine tphysac_t_update_codon
  end interface

  call tphysac_t_update_select_impl()

  if (use_native_tphysac_t_update_impl) then
     call tphysac_t_update_native(ncol, pcols_local, pver_local, ztodt, state_t, tini, tend_dtdt, dtcore, tmp_t)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)

  call tphysac_t_update_codon(ncol_c, pcols_c, pver_c, ztodt, &
       c_loc(state_t), c_loc(tini), c_loc(tend_dtdt), c_loc(dtcore), c_loc(tmp_t))

end subroutine tphysac_t_update

!=======================================================================

subroutine tphysac_t_update_native(ncol, pcols_local, pver_local, ztodt, state_t, tini, tend_dtdt, dtcore, tmp_t)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(inout) :: state_t(pcols_local, pver_local)
  real(r8), intent(in)    :: tini(pcols_local, pver_local)
  real(r8), intent(in)    :: tend_dtdt(pcols_local, pver_local)
  real(r8), intent(inout) :: dtcore(pcols_local, pver_local)
  real(r8), intent(inout) :: tmp_t(pcols_local, pver_local)

  integer :: k

  tmp_t(:ncol,:pver_local) = state_t(:ncol,:pver_local)
  state_t(:ncol,:pver_local) = tini(:ncol,:pver_local) + ztodt*tend_dtdt(:ncol,:pver_local)

  do k = 1, pver_local
     dtcore(:ncol,k) = state_t(:ncol,k)
  end do

end subroutine tphysac_t_update_native

!=======================================================================

subroutine tphysac_q_snapshot_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysac_q_snapshot_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSAC_Q_SNAPSHOT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysac_q_snapshot_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysac_q_snapshot_impl = .false.
  end if

  tphysac_q_snapshot_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysac_q_snapshot_impl) then
        write(iulog,*) 'tphysac_q_snapshot implementation = native'
     else
        write(iulog,*) 'tphysac_q_snapshot implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysac_q_snapshot_select_impl

!=======================================================================

subroutine tphysac_q_snapshot(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
     state_q, tmp_q, tmp_cldliq, tmp_cldice)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: ixcldliq
  integer, intent(in) :: ixcldice
  real(r8), intent(in), target :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout), target :: tmp_q(pcols_local, pver_local)
  real(r8), intent(inout), target :: tmp_cldliq(pcols_local, pver_local)
  real(r8), intent(inout), target :: tmp_cldice(pcols_local, pver_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: pcnst_c
  integer(c_int64_t) :: ixcldliq_c
  integer(c_int64_t) :: ixcldice_c

  interface
     subroutine tphysac_q_snapshot_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, &
          state_q_p, tmp_q_p, tmp_cldliq_p, tmp_cldice_p) bind(c, name="tphysac_q_snapshot_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       integer(c_int64_t), value :: ixcldliq_c
       integer(c_int64_t), value :: ixcldice_c
       type(c_ptr), value :: state_q_p
       type(c_ptr), value :: tmp_q_p
       type(c_ptr), value :: tmp_cldliq_p
       type(c_ptr), value :: tmp_cldice_p
     end subroutine tphysac_q_snapshot_codon
  end interface

  call tphysac_q_snapshot_select_impl()

  if (use_native_tphysac_q_snapshot_impl) then
     call tphysac_q_snapshot_native(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
          state_q, tmp_q, tmp_cldliq, tmp_cldice)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)
  ixcldliq_c = int(ixcldliq, c_int64_t)
  ixcldice_c = int(ixcldice, c_int64_t)

  call tphysac_q_snapshot_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, &
       c_loc(state_q), c_loc(tmp_q), c_loc(tmp_cldliq), c_loc(tmp_cldice))

end subroutine tphysac_q_snapshot

!=======================================================================

subroutine tphysac_q_snapshot_native(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
     state_q, tmp_q, tmp_cldliq, tmp_cldice)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: ixcldliq
  integer, intent(in) :: ixcldice
  real(r8), intent(in) :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout) :: tmp_q(pcols_local, pver_local)
  real(r8), intent(inout) :: tmp_cldliq(pcols_local, pver_local)
  real(r8), intent(inout) :: tmp_cldice(pcols_local, pver_local)

  tmp_q     (:ncol,:pver_local) = state_q(:ncol,:pver_local,1)
  tmp_cldliq(:ncol,:pver_local) = state_q(:ncol,:pver_local,ixcldliq)
  tmp_cldice(:ncol,:pver_local) = state_q(:ncol,:pver_local,ixcldice)

end subroutine tphysac_q_snapshot_native

!=======================================================================

subroutine tphysbc_qini_snapshot_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_qini_snapshot_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_QINI_SNAPSHOT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_qini_snapshot_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_qini_snapshot_impl = .false.
  end if

  tphysbc_qini_snapshot_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_qini_snapshot_impl) then
        write(iulog,*) 'tphysbc_qini_snapshot implementation = native'
     else
        write(iulog,*) 'tphysbc_qini_snapshot implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_qini_snapshot_select_impl

!=======================================================================

subroutine tphysbc_qini_snapshot(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
     state_q, qini, cldliqini, cldiceini)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: ixcldliq
  integer, intent(in) :: ixcldice
  real(r8), intent(in), target :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout), target :: qini(pcols_local, pver_local)
  real(r8), intent(inout), target :: cldliqini(pcols_local, pver_local)
  real(r8), intent(inout), target :: cldiceini(pcols_local, pver_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: pcnst_c
  integer(c_int64_t) :: ixcldliq_c
  integer(c_int64_t) :: ixcldice_c

  interface
     subroutine tphysbc_qini_snapshot_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, &
          state_q_p, qini_p, cldliqini_p, cldiceini_p) bind(c, name="tphysbc_qini_snapshot_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       integer(c_int64_t), value :: ixcldliq_c
       integer(c_int64_t), value :: ixcldice_c
       type(c_ptr), value :: state_q_p
       type(c_ptr), value :: qini_p
       type(c_ptr), value :: cldliqini_p
       type(c_ptr), value :: cldiceini_p
     end subroutine tphysbc_qini_snapshot_codon
  end interface

  call tphysbc_qini_snapshot_select_impl()

  if (use_native_tphysbc_qini_snapshot_impl) then
     call tphysbc_qini_snapshot_native(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
          state_q, qini, cldliqini, cldiceini)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)
  ixcldliq_c = int(ixcldliq, c_int64_t)
  ixcldice_c = int(ixcldice, c_int64_t)

  call tphysbc_qini_snapshot_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, &
       c_loc(state_q), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini))

end subroutine tphysbc_qini_snapshot

!=======================================================================

subroutine tphysbc_qini_snapshot_native(ncol, pcols_local, pver_local, pcnst_local, ixcldliq, ixcldice, &
     state_q, qini, cldliqini, cldiceini)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: ixcldliq
  integer, intent(in) :: ixcldice
  real(r8), intent(in) :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout) :: qini(pcols_local, pver_local)
  real(r8), intent(inout) :: cldliqini(pcols_local, pver_local)
  real(r8), intent(inout) :: cldiceini(pcols_local, pver_local)

  qini     (:ncol,:pver_local) = state_q(:ncol,:pver_local,1)
  cldliqini(:ncol,:pver_local) = state_q(:ncol,:pver_local,ixcldliq)
  cldiceini(:ncol,:pver_local) = state_q(:ncol,:pver_local,ixcldice)

end subroutine tphysbc_qini_snapshot_native

!=======================================================================

subroutine tphysbc_dadadj_input_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_dadadj_input_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_DADADJ_INPUT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_dadadj_input_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_dadadj_input_impl = .false.
  end if

  tphysbc_dadadj_input_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_dadadj_input_impl) then
        write(iulog,*) 'tphysbc_dadadj_input implementation = native'
     else
        write(iulog,*) 'tphysbc_dadadj_input implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_dadadj_input_select_impl

!=======================================================================

subroutine tphysbc_dadadj_input(ncol, pcols_local, pver_local, pcnst_local, state_t, state_q, ptend_s, ptend_q)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(in), target :: state_t(pcols_local, pver_local)
  real(r8), intent(in), target :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout), target :: ptend_s(pcols_local, pver_local)
  real(r8), intent(inout), target :: ptend_q(pcols_local, pver_local, pcnst_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: pcnst_c

  interface
     subroutine tphysbc_dadadj_input_codon(ncol_c, pcols_c, pver_c, pcnst_c, state_t_p, state_q_p, ptend_s_p, ptend_q_p) &
          bind(c, name="tphysbc_dadadj_input_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       type(c_ptr), value :: state_t_p
       type(c_ptr), value :: state_q_p
       type(c_ptr), value :: ptend_s_p
       type(c_ptr), value :: ptend_q_p
     end subroutine tphysbc_dadadj_input_codon
  end interface

  call tphysbc_dadadj_input_select_impl()

  if (use_native_tphysbc_dadadj_input_impl) then
     call tphysbc_dadadj_input_native(ncol, pcols_local, pver_local, pcnst_local, state_t, state_q, ptend_s, ptend_q)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)

  call tphysbc_dadadj_input_codon(ncol_c, pcols_c, pver_c, pcnst_c, &
       c_loc(state_t), c_loc(state_q), c_loc(ptend_s), c_loc(ptend_q))

end subroutine tphysbc_dadadj_input

!=======================================================================

subroutine tphysbc_dadadj_input_native(ncol, pcols_local, pver_local, pcnst_local, state_t, state_q, ptend_s, ptend_q)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(in) :: state_t(pcols_local, pver_local)
  real(r8), intent(in) :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout) :: ptend_s(pcols_local, pver_local)
  real(r8), intent(inout) :: ptend_q(pcols_local, pver_local, pcnst_local)

  ptend_s(:ncol,:pver_local)   = state_t(:ncol,:pver_local)
  ptend_q(:ncol,:pver_local,1) = state_q(:ncol,:pver_local,1)

end subroutine tphysbc_dadadj_input_native

!=======================================================================

subroutine tphysbc_dadadj_output_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_dadadj_output_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_DADADJ_OUTPUT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_dadadj_output_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_dadadj_output_impl = .false.
  end if

  tphysbc_dadadj_output_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_dadadj_output_impl) then
        write(iulog,*) 'tphysbc_dadadj_output implementation = native'
     else
        write(iulog,*) 'tphysbc_dadadj_output implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_dadadj_output_select_impl

!=======================================================================

subroutine tphysbc_dadadj_output(ncol, pcols_local, pver_local, pcnst_local, ztodt, cpair_local, &
     state_t, state_q, ptend_s, ptend_q)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(in) :: cpair_local
  real(r8), intent(in), target :: state_t(pcols_local, pver_local)
  real(r8), intent(in), target :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout), target :: ptend_s(pcols_local, pver_local)
  real(r8), intent(inout), target :: ptend_q(pcols_local, pver_local, pcnst_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: pcnst_c

  interface
     subroutine tphysbc_dadadj_output_codon(ncol_c, pcols_c, pver_c, pcnst_c, ztodt, cpair_local, &
          state_t_p, state_q_p, ptend_s_p, ptend_q_p) bind(c, name="tphysbc_dadadj_output_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       real(c_double), value :: ztodt
       real(c_double), value :: cpair_local
       type(c_ptr), value :: state_t_p
       type(c_ptr), value :: state_q_p
       type(c_ptr), value :: ptend_s_p
       type(c_ptr), value :: ptend_q_p
     end subroutine tphysbc_dadadj_output_codon
  end interface

  call tphysbc_dadadj_output_select_impl()

  if (use_native_tphysbc_dadadj_output_impl) then
     call tphysbc_dadadj_output_native(ncol, pcols_local, pver_local, pcnst_local, ztodt, cpair_local, &
          state_t, state_q, ptend_s, ptend_q)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)

  call tphysbc_dadadj_output_codon(ncol_c, pcols_c, pver_c, pcnst_c, ztodt, cpair_local, &
       c_loc(state_t), c_loc(state_q), c_loc(ptend_s), c_loc(ptend_q))

end subroutine tphysbc_dadadj_output

!=======================================================================

subroutine tphysbc_dadadj_output_native(ncol, pcols_local, pver_local, pcnst_local, ztodt, cpair_local, &
     state_t, state_q, ptend_s, ptend_q)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(in) :: cpair_local
  real(r8), intent(in) :: state_t(pcols_local, pver_local)
  real(r8), intent(in) :: state_q(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout) :: ptend_s(pcols_local, pver_local)
  real(r8), intent(inout) :: ptend_q(pcols_local, pver_local, pcnst_local)

  ptend_s(:ncol,:)   = (ptend_s(:ncol,:)   - state_t(:ncol,:))/ztodt * cpair_local
  ptend_q(:ncol,:,1) = (ptend_q(:ncol,:,1) - state_q(:ncol,:,1))/ztodt

end subroutine tphysbc_dadadj_output_native

!=======================================================================

subroutine tphysbc_dtcore_update_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_dtcore_update_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_DTCORE_UPDATE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_dtcore_update_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_dtcore_update_impl = .false.
  end if

  tphysbc_dtcore_update_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_dtcore_update_impl) then
        write(iulog,*) 'tphysbc_dtcore_update implementation = native'
     else
        write(iulog,*) 'tphysbc_dtcore_update implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_dtcore_update_select_impl

!=======================================================================

subroutine tphysbc_dtcore_update(ncol, pcols_local, pver_local, ztodt, tini, dtcore, tend_dtdt)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(in), target :: tini(pcols_local, pver_local)
  real(r8), intent(inout), target :: dtcore(pcols_local, pver_local)
  real(r8), intent(in), target :: tend_dtdt(pcols_local, pver_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c

  interface
     subroutine tphysbc_dtcore_update_codon(ncol_c, pcols_c, pver_c, ztodt, tini_p, dtcore_p, tend_dtdt_p) &
          bind(c, name="tphysbc_dtcore_update_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       real(c_double), value :: ztodt
       type(c_ptr), value :: tini_p
       type(c_ptr), value :: dtcore_p
       type(c_ptr), value :: tend_dtdt_p
     end subroutine tphysbc_dtcore_update_codon
  end interface

  call tphysbc_dtcore_update_select_impl()

  if (use_native_tphysbc_dtcore_update_impl) then
     call tphysbc_dtcore_update_native(ncol, pcols_local, pver_local, ztodt, tini, dtcore, tend_dtdt)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)

  call tphysbc_dtcore_update_codon(ncol_c, pcols_c, pver_c, ztodt, c_loc(tini), c_loc(dtcore), c_loc(tend_dtdt))

end subroutine tphysbc_dtcore_update

!=======================================================================

subroutine tphysbc_dtcore_update_native(ncol, pcols_local, pver_local, ztodt, tini, dtcore, tend_dtdt)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in) :: ztodt
  real(r8), intent(in) :: tini(pcols_local, pver_local)
  real(r8), intent(inout) :: dtcore(pcols_local, pver_local)
  real(r8), intent(in) :: tend_dtdt(pcols_local, pver_local)

  integer :: k

  do k = 1,pver_local
     dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend_dtdt(:ncol,k)
  end do

end subroutine tphysbc_dtcore_update_native

!=======================================================================

subroutine tphysbc_tini_copy_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_tini_copy_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_TINI_COPY_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_tini_copy_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_tini_copy_impl = .false.
  end if

  tphysbc_tini_copy_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_tini_copy_impl) then
        write(iulog,*) 'tphysbc_tini_copy implementation = native'
     else
        write(iulog,*) 'tphysbc_tini_copy implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_tini_copy_select_impl

!=======================================================================

subroutine tphysbc_tini_copy(ncol, pcols_local, pver_local, state_t, tini)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in), target :: state_t(pcols_local, pver_local)
  real(r8), intent(inout), target :: tini(pcols_local, pver_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c

  interface
     subroutine tphysbc_tini_copy_codon(ncol_c, pcols_c, pver_c, state_t_p, tini_p) &
          bind(c, name="tphysbc_tini_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       type(c_ptr), value :: state_t_p
       type(c_ptr), value :: tini_p
     end subroutine tphysbc_tini_copy_codon
  end interface

  call tphysbc_tini_copy_select_impl()

  if (use_native_tphysbc_tini_copy_impl) then
     call tphysbc_tini_copy_native(ncol, pcols_local, pver_local, state_t, tini)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)

  call tphysbc_tini_copy_codon(ncol_c, pcols_c, pver_c, c_loc(state_t), c_loc(tini))

end subroutine tphysbc_tini_copy

!=======================================================================

subroutine tphysbc_tini_copy_native(ncol, pcols_local, pver_local, state_t, tini)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  real(r8), intent(in) :: state_t(pcols_local, pver_local)
  real(r8), intent(inout) :: tini(pcols_local, pver_local)

  tini(:ncol,:pver_local) = state_t(:ncol,:pver_local)

end subroutine tphysbc_tini_copy_native

!=======================================================================

subroutine tphysbc_flx_cnd_sum_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_flx_cnd_sum_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_FLX_CND_SUM_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_flx_cnd_sum_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_flx_cnd_sum_impl = .false.
  end if

  tphysbc_flx_cnd_sum_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_flx_cnd_sum_impl) then
        write(iulog,*) 'tphysbc_flx_cnd_sum implementation = native'
     else
        write(iulog,*) 'tphysbc_flx_cnd_sum implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_flx_cnd_sum_select_impl

!=======================================================================

subroutine tphysbc_flx_cnd_sum(ncol, pcols_local, a, b, out)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(in), target :: a(pcols_local)
  real(r8), intent(in), target :: b(pcols_local)
  real(r8), intent(inout), target :: out(pcols_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c

  interface
     subroutine tphysbc_flx_cnd_sum_codon(ncol_c, pcols_c, a_p, b_p, out_p) &
          bind(c, name="tphysbc_flx_cnd_sum_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       type(c_ptr), value :: a_p
       type(c_ptr), value :: b_p
       type(c_ptr), value :: out_p
     end subroutine tphysbc_flx_cnd_sum_codon
  end interface

  call tphysbc_flx_cnd_sum_select_impl()

  if (use_native_tphysbc_flx_cnd_sum_impl) then
     call tphysbc_flx_cnd_sum_native(ncol, pcols_local, a, b, out)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)

  call tphysbc_flx_cnd_sum_codon(ncol_c, pcols_c, c_loc(a), c_loc(b), c_loc(out))

end subroutine tphysbc_flx_cnd_sum

!=======================================================================

subroutine tphysbc_flx_cnd_sum_native(ncol, pcols_local, a, b, out)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(in) :: a(pcols_local)
  real(r8), intent(in) :: b(pcols_local)
  real(r8), intent(inout) :: out(pcols_local)

  out(:ncol) = a(:ncol) + b(:ncol)

end subroutine tphysbc_flx_cnd_sum_native

!=======================================================================

subroutine tphysbc_macrop_fluxes_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_macrop_fluxes_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_MACROP_FLUXES_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_macrop_fluxes_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_macrop_fluxes_impl = .false.
  end if

  tphysbc_macrop_fluxes_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_macrop_fluxes_impl) then
        write(iulog,*) 'tphysbc_macrop_fluxes implementation = native'
     else
        write(iulog,*) 'tphysbc_macrop_fluxes implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_macrop_fluxes_select_impl

!=======================================================================

subroutine tphysbc_macrop_fluxes(mode, ncol, pcols_local, rliq, det_s, flx_cnd, flx_heat, shf)

  use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: mode
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(in), target :: rliq(pcols_local)
  real(r8), intent(in), target :: det_s(pcols_local)
  real(r8), intent(inout), target :: flx_cnd(pcols_local)
  real(r8), intent(inout), target :: flx_heat(pcols_local)
  real(r8), intent(in), target, optional :: shf(pcols_local)

  integer(c_int64_t) :: mode_c
  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  type(c_ptr) :: shf_p

  interface
     subroutine tphysbc_macrop_fluxes_codon(mode_c, ncol_c, pcols_c, rliq_p, det_s_p, flx_cnd_p, flx_heat_p, shf_p) &
          bind(c, name="tphysbc_macrop_fluxes_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: mode_c
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       type(c_ptr), value :: rliq_p
       type(c_ptr), value :: det_s_p
       type(c_ptr), value :: flx_cnd_p
       type(c_ptr), value :: flx_heat_p
       type(c_ptr), value :: shf_p
     end subroutine tphysbc_macrop_fluxes_codon
  end interface

  call tphysbc_macrop_fluxes_select_impl()

  if (use_native_tphysbc_macrop_fluxes_impl) then
     call tphysbc_macrop_fluxes_native(mode, ncol, pcols_local, rliq, det_s, flx_cnd, flx_heat, shf)
     return
  end if

  mode_c = int(mode, c_int64_t)
  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)

  if (present(shf)) then
     shf_p = c_loc(shf)
  else
     shf_p = c_null_ptr
  end if

  call tphysbc_macrop_fluxes_codon(mode_c, ncol_c, pcols_c, c_loc(rliq), c_loc(det_s), c_loc(flx_cnd), c_loc(flx_heat), shf_p)

end subroutine tphysbc_macrop_fluxes

!=======================================================================

subroutine tphysbc_macrop_fluxes_native(mode, ncol, pcols_local, rliq, det_s, flx_cnd, flx_heat, shf)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: mode
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(in) :: rliq(pcols_local)
  real(r8), intent(in) :: det_s(pcols_local)
  real(r8), intent(inout) :: flx_cnd(pcols_local)
  real(r8), intent(inout) :: flx_heat(pcols_local)
  real(r8), intent(in), optional :: shf(pcols_local)

  flx_cnd(:ncol) = -1._r8*rliq(:ncol)

  if (mode == 1) then
     flx_heat(:ncol) = det_s(:ncol)
  else
     flx_heat(:ncol) = shf(:ncol) + det_s(:ncol)
  end if

end subroutine tphysbc_macrop_fluxes_native

!=======================================================================

subroutine tphysbc_radheat_flx_net_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_radheat_flx_net_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_RADHEAT_FLX_NET_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_radheat_flx_net_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_radheat_flx_net_impl = .false.
  end if

  tphysbc_radheat_flx_net_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_radheat_flx_net_impl) then
        write(iulog,*) 'tphysbc_radheat_flx_net implementation = native'
     else
        write(iulog,*) 'tphysbc_radheat_flx_net implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_radheat_flx_net_select_impl

!=======================================================================

subroutine tphysbc_radheat_flx_net(ncol, pcols_local, tend_flx_net, net_flx)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(inout), target :: tend_flx_net(pcols_local)
  real(r8), intent(in), target :: net_flx(pcols_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c

  interface
     subroutine tphysbc_radheat_flx_net_codon(ncol_c, pcols_c, tend_flx_net_p, net_flx_p) &
          bind(c, name="tphysbc_radheat_flx_net_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       type(c_ptr), value :: tend_flx_net_p
       type(c_ptr), value :: net_flx_p
     end subroutine tphysbc_radheat_flx_net_codon
  end interface

  call tphysbc_radheat_flx_net_select_impl()

  if (use_native_tphysbc_radheat_flx_net_impl) then
     call tphysbc_radheat_flx_net_native(ncol, pcols_local, tend_flx_net, net_flx)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)

  call tphysbc_radheat_flx_net_codon(ncol_c, pcols_c, c_loc(tend_flx_net), c_loc(net_flx))

end subroutine tphysbc_radheat_flx_net

!=======================================================================

subroutine tphysbc_radheat_flx_net_native(ncol, pcols_local, tend_flx_net, net_flx)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  real(r8), intent(inout) :: tend_flx_net(pcols_local)
  real(r8), intent(in) :: net_flx(pcols_local)

  tend_flx_net(:ncol) = net_flx(:ncol)

end subroutine tphysbc_radheat_flx_net_native

!=======================================================================

subroutine tphysbc_zero_buffers_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_zero_buffers_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_ZERO_BUFFERS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_zero_buffers_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_zero_buffers_impl = .false.
  end if

  tphysbc_zero_buffers_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_zero_buffers_impl) then
        write(iulog,*) 'tphysbc_zero_buffers implementation = native'
     else
        write(iulog,*) 'tphysbc_zero_buffers implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_zero_buffers_select_impl

!=======================================================================

subroutine tphysbc_zero_buffers(pcols_local, pcnst_local, zero_sc_len, zero_tracers, zero_sc)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: zero_sc_len
  real(r8), intent(inout), target :: zero_tracers(pcols_local, pcnst_local)
  real(r8), intent(inout), target :: zero_sc(zero_sc_len)

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pcnst_c
  integer(c_int64_t) :: zero_sc_len_c

  interface
     subroutine tphysbc_zero_buffers_codon(pcols_c, pcnst_c, zero_sc_len_c, zero_tracers_p, zero_sc_p) &
          bind(c, name="tphysbc_zero_buffers_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pcnst_c
       integer(c_int64_t), value :: zero_sc_len_c
       type(c_ptr), value :: zero_tracers_p
       type(c_ptr), value :: zero_sc_p
     end subroutine tphysbc_zero_buffers_codon
  end interface

  call tphysbc_zero_buffers_select_impl()

  if (use_native_tphysbc_zero_buffers_impl) then
     call tphysbc_zero_buffers_native(pcols_local, pcnst_local, zero_sc_len, zero_tracers, zero_sc)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)
  zero_sc_len_c = int(zero_sc_len, c_int64_t)

  call tphysbc_zero_buffers_codon(pcols_c, pcnst_c, zero_sc_len_c, c_loc(zero_tracers), c_loc(zero_sc))

end subroutine tphysbc_zero_buffers

!=======================================================================

subroutine tphysbc_zero_buffers_native(pcols_local, pcnst_local, zero_sc_len, zero_tracers, zero_sc)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: zero_sc_len
  real(r8), intent(inout) :: zero_tracers(pcols_local, pcnst_local)
  real(r8), intent(inout) :: zero_sc(zero_sc_len)

  zero_tracers(:,:) = 0._r8
  zero_sc(:) = 0._r8

end subroutine tphysbc_zero_buffers_native

!=======================================================================

subroutine tphysbc_trace_water_clip_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_trace_water_clip_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_TRACE_WATER_CLIP_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_trace_water_clip_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_trace_water_clip_impl = .false.
  end if

  tphysbc_trace_water_clip_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_trace_water_clip_impl) then
        write(iulog,*) 'tphysbc_trace_water_clip implementation = native'
     else
        write(iulog,*) 'tphysbc_trace_water_clip implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_trace_water_clip_select_impl

!=======================================================================

subroutine tphysbc_trace_water_clip(lchnk, ncol, pcols_local, pver_local, pcnst_local, state_q)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: qmin
  use water_tracer_vars, only: wtrc_iatype, wtrc_nwset, iwspec, wisotope
  use water_tracers, only: wtrc_get_rstd, wtrc_is_tagged
  use water_types, only: pwtype

  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(inout), target :: state_q(pcols_local, pver_local, pcnst_local)

  integer :: m, p, n
  integer(c_int64_t), target :: wtrc_iatype64(wtrc_nwset, pwtype)
  integer(c_int64_t), target :: tagged64(pcnst_local)
  real(r8), target :: rstd_by_constituent(pcnst_local)

  interface
     subroutine tphysbc_trace_water_clip_codon(ncol_c, pcols_c, pver_c, pcnst_c, pwtype_c, &
          wtrc_nwset_c, wisotope_on_c, state_q_p, wtrc_iatype_p, tagged_p, rstd_p) &
          bind(c, name="tphysbc_trace_water_clip_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       integer(c_int64_t), value :: pwtype_c
       integer(c_int64_t), value :: wtrc_nwset_c
       integer(c_int64_t), value :: wisotope_on_c
       type(c_ptr), value :: state_q_p
       type(c_ptr), value :: wtrc_iatype_p
       type(c_ptr), value :: tagged_p
       type(c_ptr), value :: rstd_p
     end subroutine tphysbc_trace_water_clip_codon
  end interface

  call tphysbc_trace_water_clip_select_impl()

  if (use_native_tphysbc_trace_water_clip_impl) then
     call tphysbc_trace_water_clip_native(lchnk, ncol, pcols_local, pver_local, pcnst_local, state_q)
     return
  end if

  do p = 1, pwtype
     do m = 1, wtrc_nwset
        wtrc_iatype64(m, p) = int(wtrc_iatype(m, p), c_int64_t)
     end do
  end do

  do m = 1, pcnst_local
     tagged64(m) = merge(1_c_int64_t, 0_c_int64_t, wtrc_is_tagged(m))
     rstd_by_constituent(m) = wtrc_get_rstd(iwspec(m))
  end do

  call tphysbc_trace_water_clip_codon( &
       int(ncol, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
       int(pcnst_local, c_int64_t), int(pwtype, c_int64_t), int(wtrc_nwset, c_int64_t), &
       merge(1_c_int64_t, 0_c_int64_t, wisotope), c_loc(state_q), c_loc(wtrc_iatype64), &
       c_loc(tagged64), c_loc(rstd_by_constituent) &
  )

  do p = 1, pwtype
     do m = 2, wtrc_nwset
        n = wtrc_iatype(m, p)
        call qneg3('wiso', lchnk, ncol, pcols_local, pver_local, n, n, qmin(n), state_q(1,1,n))
     end do
  end do

end subroutine tphysbc_trace_water_clip

!=======================================================================

subroutine tphysbc_trace_water_clip_native(lchnk, ncol, pcols_local, pver_local, pcnst_local, state_q)

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: qmin
  use water_tracer_vars, only: wtrc_iatype, wtrc_nwset, iwspec, wisotope
  use water_tracers, only: wtrc_get_rstd, wtrc_is_tagged
  use water_types, only: pwtype

  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(inout) :: state_q(pcols_local, pver_local, pcnst_local)

  integer :: i, k, m, n, p

  do i = 1, ncol
     do k = 1, pver_local
        do p = 1, pwtype
           do m = 2, wtrc_nwset
              if (wisotope) then
                 if (state_q(i,k,wtrc_iatype(m,p)) .gt. 1.5_r8*state_q(i,k,wtrc_iatype(1,p))) then
                    state_q(i,k,wtrc_iatype(m,p)) = state_q(i,k,wtrc_iatype(1,p))
                 end if
              else
                 if (state_q(i,k,wtrc_iatype(m,p)) .gt. state_q(i,k,wtrc_iatype(1,p))) then
                    if (wtrc_is_tagged(wtrc_iatype(m,p))) then
                       state_q(i,k,wtrc_iatype(m,p)) = state_q(i,k,wtrc_iatype(1,p))
                    else
                       state_q(i,k,wtrc_iatype(m,p)) = wtrc_get_rstd(iwspec(wtrc_iatype(m,p)))* &
                            state_q(i,k,wtrc_iatype(1,p))
                    end if
                 end if
              end if
           end do
        end do
     end do
  end do

  do p = 1, pwtype
     do m = 2, wtrc_nwset
        n = wtrc_iatype(m, p)
        call qneg3('wiso', lchnk, ncol, pcols_local, pver_local, n, n, qmin(n), state_q(1,1,n))
     end do
  end do

end subroutine tphysbc_trace_water_clip_native

!=======================================================================

subroutine tphysbc_dadadj_lq_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_dadadj_lq_init_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_DADADJ_LQ_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_dadadj_lq_init_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_dadadj_lq_init_impl = .false.
  end if

  tphysbc_dadadj_lq_init_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_dadadj_lq_init_impl) then
        write(iulog,*) 'tphysbc_dadadj_lq_init implementation = native'
     else
        write(iulog,*) 'tphysbc_dadadj_lq_init implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_dadadj_lq_init_select_impl

!=======================================================================

subroutine tphysbc_dadadj_lq_init(pcnst_local, lq)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  integer, intent(in) :: pcnst_local
  logical, intent(out) :: lq(pcnst_local)

  integer :: m
  integer(c_int64_t), target :: lq_mask(pcnst_local)

  interface
     subroutine tphysbc_dadadj_lq_init_codon(pcnst_c, lq_mask_p) bind(c, name="tphysbc_dadadj_lq_init_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcnst_c
       type(c_ptr), value :: lq_mask_p
     end subroutine tphysbc_dadadj_lq_init_codon
  end interface

  call tphysbc_dadadj_lq_init_select_impl()

  if (use_native_tphysbc_dadadj_lq_init_impl) then
     call tphysbc_dadadj_lq_init_native(pcnst_local, lq)
     return
  end if

  call tphysbc_dadadj_lq_init_codon(int(pcnst_local, c_int64_t), c_loc(lq_mask))

  do m = 1, pcnst_local
     lq(m) = lq_mask(m) /= 0_c_int64_t
  end do

end subroutine tphysbc_dadadj_lq_init

!=======================================================================

subroutine tphysbc_dadadj_lq_init_native(pcnst_local, lq)

  integer, intent(in) :: pcnst_local
  logical, intent(out) :: lq(pcnst_local)

  lq(:) = .FALSE.
  lq(1) = .TRUE.

end subroutine tphysbc_dadadj_lq_init_native

!=======================================================================

subroutine phys_inidat_batch_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('PHYS_INIDAT_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine phys_inidat_batch_append_proof

!=======================================================================

subroutine phys_inidat_batch_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_batch_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_BATCH_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_batch_impl = .false.
  end if

  phys_inidat_batch_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_batch_impl) then
        write(iulog,*) 'phys_inidat_batch implementation = native'
        call phys_inidat_batch_append_proof('phys_inidat_batch selector entered implementation = native')
     else
        write(iulog,*) 'phys_inidat_batch implementation = codon'
        call phys_inidat_batch_append_proof('phys_inidat_batch selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_batch_select_impl

!=======================================================================

subroutine phys_inidat_batch_log_entered()

  if (phys_inidat_batch_entered_logged) return
  phys_inidat_batch_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'phys_inidat_batch entered (qpert/pblh/tpert/cush/tke/kvm/kvh/water/cloud/tbot direct = codon)'
     call phys_inidat_batch_append_proof( &
          'phys_inidat_batch entered (qpert/pblh/tpert/cush/tke/kvm/kvh/water/cloud/tbot direct = codon)')
     call flush(iulog)
  end if

end subroutine phys_inidat_batch_log_entered

!=======================================================================

subroutine phys_inidat_qpert_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_qpert_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_QPERT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_qpert_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_qpert_default_impl = .false.
  end if

  phys_inidat_qpert_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_qpert_default_impl) then
        write(iulog,*) 'phys_inidat_qpert_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_qpert_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_qpert_default_select_impl

!=======================================================================

subroutine phys_inidat_qpert_default(pcols_local, chunk_count_local, found_local, tptr)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr(pcols_local, chunk_count_local)

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_qpert_default_codon(pcols_c, chunk_count_c, found_c, tptr_p) &
          bind(c, name="phys_inidat_batch_qpert_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr_p
     end subroutine phys_inidat_batch_qpert_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_qpert_default_native(pcols_local, chunk_count_local, found_local, tptr)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_qpert_default_codon(pcols_c, chunk_count_c, found_c, c_loc(tptr))

end subroutine phys_inidat_qpert_default

!=======================================================================

subroutine phys_inidat_qpert_default_native(pcols_local, chunk_count_local, found_local, tptr)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr(pcols_local, chunk_count_local)

  if (.not. found_local) then
     tptr = 0._r8
  end if

end subroutine phys_inidat_qpert_default_native

!=======================================================================

subroutine phys_inidat_qpert_expand_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_qpert_expand_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_QPERT_EXPAND_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_qpert_expand_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_qpert_expand_impl = .false.
  end if

  phys_inidat_qpert_expand_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_qpert_expand_impl) then
        write(iulog,*) 'phys_inidat_qpert_expand implementation = native'
     else
        write(iulog,*) 'phys_inidat_qpert_expand implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_qpert_expand_select_impl

!=======================================================================

subroutine phys_inidat_qpert_expand(pcols_local, pcnst_local, chunk_count_local, tptr, tptr3d_2)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: chunk_count_local
  real(r8), intent(in), target :: tptr(pcols_local, chunk_count_local)
  real(r8), intent(inout), target :: tptr3d_2(pcols_local, pcnst_local, chunk_count_local)

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pcnst_c
  integer(c_int64_t) :: chunk_count_c

  interface
     subroutine phys_inidat_batch_qpert_expand_codon(pcols_c, pcnst_c, chunk_count_c, tptr_p, tptr3d_2_p) &
          bind(c, name="phys_inidat_batch_qpert_expand_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pcnst_c
       integer(c_int64_t), value :: chunk_count_c
       type(c_ptr), value :: tptr_p
       type(c_ptr), value :: tptr3d_2_p
     end subroutine phys_inidat_batch_qpert_expand_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_qpert_expand_native(pcols_local, pcnst_local, chunk_count_local, tptr, tptr3d_2)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_qpert_expand_codon(pcols_c, pcnst_c, chunk_count_c, c_loc(tptr), c_loc(tptr3d_2))

end subroutine phys_inidat_qpert_expand

!=======================================================================

subroutine phys_inidat_qpert_expand_native(pcols_local, pcnst_local, chunk_count_local, tptr, tptr3d_2)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pcnst_local
  integer, intent(in) :: chunk_count_local
  real(r8), intent(in) :: tptr(pcols_local, chunk_count_local)
  real(r8), intent(inout) :: tptr3d_2(pcols_local, pcnst_local, chunk_count_local)

  tptr3d_2 = 0._r8
  tptr3d_2(:,1,:) = tptr(:,:)

end subroutine phys_inidat_qpert_expand_native

!=======================================================================

subroutine phys_inidat_pblh_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_pblh_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_PBLH_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_pblh_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_pblh_default_impl = .false.
  end if

  phys_inidat_pblh_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_pblh_default_impl) then
        write(iulog,*) 'phys_inidat_pblh_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_pblh_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_pblh_default_select_impl

!=======================================================================

subroutine phys_inidat_pblh_default(pcols_local, chunk_count_local, found_local, tptr)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr(pcols_local, chunk_count_local)

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_pblh_default_codon(pcols_c, chunk_count_c, found_c, tptr_p) &
          bind(c, name="phys_inidat_batch_pblh_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr_p
     end subroutine phys_inidat_batch_pblh_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_pblh_default_native(pcols_local, chunk_count_local, found_local, tptr)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_pblh_default_codon(pcols_c, chunk_count_c, found_c, c_loc(tptr))

end subroutine phys_inidat_pblh_default

!=======================================================================

subroutine phys_inidat_pblh_default_native(pcols_local, chunk_count_local, found_local, tptr)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr(pcols_local, chunk_count_local)

  if (.not. found_local) then
     tptr(:,:) = 0._r8
  end if

end subroutine phys_inidat_pblh_default_native

!=======================================================================

subroutine phys_inidat_tpert_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_tpert_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_TPERT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_tpert_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_tpert_default_impl = .false.
  end if

  phys_inidat_tpert_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_tpert_default_impl) then
        write(iulog,*) 'phys_inidat_tpert_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_tpert_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_tpert_default_select_impl

!=======================================================================

subroutine phys_inidat_tpert_default(pcols_local, chunk_count_local, found_local, tptr)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr(pcols_local, chunk_count_local)

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_tpert_default_codon(pcols_c, chunk_count_c, found_c, tptr_p) &
          bind(c, name="phys_inidat_batch_tpert_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr_p
     end subroutine phys_inidat_batch_tpert_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_tpert_default_native(pcols_local, chunk_count_local, found_local, tptr)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_tpert_default_codon(pcols_c, chunk_count_c, found_c, c_loc(tptr))

end subroutine phys_inidat_tpert_default

!=======================================================================

subroutine phys_inidat_tpert_default_native(pcols_local, chunk_count_local, found_local, tptr)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr(pcols_local, chunk_count_local)

  if (.not. found_local) then
     tptr(:,:) = 0._r8
  end if

end subroutine phys_inidat_tpert_default_native

!=======================================================================

subroutine phys_inidat_cush_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_cush_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_CUSH_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_cush_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_cush_default_impl = .false.
  end if

  phys_inidat_cush_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_cush_default_impl) then
        write(iulog,*) 'phys_inidat_cush_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_cush_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_cush_default_select_impl

!=======================================================================

subroutine phys_inidat_cush_default(pcols_local, chunk_count_local, found_local, tptr, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr(pcols_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_cush_default_codon(pcols_c, chunk_count_c, found_c, tptr_p, default_value) &
          bind(c, name="phys_inidat_batch_cush_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_cush_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_cush_default_native(pcols_local, chunk_count_local, found_local, tptr, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_cush_default_codon(pcols_c, chunk_count_c, found_c, c_loc(tptr), default_value)

end subroutine phys_inidat_cush_default

!=======================================================================

subroutine phys_inidat_cush_default_native(pcols_local, chunk_count_local, found_local, tptr, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr(pcols_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr = default_value
  end if

end subroutine phys_inidat_cush_default_native

!=======================================================================

subroutine phys_inidat_tke_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_tke_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_TKE_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_tke_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_tke_default_impl = .false.
  end if

  phys_inidat_tke_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_tke_default_impl) then
        write(iulog,*) 'phys_inidat_tke_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_tke_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_tke_default_select_impl

!=======================================================================

subroutine phys_inidat_tke_default(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pverp_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_tke_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, tptr3d_p, default_value) &
          bind(c, name="phys_inidat_batch_tke_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pverp_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr3d_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_tke_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_tke_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pverp_c = int(pverp_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_tke_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, c_loc(tptr3d), default_value)

end subroutine phys_inidat_tke_default

!=======================================================================

subroutine phys_inidat_tke_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr3d = default_value
  end if

end subroutine phys_inidat_tke_default_native

!=======================================================================

subroutine phys_inidat_kvm_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_kvm_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_KVM_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_kvm_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_kvm_default_impl = .false.
  end if

  phys_inidat_kvm_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_kvm_default_impl) then
        write(iulog,*) 'phys_inidat_kvm_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_kvm_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_kvm_default_select_impl

!=======================================================================

subroutine phys_inidat_kvm_default(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pverp_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_kvm_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, tptr3d_p, default_value) &
          bind(c, name="phys_inidat_batch_kvm_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pverp_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr3d_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_kvm_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_kvm_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pverp_c = int(pverp_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_kvm_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, c_loc(tptr3d), default_value)

end subroutine phys_inidat_kvm_default

!=======================================================================

subroutine phys_inidat_kvm_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr3d = default_value
  end if

end subroutine phys_inidat_kvm_default_native

!=======================================================================

subroutine phys_inidat_kvh_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_kvh_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_KVH_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_kvh_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_kvh_default_impl = .false.
  end if

  phys_inidat_kvh_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_kvh_default_impl) then
        write(iulog,*) 'phys_inidat_kvh_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_kvh_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_kvh_default_select_impl

!=======================================================================

subroutine phys_inidat_kvh_default(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pverp_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_kvh_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, tptr3d_p, default_value) &
          bind(c, name="phys_inidat_batch_kvh_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pverp_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr3d_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_kvh_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_kvh_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pverp_c = int(pverp_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_kvh_default_codon(pcols_c, pverp_c, chunk_count_c, found_c, c_loc(tptr3d), default_value)

end subroutine phys_inidat_kvh_default

!=======================================================================

subroutine phys_inidat_kvh_default_native(pcols_local, pverp_local, chunk_count_local, found_local, tptr3d, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pverp_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr3d(pcols_local, pverp_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr3d = default_value
  end if

end subroutine phys_inidat_kvh_default_native

!=======================================================================

subroutine phys_inidat_qcwat_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_qcwat_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_QCWAT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_qcwat_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_qcwat_default_impl = .false.
  end if

  phys_inidat_qcwat_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_qcwat_default_impl) then
        write(iulog,*) 'phys_inidat_qcwat_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_qcwat_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_qcwat_default_select_impl

!=======================================================================

subroutine phys_inidat_qcwat_default(primary_found_local, fallback_found_local, init_source_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: fallback_found_local
  integer, intent(out) :: init_source_local

  integer(c_int64_t) :: primary_found_c
  integer(c_int64_t) :: fallback_found_c
  integer(c_int64_t), target :: init_source_c

  interface
     subroutine phys_inidat_batch_qcwat_default_codon(primary_found_c, fallback_found_c, init_source_p) &
          bind(c, name="phys_inidat_batch_qcwat_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: primary_found_c
       integer(c_int64_t), value :: fallback_found_c
       type(c_ptr), value :: init_source_p
     end subroutine phys_inidat_batch_qcwat_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_qcwat_default_native(primary_found_local, fallback_found_local, init_source_local)
     return
  end if

  primary_found_c = 0_c_int64_t
  if (primary_found_local) primary_found_c = 1_c_int64_t
  fallback_found_c = 0_c_int64_t
  if (fallback_found_local) fallback_found_c = 1_c_int64_t
  init_source_c = 0_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_qcwat_default_codon(primary_found_c, fallback_found_c, c_loc(init_source_c))
  init_source_local = int(init_source_c)

end subroutine phys_inidat_qcwat_default

!=======================================================================

subroutine phys_inidat_qcwat_default_native(primary_found_local, fallback_found_local, init_source_local)

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: fallback_found_local
  integer, intent(out) :: init_source_local

  if (primary_found_local) then
     init_source_local = 1
  else if (fallback_found_local) then
     init_source_local = 2
  else
     init_source_local = 0
  end if

end subroutine phys_inidat_qcwat_default_native

!=======================================================================

subroutine phys_inidat_iccwat_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_iccwat_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_ICCWAT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_iccwat_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_iccwat_default_impl = .false.
  end if

  phys_inidat_iccwat_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_iccwat_default_impl) then
        write(iulog,*) 'phys_inidat_iccwat_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_iccwat_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_iccwat_default_select_impl

!=======================================================================

subroutine phys_inidat_iccwat_default(primary_found_local, fallback_found_local, init_source_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: fallback_found_local
  integer, intent(out) :: init_source_local

  integer(c_int64_t) :: primary_found_c
  integer(c_int64_t) :: fallback_found_c
  integer(c_int64_t), target :: init_source_c

  interface
     subroutine phys_inidat_batch_iccwat_default_codon(primary_found_c, fallback_found_c, init_source_p) &
          bind(c, name="phys_inidat_batch_iccwat_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: primary_found_c
       integer(c_int64_t), value :: fallback_found_c
       type(c_ptr), value :: init_source_p
     end subroutine phys_inidat_batch_iccwat_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_iccwat_default_native(primary_found_local, fallback_found_local, init_source_local)
     return
  end if

  primary_found_c = 0_c_int64_t
  if (primary_found_local) primary_found_c = 1_c_int64_t
  fallback_found_c = 0_c_int64_t
  if (fallback_found_local) fallback_found_c = 1_c_int64_t
  init_source_c = 0_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_iccwat_default_codon(primary_found_c, fallback_found_c, c_loc(init_source_c))
  init_source_local = int(init_source_c)

end subroutine phys_inidat_iccwat_default

!=======================================================================

subroutine phys_inidat_iccwat_default_native(primary_found_local, fallback_found_local, init_source_local)

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: fallback_found_local
  integer, intent(out) :: init_source_local

  if (primary_found_local) then
     init_source_local = 1
  else if (fallback_found_local) then
     init_source_local = 2
  else
     init_source_local = 0
  end if

end subroutine phys_inidat_iccwat_default_native

!=======================================================================

subroutine phys_inidat_lcwat_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_lcwat_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_LCWAT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_lcwat_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_lcwat_default_impl = .false.
  end if

  phys_inidat_lcwat_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_lcwat_default_impl) then
        write(iulog,*) 'phys_inidat_lcwat_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_lcwat_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_lcwat_default_select_impl

!=======================================================================

subroutine phys_inidat_lcwat_default(primary_found_local, cldice_found_local, cldliq_found_local, init_source_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: cldice_found_local
  logical, intent(in) :: cldliq_found_local
  integer, intent(out) :: init_source_local

  integer(c_int64_t) :: primary_found_c
  integer(c_int64_t) :: cldice_found_c
  integer(c_int64_t) :: cldliq_found_c
  integer(c_int64_t), target :: init_source_c

  interface
     subroutine phys_inidat_batch_lcwat_default_codon(primary_found_c, cldice_found_c, cldliq_found_c, init_source_p) &
          bind(c, name="phys_inidat_batch_lcwat_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: primary_found_c
       integer(c_int64_t), value :: cldice_found_c
       integer(c_int64_t), value :: cldliq_found_c
       type(c_ptr), value :: init_source_p
     end subroutine phys_inidat_batch_lcwat_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_lcwat_default_native(primary_found_local, cldice_found_local, cldliq_found_local, init_source_local)
     return
  end if

  primary_found_c = 0_c_int64_t
  if (primary_found_local) primary_found_c = 1_c_int64_t
  cldice_found_c = 0_c_int64_t
  if (cldice_found_local) cldice_found_c = 1_c_int64_t
  cldliq_found_c = 0_c_int64_t
  if (cldliq_found_local) cldliq_found_c = 1_c_int64_t
  init_source_c = 0_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_lcwat_default_codon(primary_found_c, cldice_found_c, cldliq_found_c, c_loc(init_source_c))
  init_source_local = int(init_source_c)

end subroutine phys_inidat_lcwat_default

!=======================================================================

subroutine phys_inidat_lcwat_default_native(primary_found_local, cldice_found_local, cldliq_found_local, init_source_local)

  logical, intent(in) :: primary_found_local
  logical, intent(in) :: cldice_found_local
  logical, intent(in) :: cldliq_found_local
  integer, intent(out) :: init_source_local

  if (primary_found_local) then
     init_source_local = 1
  else if (cldice_found_local .and. cldliq_found_local) then
     init_source_local = 2
  else if (cldice_found_local) then
     init_source_local = 3
  else if (cldliq_found_local) then
     init_source_local = 4
  else
     init_source_local = 0
  end if

end subroutine phys_inidat_lcwat_default_native

!=======================================================================

subroutine phys_inidat_tcwat_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_tcwat_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_TCWAT_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_tcwat_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_tcwat_default_impl = .false.
  end if

  phys_inidat_tcwat_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_tcwat_default_impl) then
        write(iulog,*) 'phys_inidat_tcwat_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_tcwat_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_tcwat_default_select_impl

!=======================================================================

subroutine phys_inidat_tcwat_default(primary_found_local, init_source_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: primary_found_local
  integer, intent(out) :: init_source_local

  integer(c_int64_t) :: primary_found_c
  integer(c_int64_t), target :: init_source_c

  interface
     subroutine phys_inidat_batch_tcwat_default_codon(primary_found_c, init_source_p) &
          bind(c, name="phys_inidat_batch_tcwat_default_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: primary_found_c
       type(c_ptr), value :: init_source_p
     end subroutine phys_inidat_batch_tcwat_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_tcwat_default_native(primary_found_local, init_source_local)
     return
  end if

  primary_found_c = 0_c_int64_t
  if (primary_found_local) primary_found_c = 1_c_int64_t
  init_source_c = 0_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_tcwat_default_codon(primary_found_c, c_loc(init_source_c))
  init_source_local = int(init_source_c)

end subroutine phys_inidat_tcwat_default

!=======================================================================

subroutine phys_inidat_tcwat_default_native(primary_found_local, init_source_local)

  logical, intent(in) :: primary_found_local
  integer, intent(out) :: init_source_local

  if (primary_found_local) then
     init_source_local = 1
  else
     init_source_local = 2
  end if

end subroutine phys_inidat_tcwat_default_native

!=======================================================================

subroutine phys_inidat_cloud_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_cloud_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_CLOUD_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_cloud_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_cloud_default_impl = .false.
  end if

  phys_inidat_cloud_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_cloud_default_impl) then
        write(iulog,*) 'phys_inidat_cloud_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_cloud_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_cloud_default_select_impl

!=======================================================================

subroutine phys_inidat_cloud_default(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr3d(pcols_local, pver_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_cloud_default_codon(pcols_c, pver_c, chunk_count_c, found_c, tptr3d_p, default_value) &
          bind(c, name="phys_inidat_batch_cloud_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr3d_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_cloud_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_cloud_default_native(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_cloud_default_codon(pcols_c, pver_c, chunk_count_c, found_c, c_loc(tptr3d), default_value)

end subroutine phys_inidat_cloud_default

!=======================================================================

subroutine phys_inidat_cloud_default_native(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr3d(pcols_local, pver_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr3d = default_value
  end if

end subroutine phys_inidat_cloud_default_native

!=======================================================================

subroutine phys_inidat_concld_default_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_concld_default_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_CONCLD_DEFAULT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_concld_default_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_concld_default_impl = .false.
  end if

  phys_inidat_concld_default_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_concld_default_impl) then
        write(iulog,*) 'phys_inidat_concld_default implementation = native'
     else
        write(iulog,*) 'phys_inidat_concld_default implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_concld_default_select_impl

!=======================================================================

subroutine phys_inidat_concld_default(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout), target :: tptr3d(pcols_local, pver_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: chunk_count_c
  integer(c_int64_t) :: found_c

  interface
     subroutine phys_inidat_batch_concld_default_codon(pcols_c, pver_c, chunk_count_c, found_c, tptr3d_p, default_value) &
          bind(c, name="phys_inidat_batch_concld_default_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: chunk_count_c
       integer(c_int64_t), value :: found_c
       type(c_ptr), value :: tptr3d_p
       real(c_double), value :: default_value
     end subroutine phys_inidat_batch_concld_default_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_concld_default_native(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  chunk_count_c = int(chunk_count_local, c_int64_t)
  found_c = 0_c_int64_t
  if (found_local) found_c = 1_c_int64_t

  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_concld_default_codon(pcols_c, pver_c, chunk_count_c, found_c, c_loc(tptr3d), default_value)

end subroutine phys_inidat_concld_default

!=======================================================================

subroutine phys_inidat_concld_default_native(pcols_local, pver_local, chunk_count_local, found_local, tptr3d, default_value)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: chunk_count_local
  logical, intent(in) :: found_local
  real(r8), intent(inout) :: tptr3d(pcols_local, pver_local, chunk_count_local)
  real(r8), intent(in) :: default_value

  if (.not. found_local) then
     tptr3d = default_value
  end if

end subroutine phys_inidat_concld_default_native

!=======================================================================

subroutine phys_inidat_tbot_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (phys_inidat_tbot_init_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('PHYS_INIDAT_TBOT_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_phys_inidat_tbot_init_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_phys_inidat_tbot_init_impl = .false.
  end if

  phys_inidat_tbot_init_impl_selected = .true.

  if (masterproc) then
     if (use_native_phys_inidat_tbot_init_impl) then
        write(iulog,*) 'phys_inidat_tbot_init implementation = native'
     else
        write(iulog,*) 'phys_inidat_tbot_init implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine phys_inidat_tbot_init_select_impl

!=======================================================================

subroutine phys_inidat_tbot_init(pcols_local, tbot, posinf_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  real(r8), intent(inout), target :: tbot(pcols_local)
  real(r8), intent(in) :: posinf_local

  integer(c_int64_t) :: pcols_c

  interface
     subroutine phys_inidat_batch_tbot_init_codon(pcols_c, tbot_p, posinf_local) bind(c, name="phys_inidat_batch_tbot_init_codon")
       use iso_c_binding, only: c_int64_t, c_double, c_ptr
       integer(c_int64_t), value :: pcols_c
       type(c_ptr), value :: tbot_p
       real(c_double), value :: posinf_local
     end subroutine phys_inidat_batch_tbot_init_codon
  end interface

  call phys_inidat_batch_select_impl()

  if (use_native_phys_inidat_batch_impl) then
     call phys_inidat_tbot_init_native(pcols_local, tbot, posinf_local)
     return
  end if

  pcols_c = int(pcols_local, c_int64_t)
  call phys_inidat_batch_log_entered()
  call phys_inidat_batch_tbot_init_codon(pcols_c, c_loc(tbot), posinf_local)

end subroutine phys_inidat_tbot_init

!=======================================================================

subroutine phys_inidat_tbot_init_native(pcols_local, tbot, posinf_local)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: pcols_local
  real(r8), intent(inout) :: tbot(pcols_local)
  real(r8), intent(in) :: posinf_local

  tbot(:) = posinf_local

end subroutine phys_inidat_tbot_init_native

!=======================================================================

subroutine tphysbc_init_fields_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tphysbc_init_fields_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TPHYSBC_INIT_FIELDS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tphysbc_init_fields_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tphysbc_init_fields_impl = .false.
  end if

  tphysbc_init_fields_impl_selected = .true.

  if (masterproc) then
     if (use_native_tphysbc_init_fields_impl) then
        write(iulog,*) 'tphysbc_init_fields implementation = native'
     else
        write(iulog,*) 'tphysbc_init_fields implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine tphysbc_init_fields_select_impl

!=======================================================================

subroutine tphysbc_init_fields(ncol, pcols_local, pver_local, pcnst_local, fracis, tend_dtdt, tend_dudt, tend_dvdt)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(inout), target :: fracis(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout), target :: tend_dtdt(pcols_local, pver_local)
  real(r8), intent(inout), target :: tend_dudt(pcols_local, pver_local)
  real(r8), intent(inout), target :: tend_dvdt(pcols_local, pver_local)

  integer(c_int64_t) :: ncol_c
  integer(c_int64_t) :: pcols_c
  integer(c_int64_t) :: pver_c
  integer(c_int64_t) :: pcnst_c

  interface
     subroutine tphysbc_init_fields_codon(ncol_c, pcols_c, pver_c, pcnst_c, fracis_p, tend_dtdt_p, tend_dudt_p, tend_dvdt_p) &
          bind(c, name="tphysbc_init_fields_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       integer(c_int64_t), value :: pcols_c
       integer(c_int64_t), value :: pver_c
       integer(c_int64_t), value :: pcnst_c
       type(c_ptr), value :: fracis_p
       type(c_ptr), value :: tend_dtdt_p
       type(c_ptr), value :: tend_dudt_p
       type(c_ptr), value :: tend_dvdt_p
     end subroutine tphysbc_init_fields_codon
  end interface

  call tphysbc_init_fields_select_impl()

  if (use_native_tphysbc_init_fields_impl) then
     call tphysbc_init_fields_native(ncol, pcols_local, pver_local, pcnst_local, fracis, tend_dtdt, tend_dudt, tend_dvdt)
     return
  end if

  ncol_c = int(ncol, c_int64_t)
  pcols_c = int(pcols_local, c_int64_t)
  pver_c = int(pver_local, c_int64_t)
  pcnst_c = int(pcnst_local, c_int64_t)

  call tphysbc_init_fields_codon(ncol_c, pcols_c, pver_c, pcnst_c, &
       c_loc(fracis), c_loc(tend_dtdt), c_loc(tend_dudt), c_loc(tend_dvdt))

end subroutine tphysbc_init_fields

!=======================================================================

subroutine tphysbc_init_fields_native(ncol, pcols_local, pver_local, pcnst_local, fracis, tend_dtdt, tend_dudt, tend_dvdt)

  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: ncol
  integer, intent(in) :: pcols_local
  integer, intent(in) :: pver_local
  integer, intent(in) :: pcnst_local
  real(r8), intent(inout) :: fracis(pcols_local, pver_local, pcnst_local)
  real(r8), intent(inout) :: tend_dtdt(pcols_local, pver_local)
  real(r8), intent(inout) :: tend_dudt(pcols_local, pver_local)
  real(r8), intent(inout) :: tend_dvdt(pcols_local, pver_local)

  fracis(:ncol,:,1:pcnst_local) = 1._r8
  tend_dtdt(:ncol,:pver_local) = 0._r8
  tend_dudt(:ncol,:pver_local) = 0._r8
  tend_dvdt(:ncol,:pver_local) = 0._r8

end subroutine tphysbc_init_fields_native

end module physpkg

module cldwat2m_macro

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ap_mmacro_pcond_scheme, only: &
       mmacro_pcond_init, &
       mmacro_pcond => mmacro_pcond_run
  use cloud_microphysics_host_hooks, only: &
       register_cloud_microphysics_host_hooks

  implicit none
  private

  public :: ini_macro
  public :: mmacro_pcond

contains

  subroutine ini_macro(rhminl_opt, rhmini_opt)
    use cam_abortutils, only: endrun
    use cam_history, only: addfld, phys_decomp
    use cam_logfile, only: iulog
    use cldfrc2m, only: rhmini_const, rhmaxi_const
    use cloud_fraction, only: cldfrc_getparams
    use constituents, only: qmin, cnst_get_ind
    use physconst, only: cpair, epsilo, gravit, h2otrip, latice, latvap, &
         rair, rh2o, tmelt
    use ppgrid, only: pcols, pver, pverp
    use ref_pres, only: top_lev => trop_cloud_top_lev
    use spmd_utils, only: masterproc

    integer, intent(in) :: rhminl_opt, rhmini_opt

    character(len=256) :: errmsg
    integer :: errflg, iceopt, ixcldice, ixcldliq
    real(r8) :: icecrit, premib, premit
    real(r8) :: rhminh, rhminl, rhminl_adj_land

    call cldfrc_getparams(rhminl_out=rhminl, &
         rhminl_adj_land_out=rhminl_adj_land, rhminh_out=rhminh, &
         premit_out=premit, premib_out=premib, iceopt_out=iceopt, &
         icecrit_out=icecrit)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)

    call register_cloud_microphysics_host_hooks(cam_timer_start, &
         cam_timer_stop, cam_outfld_real2d, cam_history_active, &
         cam_endrun)

    call mmacro_pcond_init(rhminl_opt, rhmini_opt, pcols, pver, pverp, &
         top_lev, iulog, cpair, latvap, latice, rh2o, gravit, rair, &
         qmin(1), qmin(ixcldliq), qmin(ixcldice), rhmini_const, &
         rhmaxi_const, rhminl, rhminl_adj_land, rhminh, premit, premib, &
         iceopt, icecrit, epsilo, tmelt, h2otrip, errmsg, errflg)
    if (errflg /= 0) call endrun(trim(errmsg))

    if (masterproc) then
       write(iulog,*) 'Park Macrophysics Parameters'
       write(iulog,*) '  rhminl          = ', rhminl
       write(iulog,*) '  rhminl_adj_land = ', rhminl_adj_land
       write(iulog,*) '  rhminh          = ', rhminh
       write(iulog,*) '  premit          = ', premit
       write(iulog,*) '  premib          = ', premib
       write(iulog,*) '  i_rhminl        = ', rhminl_opt
       write(iulog,*) '  i_rhmini        = ', rhmini_opt
    end if

    call addfld('RHMIN_LIQ', 'fraction', pver, 'A', &
         'Default critical RH for liquid-stratus', phys_decomp)
    call addfld('RHMIN_ICE', 'fraction', pver, 'A', &
         'Default critical RH for    ice-stratus', phys_decomp)
    call addfld('DRHMINPBL_LIQ', 'fraction', pver, 'A', &
         'Drop of liquid-stratus critical RH by PBL turbulence', phys_decomp)
    call addfld('DRHMINPBL_ICE', 'fraction', pver, 'A', &
         'Drop of    ice-stratus critical RH by PBL turbulence', phys_decomp)
    call addfld('DRHMINDET_LIQ', 'fraction', pver, 'A', &
         'Drop of liquid-stratus critical RH by convective detrainment', phys_decomp)
    call addfld('DRHMINDET_ICE', 'fraction', pver, 'A', &
         'Drop of    ice-stratus critical RH by convective detrainment', phys_decomp)
  end subroutine ini_macro

  subroutine cam_timer_start(timer_name)
    use perf_mod, only: t_startf
    character(len=*), intent(in) :: timer_name

    call t_startf(timer_name)
  end subroutine cam_timer_start

  subroutine cam_timer_stop(timer_name)
    use perf_mod, only: t_stopf
    character(len=*), intent(in) :: timer_name

    call t_stopf(timer_name)
  end subroutine cam_timer_stop

  subroutine cam_outfld_real2d(field_name, field, dim1, lchnk)
    use cam_history, only: outfld
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: dim1, lchnk

    call outfld(field_name, field, dim1, lchnk)
  end subroutine cam_outfld_real2d

  logical function cam_history_active(field_name)
    use cam_history, only: hist_fld_active
    character(len=*), intent(in) :: field_name

    cam_history_active = hist_fld_active(field_name)
  end function cam_history_active

  subroutine cam_endrun(message)
    use cam_abortutils, only: endrun
    character(len=*), intent(in) :: message

    call endrun(message)
  end subroutine cam_endrun

end module cldwat2m_macro

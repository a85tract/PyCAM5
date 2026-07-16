module uwshcu

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ap_uwshcu_processes_scheme, only: scheme_set_rpen => uwshcu_set_rpen, &
       scheme_init_uwshcu => init_uwshcu, &
       compute_uwshcu_run, &
       compute_uwshcu_inv => compute_uwshcu_inv_run
  use ap_water_isotope_fractionation, only: &
       register_water_isotope_alpha_hook
  use ap_saturation_table, only: register_qsat_table_hook, &
       register_findsp_table_hook
  use atmos_phys_history_hooks, only: register_atmos_phys_history_hooks
  use cam_abortutils, only: endrun
  use cam_history, only: addfld, outfld, phys_decomp
  use perf_mod, only: t_startf, t_stopf
  use ppgrid, only: pver, pverp
  use physconst, only: epsilo, rh2o, rhoh2o, tmelt, h2otrip
  use wv_saturation, only: wv_sat_export_table, qsat, findsp_vc
  use constituents, only: pcnst, qmin, cnst_get_ind, cnst_get_type_byind
  use water_tracer_vars, only: trace_water, wisotope, wtrc_nwset, &
       wtrc_iatype, iwspec, wtrc_qmin, wtrc_alpha_kinetic, &
       wtrc_fixed_alpha
  use water_tracers, only: wtrc_is_wtrc, wtrc_get_rstd
  use water_types, only: iwtvap, iwtliq, iwtice, iwtstrain, iwtcvrain, &
       wtype_get_alpha
  use water_isotopes, only: pwtspec, isph2o, isphdo

  implicit none
  private

  public :: uwshcu_readnl
  public :: init_uwshcu
  public :: compute_uwshcu
  public :: compute_uwshcu_inv

contains

  subroutine uwshcu_readnl(nlfile)
    use namelist_utils, only: find_group_name
    use units, only: getunit, freeunit
    use spmd_utils, only: masterproc
    use mpishorthand

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr
    real(r8) :: uwshcu_rpen = huge(1._r8)
    namelist /uwshcu_nl/ uwshcu_rpen

    if (masterproc) then
       unitn = getunit()
       open(unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'uwshcu_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, uwshcu_nl, iostat=ierr)
          if (ierr /= 0) call endrun('uwshcu_readnl:: ERROR reading namelist')
       end if
       close(unitn)
       call freeunit(unitn)
    end if
#ifdef SPMD
    call mpibcast(uwshcu_rpen, 1, mpir8, 0, mpicom)
#endif
    call scheme_set_rpen(uwshcu_rpen)
  end subroutine uwshcu_readnl

  subroutine init_uwshcu(kind, xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in)
    integer, intent(in) :: kind
    real(r8), intent(in) :: xlv_in, cp_in, xlf_in, zvir_in, r_in, g_in, ep2_in
    real(r8), allocatable :: sat_table(:)
    logical :: constituent_is_wet(pcnst), water_tracer_mask(pcnst)
    real(r8) :: species_rstd(pwtspec)
    integer :: ixnumliq, ixnumice, ixcldliq, ixcldice, m
    character(len=512) :: errmsg
    integer :: errflg

    call cnst_get_ind('NUMLIQ', ixnumliq)
    call cnst_get_ind('NUMICE', ixnumice)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    do m = 1, pcnst
       constituent_is_wet(m) = cnst_get_type_byind(m) == 'wet'
       water_tracer_mask(m) = wtrc_is_wtrc(m)
    end do
    do m = 1, pwtspec
       species_rstd(m) = wtrc_get_rstd(m)
    end do
    call register_water_isotope_alpha_hook(wtype_get_alpha)
    call register_qsat_table_hook(cam_qsat)
    call register_findsp_table_hook(cam_findsp)
    call register_atmos_phys_history_hooks(cam_outfld_1d, cam_outfld_2d)
    call register_uwshcu_history_fields()
    call wv_sat_export_table(sat_table)
    call scheme_init_uwshcu(kind, xlv_in, cp_in, xlf_in, zvir_in, r_in, &
         g_in, ep2_in, sat_table, epsilo, rh2o, tmelt, qmin, &
         constituent_is_wet, ixnumliq, ixnumice, ixcldliq, ixcldice, &
         trace_water, wisotope, wtrc_nwset, wtrc_iatype, iwspec, &
         water_tracer_mask, wtrc_qmin, iwtvap, iwtliq, iwtice, &
         iwtstrain, iwtcvrain, isph2o, isphdo, h2otrip, rhoh2o, &
         wtrc_alpha_kinetic, wtrc_fixed_alpha, species_rstd, errmsg, errflg)
    if (errflg /= 0) call endrun(trim(errmsg))
  end subroutine init_uwshcu

  subroutine compute_uwshcu(mix, mkx, iend, ncnst, dt, ps0_in, zs0_in, &
       p0_in, z0_in, dp0_in, u0_in, v0_in, qv0_in, ql0_in, qi0_in, t0_in, &
       s0_in, tr0_in, tke_in, cldfrct_in, concldfrct_in, pblh_in, &
       cush_inout, umf_out, slflx_out, qtflx_out, flxprc1_out, flxsnow1_out, &
       qvten_out, qlten_out, qiten_out, sten_out, uten_out, vten_out, &
       trten_out, qrten_out, qsten_out, precip_out, snow_out, evapc_out, &
       cufrc_out, qcu_out, qlu_out, qiu_out, cbmf_out, qc_out, rliq_out, &
       cnt_out, cnb_out, lchnk, dpdry0_in, wtprec_out, wtsnow_out, wtqc_out)

    integer, intent(in) :: mix, mkx, iend, ncnst, lchnk
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: ps0_in(:,0:), zs0_in(:,0:), p0_in(:,:), z0_in(:,:)
    real(r8), intent(in) :: dp0_in(:,:), u0_in(:,:), v0_in(:,:), qv0_in(:,:)
    real(r8), intent(in) :: ql0_in(:,:), qi0_in(:,:), t0_in(:,:), s0_in(:,:)
    real(r8), intent(in) :: tr0_in(:,:,:), tke_in(:,0:), cldfrct_in(:,:)
    real(r8), intent(in) :: concldfrct_in(:,:), pblh_in(:), dpdry0_in(:,:)
    real(r8), intent(inout) :: cush_inout(:)
    real(r8), intent(out) :: umf_out(:,0:), slflx_out(:,0:), qtflx_out(:,0:)
    real(r8), intent(out) :: flxprc1_out(:,0:), flxsnow1_out(:,0:)
    real(r8), intent(out) :: qvten_out(:,:), qlten_out(:,:), qiten_out(:,:)
    real(r8), intent(out) :: sten_out(:,:), uten_out(:,:), vten_out(:,:)
    real(r8), intent(out) :: trten_out(:,:,:), qrten_out(:,:), qsten_out(:,:)
    real(r8), intent(out) :: precip_out(:), snow_out(:), evapc_out(:,:)
    real(r8), intent(out) :: cufrc_out(:,:), qcu_out(:,:), qlu_out(:,:), qiu_out(:,:)
    real(r8), intent(out) :: cbmf_out(:), qc_out(:,:), rliq_out(:), cnt_out(:), cnb_out(:)
    real(r8), intent(out) :: wtprec_out(:,:), wtsnow_out(:,:), wtqc_out(:,:,:)
    character(len=512) :: errmsg
    integer :: errflg

    call t_startf('ap_compute_uwshcu_run')
    call compute_uwshcu_run(mix, mkx, iend, ncnst, dt, ps0_in, zs0_in, &
         p0_in, z0_in, dp0_in, u0_in, v0_in, qv0_in, ql0_in, qi0_in, t0_in, &
         s0_in, tr0_in, tke_in, cldfrct_in, concldfrct_in, pblh_in, &
         cush_inout, umf_out, slflx_out, qtflx_out, flxprc1_out, flxsnow1_out, &
         qvten_out, qlten_out, qiten_out, sten_out, uten_out, vten_out, &
         trten_out, qrten_out, qsten_out, precip_out, snow_out, evapc_out, &
         cufrc_out, qcu_out, qlu_out, qiu_out, cbmf_out, qc_out, rliq_out, &
         cnt_out, cnb_out, lchnk, dpdry0_in, wtprec_out, wtsnow_out, wtqc_out, &
         errmsg, errflg)
    call t_stopf('ap_compute_uwshcu_run')
    if (errflg /= 0) call endrun(trim(errmsg))
  end subroutine compute_uwshcu

  subroutine register_uwshcu_history_fields()
    call addfld( 'qtflx_Cu'       , 'kg/m2/s' , pverp , 'A' , 'Convective qt flux'                                  , phys_decomp )
    call addfld( 'slflx_Cu'       , 'J/m2/s'  , pverp , 'A' , 'Convective sl flux'                                  , phys_decomp )
    call addfld( 'uflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  u flux'                                  , phys_decomp )
    call addfld( 'vflx_Cu'        , 'kg/m/s2' , pverp , 'A' , 'Convective  v flux'                                  , phys_decomp )
    call addfld( 'qtten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qt tendency by convection'                           , phys_decomp )
    call addfld( 'slten_Cu'       , 'J/kg/s'  , pver  , 'A' , 'sl tendency by convection'                           , phys_decomp )
    call addfld( 'uten_Cu'        , 'm/s2'    , pver  , 'A' , ' u tendency by convection'                           , phys_decomp )
    call addfld( 'vten_Cu'        , 'm/s2'    , pver  , 'A' , ' v tendency by convection'                           , phys_decomp )
    call addfld( 'qvten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qv tendency by convection'                           , phys_decomp )
    call addfld( 'qlten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'ql tendency by convection'                           , phys_decomp )
    call addfld( 'qiten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'qi tendency by convection'                           , phys_decomp )
    call addfld( 'cbmf_Cu'        , 'kg/m2/s' , 1     , 'A' , 'Cumulus base mass flux'                              , phys_decomp )
    call addfld( 'ufrcinvbase_Cu' , 'fraction', 1     , 'A' , 'Cumulus fraction at PBL top'                         , phys_decomp )
    call addfld( 'ufrclcl_Cu'     , 'fraction', 1     , 'A' , 'Cumulus fraction at LCL'                             , phys_decomp )
    call addfld( 'winvbase_Cu'    , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at PBL top'                , phys_decomp )
    call addfld( 'wlcl_Cu'        , 'm/s'     , 1     , 'A' , 'Cumulus vertical velocity at LCL'                    , phys_decomp )
    call addfld( 'plcl_Cu'        , 'Pa'      , 1     , 'A' , 'LCL of source air'                                   , phys_decomp )
    call addfld( 'pinv_Cu'        , 'Pa'      , 1     , 'A' , 'PBL top pressure'                                    , phys_decomp )
    call addfld( 'plfc_Cu'        , 'Pa'      , 1     , 'A' , 'LFC of source air'                                   , phys_decomp )
    call addfld( 'pbup_Cu'        , 'Pa'      , 1     , 'A' , 'Highest interface level of positive cumulus buoyancy', phys_decomp )
    call addfld( 'ppen_Cu'        , 'Pa'      , 1     , 'A' , 'Highest level where cumulus w is 0'                  , phys_decomp )
    call addfld( 'qtsrc_Cu'       , 'kg/kg'   , 1     , 'A' , 'Cumulus source air qt'                               , phys_decomp )
    call addfld( 'thlsrc_Cu'      , 'K'       , 1     , 'A' , 'Cumulus source air thl'                              , phys_decomp )
    call addfld( 'thvlsrc_Cu'     , 'K'       , 1     , 'A' , 'Cumulus source air thvl'                             , phys_decomp )
    call addfld( 'emfkbup_Cu'     , 'kg/m2/s' , 1     , 'A' , 'Penetrative mass flux at kbup'                       , phys_decomp )
    call addfld( 'cin_Cu'         , 'J/kg'    , 1     , 'A' , 'CIN upto LFC'                                        , phys_decomp )
    call addfld( 'cinlcl_Cu'      , 'J/kg'    , 1     , 'A' , 'CIN upto LCL'                                        , phys_decomp )
    call addfld( 'cbmflimit_Cu'   , 'kg/m2/s' , 1     , 'A' , 'cbmflimiter'                                         , phys_decomp )
    call addfld( 'tkeavg_Cu'      , 'm2/s2'   , 1     , 'A' , 'Average tke within PBL for convection scheme'        , phys_decomp )
    call addfld( 'zinv_Cu'        , 'm'       , 1     , 'A' , 'PBL top height'                                      , phys_decomp )
    call addfld( 'rcwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP+IWP'                                     , phys_decomp )
    call addfld( 'rlwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus LWP'                                         , phys_decomp )
    call addfld( 'riwp_Cu'        , 'kg/m2'   , 1     , 'A' , 'Cumulus IWP'                                         , phys_decomp )
    call addfld( 'tophgt_Cu'      , 'm'       , 1     , 'A' , 'Cumulus top height'                                  , phys_decomp )
    call addfld( 'wu_Cu'          , 'm/s'     , pverp , 'A' , 'Convective updraft vertical velocity'                , phys_decomp )
    call addfld( 'ufrc_Cu'        , 'fraction', pverp , 'A' , 'Convective updraft fractional area'                  , phys_decomp )
    call addfld( 'qtu_Cu'         , 'kg/kg'   , pverp , 'A' , 'Cumulus updraft qt'                                  , phys_decomp )
    call addfld( 'thlu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thl'                                 , phys_decomp )
    call addfld( 'thvu_Cu'        , 'K'       , pverp , 'A' , 'Cumulus updraft thv'                                 , phys_decomp )
    call addfld( 'uu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft uwnd'                                , phys_decomp )
    call addfld( 'vu_Cu'          , 'm/s'     , pverp , 'A' , 'Cumulus updraft vwnd'                                , phys_decomp )
    call addfld( 'qtu_emf_Cu'     , 'kg/kg'   , pverp , 'A' , 'qt of penatratively entrained air'                   , phys_decomp )
    call addfld( 'thlu_emf_Cu'    , 'K'       , pverp , 'A' , 'thl of penatratively entrained air'                  , phys_decomp )
    call addfld( 'uu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'uwnd of penatratively entrained air'                 , phys_decomp )
    call addfld( 'vu_emf_Cu'      , 'm/s'     , pverp , 'A' , 'vwnd of penatratively entrained air'                 , phys_decomp )
    call addfld( 'umf_Cu'         , 'kg/m2/s' , pverp , 'A' , 'Cumulus updraft mass flux'                           , phys_decomp )
    call addfld( 'uemf_Cu'        , 'kg/m2/s' , pverp , 'A' , 'Cumulus net ( updraft + entrainment ) mass flux'     , phys_decomp )
    call addfld( 'qcu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC+IWC'                             , phys_decomp )
    call addfld( 'qlu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft LWC'                                 , phys_decomp )
    call addfld( 'qiu_Cu'         , 'kg/kg'   , pver  , 'A' , 'Cumulus updraft IWC'                                 , phys_decomp )
    call addfld( 'cufrc_Cu'       , 'fraction', pver  , 'A' , 'Cumulus cloud fraction'                              , phys_decomp )
    call addfld( 'fer_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional entrainment rate'         , phys_decomp )
    call addfld( 'fdr_Cu'         , '1/m'     , pver  , 'A' , 'Cumulus lateral fractional detrainment Rate'         , phys_decomp )
    call addfld( 'dwten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus cloud water to env.'      , phys_decomp )
    call addfld( 'diten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Expellsion rate of cumulus ice water to env.'        , phys_decomp )
    call addfld( 'qrten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of rain by cumulus'                  , phys_decomp )
    call addfld( 'qsten_Cu'       , 'kg/kg/s' , pver  , 'A' , 'Production rate of snow by cumulus'                  , phys_decomp )
    call addfld( 'flxrain_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Rain flux induced by Cumulus'                        , phys_decomp )
    call addfld( 'flxsnow_Cu'     , 'kg/m2/s' , pverp , 'A' , 'Snow flux induced by Cumulus'                        , phys_decomp )
    call addfld( 'ntraprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of rain by Cumulus'              , phys_decomp )
    call addfld( 'ntsnprd_Cu'     , 'kg/kg/s' , pver  , 'A' , 'Net production rate of snow by Cumulus'              , phys_decomp )
    call addfld( 'excessu_Cu'     , 'no'      , pver  , 'A' , 'Updraft saturation excess'                           , phys_decomp )
    call addfld( 'excess0_Cu'     , 'no'      , pver  , 'A' , 'Environmental saturation excess'                     , phys_decomp )
    call addfld( 'xc_Cu'          , 'no'      , pver  , 'A' , 'Critical mixing ratio'                               , phys_decomp )
    call addfld( 'aquad_Cu'       , 'no'      , pver  , 'A' , 'aquad'                                               , phys_decomp )
    call addfld( 'bquad_Cu'       , 'no'      , pver  , 'A' , 'bquad'                                               , phys_decomp )
    call addfld( 'cquad_Cu'       , 'no'      , pver  , 'A' , 'cquad'                                               , phys_decomp )
    call addfld( 'bogbot_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the bottom interface'              , phys_decomp )
    call addfld( 'bogtop_Cu'      , 'no'      , pver  , 'A' , 'Cloud buoyancy at the top interface'                 , phys_decomp )
    call addfld('exit_UWCu_Cu'    , 'no'      , 1     , 'A' , 'exit_UWCu'                                           , phys_decomp )
    call addfld('exit_conden_Cu'  , 'no'      , 1     , 'A' , 'exit_conden'                                         , phys_decomp )
    call addfld('exit_klclmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klclmkx'                                        , phys_decomp )
    call addfld('exit_klfcmkx_Cu' , 'no'      , 1     , 'A' , 'exit_klfcmkx'                                        , phys_decomp )
    call addfld('exit_ufrc_Cu'    , 'no'      , 1     , 'A' , 'exit_ufrc'                                           , phys_decomp )
    call addfld('exit_wtw_Cu'     , 'no'      , 1     , 'A' , 'exit_wtw'                                            , phys_decomp )
    call addfld('exit_drycore_Cu' , 'no'      , 1     , 'A' , 'exit_drycore'                                        , phys_decomp )
    call addfld('exit_wu_Cu'      , 'no'      , 1     , 'A' , 'exit_wu'                                             , phys_decomp )
    call addfld('exit_cufilter_Cu', 'no'      , 1     , 'A' , 'exit_cufilter'                                       , phys_decomp )
    call addfld('exit_kinv1_Cu'   , 'no'      , 1     , 'A' , 'exit_kinv1'                                          , phys_decomp )
    call addfld('exit_rei_Cu'     , 'no'      , 1     , 'A' , 'exit_rei'                                            , phys_decomp )
    call addfld('limit_shcu_Cu'   , 'no'      , 1     , 'A' , 'limit_shcu'                                          , phys_decomp )
    call addfld('limit_negcon_Cu' , 'no'      , 1     , 'A' , 'limit_negcon'                                        , phys_decomp )
    call addfld('limit_ufrc_Cu'   , 'no'      , 1     , 'A' , 'limit_ufrc'                                          , phys_decomp )
    call addfld('limit_ppen_Cu'   , 'no'      , 1     , 'A' , 'limit_ppen'                                          , phys_decomp )
    call addfld('limit_emf_Cu'    , 'no'      , 1     , 'A' , 'limit_emf'                                           , phys_decomp )
    call addfld('limit_cinlcl_Cu' , 'no'      , 1     , 'A' , 'limit_cinlcl'                                        , phys_decomp )
    call addfld('limit_cin_Cu'    , 'no'      , 1     , 'A' , 'limit_cin'                                           , phys_decomp )
    call addfld('limit_cbmf_Cu'   , 'no'      , 1     , 'A' , 'limit_cbmf'                                          , phys_decomp )
    call addfld('limit_rei_Cu'    , 'no'      , 1     , 'A' , 'limit_rei'                                           , phys_decomp )
    call addfld('ind_delcin_Cu'   , 'no'      , 1     , 'A' , 'ind_delcin'                                          , phys_decomp )
  end subroutine register_uwshcu_history_fields

  subroutine cam_qsat(t, p, es, qs_out, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t, p
    real(r8), intent(out) :: es, qs_out
    real(r8), intent(out), optional :: gam, dqsdt, enthalpy

    call qsat(t, p, es, qs_out, gam, dqsdt, enthalpy)
  end subroutine cam_qsat

  subroutine cam_findsp(q, t, p, use_ice, tsp, qsp, errflg)
    real(r8), intent(in) :: q(:), t(:), p(:)
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp(:), qsp(:)
    integer, intent(out) :: errflg

    call findsp_vc(q, t, p, use_ice, tsp, qsp)
    errflg = 0
  end subroutine cam_findsp

  subroutine cam_outfld_1d(name, field, horizontal_size, chunk)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: horizontal_size, chunk

    call outfld(name, field, horizontal_size, chunk)
  end subroutine cam_outfld_1d

  subroutine cam_outfld_2d(name, field, horizontal_size, chunk)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: horizontal_size, chunk

    call outfld(name, field, horizontal_size, chunk)
  end subroutine cam_outfld_2d

end module uwshcu

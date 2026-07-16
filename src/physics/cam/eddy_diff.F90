module eddy_diff

  use shr_kind_mod, only : r8 => shr_kind_r8, i4 => shr_kind_i4
  use ap_compute_eddy_diff_scheme, only : scheme_init_eddy_diff => init_eddy_diff, &
       compute_eddy_diff_run
  use ap_saturation_table, only : register_qsat_table_hook
  use atmos_phys_history_hooks, only : register_atmos_phys_history_hooks
  use atmos_phys_log_hooks, only : register_atmos_phys_log_hook
  use perf_mod, only : t_startf, t_stopf
  use cam_abortutils, only : endrun
  use cam_history, only : addfld, outfld, phys_decomp
  use cam_logfile, only : iulog
  use physconst, only : epsilo, rh2o, tmelt, cpairv, rairv, rair
  use phys_control, only : waccmx_is
  use ref_pres, only : pref_mid
  use spmd_utils, only : masterproc
  use wv_saturation, only : wv_sat_export_table, qsat

  implicit none
  private

  public :: init_eddy_diff
  public :: compute_eddy_diff

contains

  subroutine init_eddy_diff(kind, pver, gravx, cpairx, rairx, zvirx, &
       latvapx, laticex, ntop_eddy, nbot_eddy, vkx, eddy_lbulk_max, &
       eddy_leng_max, eddy_max_bot_pressure, eddy_moist_entrain_a2l)

    integer, intent(in) :: kind, pver, ntop_eddy, nbot_eddy
    real(r8), intent(in) :: gravx, cpairx, rairx, zvirx
    real(r8), intent(in) :: latvapx, laticex, vkx
    real(r8), intent(in) :: eddy_lbulk_max, eddy_leng_max
    real(r8), intent(in) :: eddy_max_bot_pressure, eddy_moist_entrain_a2l

    real(r8), allocatable :: sat_table(:)
    character(len=512) :: errmsg
    integer :: errflg, k
    real(r8) :: leng_max_value

    if (kind /= r8) then
       write(iulog,*) 'wrong KIND of reals passed to init_diffusvity -- exiting.'
    end if
    call wv_sat_export_table(sat_table)
    call scheme_init_eddy_diff(kind, pver, gravx, cpairx, rairx, zvirx, &
         latvapx, laticex, ntop_eddy, nbot_eddy, vkx, eddy_lbulk_max, &
         eddy_leng_max, eddy_max_bot_pressure, eddy_moist_entrain_a2l, &
         pref_mid(1:pver), sat_table, epsilo, rh2o, tmelt, errmsg, errflg)
    if (errflg /= 0) call endrun(trim(errmsg))

    if (masterproc) then
       write(iulog,*) 'init_eddy_diff: eddy_leng_max=', eddy_leng_max, &
            ' lbulk_max=', eddy_lbulk_max
       do k = 1, pver
          leng_max_value = 40.e3_r8
          if (pref_mid(k) <= eddy_max_bot_pressure*1.D2) then
             leng_max_value = eddy_leng_max
          end if
          write(iulog,*) 'init_eddy_diff:', k, pref_mid(k), &
               'leng_max=', leng_max_value
       end do
    end if

    call register_qsat_table_hook(cam_qsat)
    call register_atmos_phys_history_hooks(cam_outfld_1d, cam_outfld_2d)
    call register_atmos_phys_log_hook(cam_log_message)
    call register_eddy_history_fields(pver)
  end subroutine init_eddy_diff

  subroutine compute_eddy_diff( lchnk, pcols, pver, ncol, t, qv, ztodt, &
       ql, qi, s, pdel, rpdel, cldn, qrl, wsedl, z, zi, pmid, pi, u, v, &
       taux, tauy, shflx, qflx, wstarent, nturb, rrho, ustar, pblh, &
       kvm_in, kvh_in, kvm_out, kvh_out, kvq, cgh, cgs, tpert, qpert, &
       wpert, tke, bprod, sprod, sfi, kvinit, tauresx, tauresy, ksrftms, &
       ipbl, kpblh, wstarPBL, tkes, went, turbtype, sm_aw )

    integer, intent(in) :: lchnk
    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncol
    integer, intent(in) :: nturb
    logical, intent(in) :: wstarent
    logical, intent(in) :: kvinit
    real(r8), intent(in) :: ztodt
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: qv(pcols,pver)
    real(r8), intent(in) :: ql(pcols,pver)
    real(r8), intent(in) :: qi(pcols,pver)
    real(r8), intent(in) :: s(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: rpdel(pcols,pver)
    real(r8), intent(in) :: cldn(pcols,pver)
    real(r8), intent(in) :: qrl(pcols,pver)
    real(r8), intent(in) :: wsedl(pcols,pver)
    real(r8), intent(in) :: z(pcols,pver)
    real(r8), intent(in) :: zi(pcols,pver+1)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pi(pcols,pver+1)
    real(r8), intent(in) :: u(pcols,pver)
    real(r8), intent(in) :: v(pcols,pver)
    real(r8), intent(in) :: taux(pcols)
    real(r8), intent(in) :: tauy(pcols)
    real(r8), intent(in) :: shflx(pcols)
    real(r8), intent(in) :: qflx(pcols)
    real(r8), intent(in) :: kvm_in(pcols,pver+1)
    real(r8), intent(in) :: kvh_in(pcols,pver+1)
    real(r8), intent(in) :: ksrftms(pcols)

    real(r8), intent(out) :: kvm_out(pcols,pver+1)
    real(r8), intent(out) :: kvh_out(pcols,pver+1)
    real(r8), intent(out) :: kvq(pcols,pver+1)
    real(r8), intent(out) :: rrho(pcols)
    real(r8), intent(out) :: ustar(pcols)
    real(r8), intent(out) :: pblh(pcols)
    real(r8), intent(out) :: cgh(pcols,pver+1)
    real(r8), intent(out) :: cgs(pcols,pver+1)
    real(r8), intent(out) :: tpert(pcols)
    real(r8), intent(out) :: qpert(pcols)
    real(r8), intent(out) :: wpert(pcols)
    real(r8), intent(out) :: tke(pcols,pver+1)
    real(r8), intent(out) :: bprod(pcols,pver+1)
    real(r8), intent(out) :: sprod(pcols,pver+1)
    real(r8), intent(out) :: sfi(pcols,pver+1)
    integer(i4), intent(out) :: turbtype(pcols,pver+1)
    real(r8), intent(out) :: sm_aw(pcols,pver+1)
    integer(i4), intent(out) :: ipbl(pcols)
    integer(i4), intent(out) :: kpblh(pcols)
    real(r8), intent(out) :: wstarPBL(pcols)
    real(r8), intent(out) :: tkes(pcols)
    real(r8), intent(out) :: went(pcols)

    real(r8), intent(inout) :: tauresx(pcols)
    real(r8), intent(inout) :: tauresy(pcols)

    character(len=512) :: errmsg
    integer :: errflg
    integer :: i, k
    real(r8) :: cpairv_col(pcols,pver)
    real(r8) :: rairi(pcols,pver+1)

    ! Keep the padded columns deterministic.  The standalone root only
    ! processes 1:ncol, but it receives explicit-shape pcols arrays.
    cpairv_col(:,:) = 0._r8
    cpairv_col(:ncol,:) = cpairv(:ncol,:,lchnk)
    rairi(:,:) = rair
    if (waccmx_is('ionosphere') .or. waccmx_is('neutral')) then
       rairi(:ncol,1) = rairv(:ncol,1,lchnk)
       do k = 2, pver
          do i = 1, ncol
             rairi(i,k) = 0.5_r8*(rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
          end do
       end do
    end if

    call t_startf('ap_compute_eddy_diff_run')
    call compute_eddy_diff_run(lchnk, pcols, pver, pver+1, ncol, t, qv, &
         ztodt, ql, qi, s, pdel, rpdel, cldn, qrl, wsedl, z, zi, pmid, &
         pi, u, v, taux, tauy, shflx, qflx, wstarent, nturb, rrho, &
         ustar, pblh, kvm_in, kvh_in, kvm_out, kvh_out, kvq, cgh, cgs, &
         tpert, qpert, wpert, tke, bprod, sprod, sfi, kvinit, tauresx, &
         tauresy, ksrftms, cpairv_col, rairi, ipbl, kpblh, wstarPBL, tkes, went, turbtype, &
         sm_aw, errmsg, errflg)
    call t_stopf('ap_compute_eddy_diff_run')
    if (errflg /= 0) call endrun(trim(errmsg))

  end subroutine compute_eddy_diff

  subroutine register_eddy_history_fields(pver)
    integer, intent(in) :: pver

    call addfld('UW_errorPBL',      'm2/s',    1,      'A',  'Error function of UW PBL',                              phys_decomp )
    call addfld('UW_n2',            's-2',     pver,   'A',  'Buoyancy Frequency, LI',                                phys_decomp )
    call addfld('UW_s2',            's-2',     pver,   'A',  'Shear Frequency, LI',                                   phys_decomp )
    call addfld('UW_ri',            'no',      pver,   'A',  'Interface Richardson Number, I',                        phys_decomp )
    call addfld('UW_sfuh',          'no',      pver,   'A',  'Upper-Half Saturation Fraction, L',                     phys_decomp )
    call addfld('UW_sflh',          'no',      pver,   'A',  'Lower-Half Saturation Fraction, L',                     phys_decomp )
    call addfld('UW_sfi',           'no',      pver+1, 'A',  'Interface Saturation Fraction, I',                      phys_decomp )
    call addfld('UW_cldn',          'no',      pver,   'A',  'Cloud Fraction, L',                                     phys_decomp )
    call addfld('UW_qrl',           'g*W/m2',  pver,   'A',  'LW cooling rate, L',                                    phys_decomp )
    call addfld('UW_ql',            'kg/kg',   pver,   'A',  'ql(LWC), L',                                            phys_decomp )
    call addfld('UW_chu',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chu, I',                          phys_decomp )
    call addfld('UW_chs',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chs, I',                          phys_decomp )
    call addfld('UW_cmu',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cmu, I',                          phys_decomp )
    call addfld('UW_cms',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cms, I',                          phys_decomp )
    call addfld('UW_tke',           'm2/s2',   pver+1, 'A',  'TKE, I',                                                phys_decomp )
    call addfld('UW_wcap',          'm2/s2',   pver+1, 'A',  'Wcap, I',                                               phys_decomp )
    call addfld('UW_bprod',         'm2/s3',   pver+1, 'A',  'Buoyancy production, I',                                phys_decomp )
    call addfld('UW_sprod',         'm2/s3',   pver+1, 'A',  'Shear production, I',                                   phys_decomp )
    call addfld('UW_kvh',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of heat, I',                           phys_decomp )
    call addfld('UW_kvm',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of uv, I',                             phys_decomp )
    call addfld('UW_pblh',          'm',       1,      'A',  'PBLH, 1',                                               phys_decomp )
    call addfld('UW_pblhp',         'Pa',      1,      'A',  'PBLH pressure, 1',                                      phys_decomp )
    call addfld('UW_tpert',         'K',       1,      'A',  'Convective T excess, 1',                                phys_decomp )
    call addfld('UW_qpert',         'kg/kg',   1,      'A',  'Convective qt excess, I',                               phys_decomp )
    call addfld('UW_wpert',         'm/s',     1,      'A',  'Convective W excess, I',                                phys_decomp )
    call addfld('UW_ustar',         'm/s',     1,      'A',  'Surface Frictional Velocity, 1',                        phys_decomp )
    call addfld('UW_tkes',          'm2/s2',   1,      'A',  'Surface TKE, 1',                                        phys_decomp )
    call addfld('UW_minpblh',       'm',       1,      'A',  'Minimum PBLH, 1',                                       phys_decomp )
    call addfld('UW_turbtype',      'no',      pver+1, 'A',  'Interface Turbulence Type, I',                          phys_decomp )
    call addfld('UW_kbase_o',       'no',      pver, 'A',  'Initial CL Base Exterbal Interface Index, CL',          phys_decomp )
    call addfld('UW_ktop_o',        'no',      pver, 'A',  'Initial Top Exterbal Interface Index, CL',              phys_decomp )
    call addfld('UW_ncvfin_o',      '#',       1,      'A',  'Initial Total Number of CL regimes, CL',                phys_decomp )
    call addfld('UW_kbase_mg',      'no',      pver, 'A',  'kbase after merging, CL',                               phys_decomp )
    call addfld('UW_ktop_mg',       'no',      pver, 'A',  'ktop after merging, CL',                                phys_decomp )
    call addfld('UW_ncvfin_mg',     '#',       1,      'A',  'ncvfin after merging, CL',                              phys_decomp )
    call addfld('UW_kbase_f',       'no',      pver, 'A',  'Final kbase with SRCL, CL',                             phys_decomp )
    call addfld('UW_ktop_f',        'no',      pver, 'A',  'Final ktop with SRCL, CL',                              phys_decomp )
    call addfld('UW_ncvfin_f',      '#',       1,      'A',  'Final ncvfin with SRCL, CL',                            phys_decomp )
    call addfld('UW_wet',           'm/s',     pver, 'A',  'Entrainment rate at CL top, CL',                        phys_decomp )
    call addfld('UW_web',           'm/s',     pver, 'A',  'Entrainment rate at CL base, CL',                       phys_decomp )
    call addfld('UW_jtbu',          'm/s2',    pver, 'A',  'Buoyancy jump across CL top, CL',                       phys_decomp )
    call addfld('UW_jbbu',          'm/s2',    pver, 'A',  'Buoyancy jump across CL base, CL',                      phys_decomp )
    call addfld('UW_evhc',          'no',      pver, 'A',  'Evaporative enhancement factor, CL',                    phys_decomp )
    call addfld('UW_jt2slv',        'J/kg',    pver, 'A',  'slv jump for evhc, CL',                                 phys_decomp )
    call addfld('UW_n2ht',          's-2',     pver, 'A',  'n2 at just below CL top interface, CL',                 phys_decomp )
    call addfld('UW_n2hb',          's-2',     pver, 'A',  'n2 at just above CL base interface',                    phys_decomp )
    call addfld('UW_lwp',           'kg/m2',   pver, 'A',  'LWP in the CL top layer, CL',                           phys_decomp )
    call addfld('UW_optdepth',      'no',      pver, 'A',  'Optical depth of the CL top layer, CL',                 phys_decomp )
    call addfld('UW_radfrac',       'no',      pver, 'A',  'Fraction of radiative cooling confined in the CL top',  phys_decomp )
    call addfld('UW_radf',          'm2/s3',   pver, 'A',  'Buoyancy production at the CL top by radf, I',          phys_decomp )
    call addfld('UW_wstar',         'm/s',     pver, 'A',  'Convective velocity, Wstar, CL',                        phys_decomp )
    call addfld('UW_wstar3fact',    'no',      pver, 'A',  'Enhancement of wstar3 due to entrainment, CL',          phys_decomp )
    call addfld('UW_ebrk',          'm2/s2',   pver, 'A',  'CL-averaged TKE, CL',                                   phys_decomp )
    call addfld('UW_wbrk',          'm2/s2',   pver, 'A',  'CL-averaged W, CL',                                     phys_decomp )
    call addfld('UW_lbrk',          'm',       pver, 'A',  'CL internal thickness, CL',                             phys_decomp )
    call addfld('UW_ricl',          'no',      pver, 'A',  'CL-averaged Ri, CL',                                    phys_decomp )
    call addfld('UW_ghcl',          'no',      pver, 'A',  'CL-averaged gh, CL',                                    phys_decomp )
    call addfld('UW_shcl',          'no',      pver, 'A',  'CL-averaged sh, CL',                                    phys_decomp )
    call addfld('UW_smcl',          'no',      pver, 'A',  'CL-averaged sm, CL',                                    phys_decomp )
    call addfld('UW_gh',            'no',      pver+1, 'A',  'gh at all interfaces, I',                               phys_decomp )
    call addfld('UW_sh',            'no',      pver+1, 'A',  'sh at all interfaces, I',                               phys_decomp )
    call addfld('UW_sm',            'no',      pver+1, 'A',  'sm at all interfaces, I',                               phys_decomp )
    call addfld('UW_ria',           'no',      pver+1, 'A',  'ri at all interfaces, I',                               phys_decomp )
    call addfld('UW_leng',          'm/s',     pver+1, 'A',  'Turbulence length scale, I',                            phys_decomp )
    call addfld('UW_wsed',          'm/s',     pver, 'A',  'Sedimentation velocity at CL top, CL',                  phys_decomp )
  end subroutine register_eddy_history_fields

  subroutine cam_qsat(t, p, es, qs_out, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t, p
    real(r8), intent(out) :: es, qs_out
    real(r8), intent(out), optional :: gam, dqsdt, enthalpy

    call qsat(t, p, es, qs_out, gam, dqsdt, enthalpy)
  end subroutine cam_qsat

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

  subroutine cam_log_message(message)
    character(len=*), intent(in) :: message

    write(iulog,*) trim(message)
  end subroutine cam_log_message

end module eddy_diff

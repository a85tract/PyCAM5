! CAM compatibility facade for the standalone RRTMG shortwave driver.
module radsw

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, pverp
  use parrrsw, only: nbndsw
  use rrtmg_state, only: rrtmg_state_t
  use ap_rad_rrtmg_sw_scheme, only: rad_rrtmg_sw_run
  use radiation_host_hooks, only: register_radiation_host_hooks

  implicit none
  private
  save

  real(r8) :: fractional_solar_irradiance(nbndsw)
  real(r8) :: solar_band_irrad(nbndsw)

  public :: radsw_init
  public :: rad_rrtmg_sw

contains

  subroutine rad_rrtmg_sw(lchnk, ncol, rrtmg_levs, r_state, &
       E_pmid, E_cld, E_aer_tau, E_aer_tau_w, E_aer_tau_w_g, E_aer_tau_w_f, &
       eccf, E_coszrs, solin, sfac, E_asdir, E_asdif, E_aldir, E_aldif, &
       qrs, qrsc, fsnt, fsntc, fsntoa, fsutoa, fsntoac, fsnirtoa, &
       fsnrtoac, fsnrtoaq, fsns, fsnsc, fsdsc, fsds, sols, soll, solsd, &
       solld, fns, fcns, Nday, Nnite, IdxDay, IdxNite, su, sd, &
       E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f, old_convert)

    use perf_mod, only: t_startf, t_stopf
    use physconst, only: cpair
    use scamMod, only: single_column, scm_crm_mode, have_asdir, asdirobs, &
         have_asdif, asdifobs, have_aldir, aldirobs, have_aldif, aldifobs

    integer, intent(in) :: lchnk, ncol, rrtmg_levs, Nday, Nnite
    integer, intent(in) :: IdxDay(:), IdxNite(:)
    type(rrtmg_state_t), intent(in) :: r_state
    real(r8), intent(in) :: E_pmid(:,:), E_cld(:,:)
    real(r8), intent(in) :: E_aer_tau(:,:,:), E_aer_tau_w(:,:,:)
    real(r8), intent(in) :: E_aer_tau_w_g(:,:,:), E_aer_tau_w_f(:,:,:)
    real(r8), intent(in) :: eccf, E_coszrs(:), sfac(:)
    real(r8), intent(in) :: E_asdir(:), E_asdif(:), E_aldir(:), E_aldif(:)
    real(r8), intent(out) :: solin(:), qrs(:,:), qrsc(:,:)
    real(r8), intent(out) :: fsnt(:), fsntc(:), fsntoa(:), fsutoa(:)
    real(r8), intent(out) :: fsntoac(:), fsnirtoa(:), fsnrtoac(:), fsnrtoaq(:)
    real(r8), intent(out) :: fsns(:), fsnsc(:), fsdsc(:), fsds(:)
    real(r8), intent(out) :: sols(:), soll(:), solsd(:), solld(:)
    real(r8), intent(out) :: fns(:,:), fcns(:,:)
    real(r8), pointer, intent(inout) :: su(:,:,:), sd(:,:,:)
    real(r8), optional, intent(in) :: E_cld_tau(:,:,:), E_cld_tau_w(:,:,:)
    real(r8), optional, intent(in) :: E_cld_tau_w_g(:,:,:), E_cld_tau_w_f(:,:,:)
    logical, optional, intent(in) :: old_convert

    real(r8) :: asdir(pcols), asdif(pcols), aldir(pcols), aldif(pcols)
    real(r8) :: fus(pcols,pverp), fds(pcols,pverp)
    real(r8) :: fusc(pcols,pverp), fdsc(pcols,pverp)

    asdir = E_asdir
    asdif = E_asdif
    aldir = E_aldir
    aldif = E_aldif
    if (scm_crm_mode) then
       if (have_asdir) asdir = asdirobs(1)
       if (have_asdif) asdif = asdifobs(1)
       if (have_aldir) aldir = aldirobs(1)
       if (have_aldif) aldif = aldifobs(1)
    end if

    call t_startf('ap_rad_rrtmg_sw_run')
    call rad_rrtmg_sw_run(pcols, pver, pverp, lchnk, ncol, rrtmg_levs, cpair, &
         single_column, scm_crm_mode, &
         r_state%pmidmb, r_state%pintmb, r_state%tlay, r_state%tlev, &
         r_state%h2ovmr, r_state%o3vmr, r_state%co2vmr, r_state%ch4vmr, &
         r_state%o2vmr, r_state%n2ovmr, solar_band_irrad, E_pmid, E_cld, &
         E_aer_tau, E_aer_tau_w, E_aer_tau_w_g, E_aer_tau_w_f, eccf, &
         E_coszrs, solin, sfac, asdir, asdif, aldir, aldif, qrs, qrsc, &
         fsnt, fsntc, fsntoa, fsutoa, fsntoac, fsnirtoa, fsnrtoac, &
         fsnrtoaq, fsns, fsnsc, fsdsc, fsds, sols, soll, solsd, solld, &
         fns, fcns, Nday, Nnite, IdxDay, IdxNite, su, sd, fus, fds, fusc, &
         fdsc, E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f, old_convert)
    call t_stopf('ap_rad_rrtmg_sw_run')
  end subroutine rad_rrtmg_sw

  subroutine radsw_init()
    use radconstants, only: get_solar_band_fraction_irrad, get_ref_solar_band_irrad
    use rrtmg_sw_init, only: rrtmg_sw_ini

    call register_radiation_host_hooks(cam_timer_start, cam_timer_stop, cam_outfld_real2d)
    call get_solar_band_fraction_irrad(fractional_solar_irradiance)
    call get_ref_solar_band_irrad(solar_band_irrad)
    call rrtmg_sw_ini
  end subroutine radsw_init

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

end module radsw

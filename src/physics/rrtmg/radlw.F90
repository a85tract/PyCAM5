! CAM compatibility facade for the decoupled RRTMG longwave driver.
! The process body, initialization routine, and saved top-level state share
! one scheme-module instance and are re-exported under the historical API.
module radlw

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, pverp
  use rrtmg_state, only: rrtmg_state_t
  use ap_rad_rrtmg_lw_scheme, only: rad_rrtmg_lw_run
  use radiation_host_hooks, only: register_radiation_host_hooks

  implicit none
  private

  public :: radlw_init
  public :: rad_rrtmg_lw

  integer :: ntoplw

contains

  subroutine rad_rrtmg_lw(lchnk, ncol, rrtmg_levs, r_state, &
       pmid, aer_lw_abs, cld, tauc_lw, qrl, qrlc, flns, flnt, &
       flnsc, flntc, flwds, flut, flutc, fnl, fcnl, fldsc, lu, ld)

    use perf_mod, only: t_startf, t_stopf
    use physconst, only: cpair
    use scamMod, only: single_column, scm_crm_mode

    integer, intent(in) :: lchnk, ncol, rrtmg_levs
    type(rrtmg_state_t), intent(in) :: r_state
    real(r8), intent(in) :: pmid(:,:), aer_lw_abs(:,:,:), cld(:,:), tauc_lw(:,:,:)
    real(r8), intent(out) :: qrl(:,:), qrlc(:,:), flns(:), flnt(:), flnsc(:), flntc(:)
    real(r8), intent(out) :: flwds(:), flut(:), flutc(:), fnl(:,:), fcnl(:,:), fldsc(:)
    real(r8), pointer, intent(inout) :: lu(:,:,:), ld(:,:,:)
    real(r8) :: ful(pcols,pverp), fsul(pcols,pverp)
    real(r8) :: fdl(pcols,pverp), fsdl(pcols,pverp)

    call t_startf('ap_rad_rrtmg_lw_run')
    call rad_rrtmg_lw_run(pcols, pver, pverp, lchnk, ncol, rrtmg_levs, ntoplw, cpair, &
         single_column, scm_crm_mode, &
         r_state%pmidmb, r_state%pintmb, r_state%tlay, r_state%tlev, &
         r_state%h2ovmr, r_state%o3vmr, r_state%co2vmr, r_state%ch4vmr, &
         r_state%o2vmr, r_state%n2ovmr, r_state%cfc11vmr, r_state%cfc12vmr, &
         r_state%cfc22vmr, r_state%ccl4vmr, pmid, aer_lw_abs, cld, tauc_lw, &
         qrl, qrlc, flns, flnt, flnsc, flntc, flwds, flut, flutc, fnl, fcnl, &
         fldsc, lu, ld, ful, fsul, fdl, fsdl)
    call t_stopf('ap_rad_rrtmg_lw_run')
  end subroutine rad_rrtmg_lw

  subroutine radlw_init()
    use ref_pres, only: pref_mid
    use spmd_utils, only: masterproc
    use cam_logfile, only: iulog
    use rrtmg_lw_init, only: rrtmg_lw_ini
    integer :: k

    call register_radiation_host_hooks(cam_timer_start, cam_timer_stop, cam_outfld_real2d)

    if (pref_mid(1) < 0.1_r8) then
       do k = 1, pver
          if (pref_mid(k) < 1._r8) ntoplw = k
       end do
    else
       ntoplw = 1
    end if
    if (masterproc) write(iulog,*) 'radlw_init: ntoplw =',ntoplw
    call rrtmg_lw_ini
  end subroutine radlw_init

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

end module radlw

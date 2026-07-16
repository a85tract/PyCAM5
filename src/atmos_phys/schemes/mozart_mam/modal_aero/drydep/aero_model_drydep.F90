! Decoupled CAM modal-aerosol dry-deposition process.
module ap_aero_model_drydep_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_calcram_scheme, only : calcram_run
  use ap_modal_aero_depvel_part_scheme, only : modal_aero_depvel_part_run
  use ap_dust_sediment_tend_scheme, only : dust_sediment_tend_run
  use ap_modal_aerosol_timer_hooks, only : modal_aerosol_timer_start, &
       modal_aerosol_timer_stop

  implicit none
  private

  public :: aero_model_drydep_run

contains

  !> \section arg_table_aero_model_drydep_run Argument Table
  !! \htmlinclude aero_model_drydep_run.html
  subroutine aero_model_drydep_run( &
       pcols, pver, pverp, pcnst, ncol, nmodes, maxspec, n_land_type, &
       dt, pi, boltz, gravit, rair, rhoh2o, landfrac, icefrac, ocnfrac, &
       obklen, ustar, ram1in, fvin, fraction_landuse, temperature, pmid, &
       pdel, pint, q, dgncur_awet, wetdens, alnsg_amode, sigmag_amode, &
       nspec_amode, numptr_amode, numptrcw_amode, lmassptr_amode, &
       lmassptrcw_amode, qcw, ptend_lq, ptend_q, fv, ram1, ddv, &
       dep_flux_is, dep_trb_is, dep_grv_is, dep_flux_cw, dep_trb_cw, &
       dep_grv_cw, aerdepdryis, aerdepdrycw, active_is, active_cw, &
       errmsg, errflg)

    integer, intent(in) :: pcols, pver, pverp, pcnst
    integer, intent(in) :: ncol, nmodes, maxspec, n_land_type
    real(r8), intent(in) :: dt, pi, boltz, gravit, rair, rhoh2o
    real(r8), intent(in) :: landfrac(pcols), icefrac(pcols), ocnfrac(pcols)
    real(r8), intent(in) :: obklen(pcols), ustar(pcols)
    real(r8), intent(in) :: ram1in(pcols), fvin(pcols)
    real(r8), intent(in) :: fraction_landuse(pcols,n_land_type)
    real(r8), intent(in) :: temperature(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver), pdel(pcols,pver)
    real(r8), intent(in) :: pint(pcols,pverp)
    real(r8), intent(in) :: q(pcols,pver,pcnst)
    real(r8), intent(in) :: dgncur_awet(pcols,pver,nmodes)
    real(r8), intent(in) :: wetdens(pcols,pver,nmodes)
    real(r8), intent(in) :: alnsg_amode(nmodes), sigmag_amode(nmodes)
    integer, intent(in) :: nspec_amode(nmodes)
    integer, intent(in) :: numptr_amode(nmodes), numptrcw_amode(nmodes)
    integer, intent(in) :: lmassptr_amode(maxspec,nmodes)
    integer, intent(in) :: lmassptrcw_amode(maxspec,nmodes)

    real(r8), intent(inout) :: qcw(pcols,pver,pcnst)
    logical, intent(inout) :: ptend_lq(pcnst)
    real(r8), intent(inout) :: ptend_q(pcols,pver,pcnst)

    real(r8), intent(out) :: fv(pcols), ram1(pcols)
    real(r8), intent(out) :: ddv(pcols,pver,pcnst)
    real(r8), intent(out) :: dep_flux_is(pcols,pcnst)
    real(r8), intent(out) :: dep_trb_is(pcols,pcnst)
    real(r8), intent(out) :: dep_grv_is(pcols,pcnst)
    real(r8), intent(out) :: dep_flux_cw(pcols,pcnst)
    real(r8), intent(out) :: dep_trb_cw(pcols,pcnst)
    real(r8), intent(out) :: dep_grv_cw(pcols,pcnst)
    real(r8), intent(out) :: aerdepdryis(pcols,pcnst)
    real(r8), intent(out) :: aerdepdrycw(pcols,pcnst)
    logical, intent(out) :: active_is(pcnst), active_cw(pcnst)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i, lphase, lspec, m, mm, jvlc
    integer :: child_errflg
    character(len=512) :: child_errmsg

    real(r8) :: rho(pcols,pver)
    real(r8) :: sflx(pcols)
    real(r8) :: dep_trb(pcols), dep_grv(pcols)
    real(r8) :: pvmzaer(pcols,pverp)
    real(r8) :: dqdt_tmp(pcols,pver)
    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)
    real(r8) :: vlc_dry(pcols,pver,4)
    real(r8) :: vlc_grv(pcols,pver,4)
    real(r8) :: vlc_trb(pcols,4)

    errmsg = ''
    errflg = 0
    ddv(:,:,:) = 0._r8
    dep_flux_is(:,:) = 0._r8
    dep_trb_is(:,:) = 0._r8
    dep_grv_is(:,:) = 0._r8
    dep_flux_cw(:,:) = 0._r8
    dep_trb_cw(:,:) = 0._r8
    dep_grv_cw(:,:) = 0._r8
    aerdepdryis(:,:) = 0._r8
    aerdepdrycw(:,:) = 0._r8
    active_is(:) = .false.
    active_cw(:) = .false.

    call calcram_run(pcols, rair, gravit, ncol, landfrac, icefrac, &
         ocnfrac, obklen, ustar, ram1in, ram1, temperature(:,pver), &
         pmid(:,pver), pdel(:,pver), fvin, fv)

    rho(:ncol,:) = pmid(:ncol,:)/(rair*temperature(:ncol,:))

    ! Cloud droplet and cloud-borne aerosol deposition velocities.
    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8
    jvlc = 3
    call modal_aero_depvel_part_run(pcols, pver, ncol, temperature, pmid, &
         ram1, fv, vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc), &
         rad_drop, dens_drop, sg_drop, 0, n_land_type, fraction_landuse, &
         pi, boltz, gravit, rair)
    jvlc = 4
    call modal_aero_depvel_part_run(pcols, pver, ncol, temperature, pmid, &
         ram1, fv, vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc), &
         rad_drop, dens_drop, sg_drop, 3, n_land_type, fraction_landuse, &
         pi, boltz, gravit, rair)

    do m = 1, nmodes
       do lphase = 1, 2
          if (lphase == 1) then
             rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,m) * &
                  exp(1.5_r8*(alnsg_amode(m)**2))
             dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
             sg_aer(1:ncol,:) = sigmag_amode(m)

             jvlc = 1
             call modal_aero_depvel_part_run(pcols, pver, ncol, temperature, &
                  pmid, ram1, fv, vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), &
                  vlc_grv(:,:,jvlc), rad_aer, dens_aer, sg_aer, 0, &
                  n_land_type, fraction_landuse, pi, boltz, gravit, rair)
             jvlc = 2
             call modal_aero_depvel_part_run(pcols, pver, ncol, temperature, &
                  pmid, ram1, fv, vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), &
                  vlc_grv(:,:,jvlc), rad_aer, dens_aer, sg_aer, 3, &
                  n_land_type, fraction_landuse, pi, boltz, gravit, rair)
          end if

          ! The original lspec=nspec+1 aerosol-water branch is preceded by
          ! an unconditional CYCLE, so the active science ends at nspec.
          do lspec = 0, nspec_amode(m)
             if (lspec == 0) then
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jvlc = 1
                else
                   mm = numptrcw_amode(m)
                   jvlc = 3
                end if
             else
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jvlc = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jvlc = 4
                end if
             end if

             if (mm <= 0) cycle

             pvmzaer(:ncol,1) = 0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             call modal_aerosol_timer_start('ap_dust_sediment_tend_run')
             if (lphase == 1) then
                ptend_lq(mm) = .true.
                ddv(:ncol,:,mm) = pvmzaer(:ncol,2:pverp)
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * &
                     rho(:ncol,:)*gravit
                call dust_sediment_tend_run(ncol, dt, pint, pmid, pdel, &
                     temperature, q(:,:,mm), pvmzaer, gravit, &
                     ptend_q(:,:,mm), sflx, child_errmsg, child_errflg)
             else
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * &
                     rho(:ncol,:)*gravit
                call dust_sediment_tend_run(ncol, dt, pint, pmid, pdel, &
                     temperature, qcw(:,:,mm), pvmzaer, gravit, &
                     dqdt_tmp, sflx, child_errmsg, child_errflg)
             end if
             call modal_aerosol_timer_stop('ap_dust_sediment_tend_run')

             if (child_errflg /= 0) then
                errmsg = child_errmsg
                errflg = child_errflg
                return
             end if

             do i = 1, ncol
                dep_trb(i) = sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i) = sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             end do

             if (lphase == 1) then
                dep_flux_is(:ncol,mm) = sflx(:ncol)
                dep_trb_is(:ncol,mm) = dep_trb(:ncol)
                dep_grv_is(:ncol,mm) = dep_grv(:ncol)
                aerdepdryis(:ncol,mm) = sflx(:ncol)
                active_is(mm) = .true.
             else
                qcw(1:ncol,:,mm) = qcw(1:ncol,:,mm) + &
                     dqdt_tmp(1:ncol,:)*dt
                dep_flux_cw(:ncol,mm) = sflx(:ncol)
                dep_trb_cw(:ncol,mm) = dep_trb(:ncol)
                dep_grv_cw(:ncol,mm) = dep_grv(:ncol)
                aerdepdrycw(:ncol,mm) = sflx(:ncol)
                active_cw(mm) = .true.
             end if
          end do
       end do
    end do

  end subroutine aero_model_drydep_run

end module ap_aero_model_drydep_scheme

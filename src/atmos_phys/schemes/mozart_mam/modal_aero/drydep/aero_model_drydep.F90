! Decoupled CAM modal-aerosol dry-deposition process.
module ap_aero_model_drydep_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8
  use constituents, only : pcnst, cnst_name
  use ppgrid, only : pcols, pver, pverp
  use camsrfexch, only : cam_in_t, cam_out_t
  use aerodep_flx, only : aerodep_flx_prescribed
  use physics_types, only : physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
  use physconst, only : gravit, rair, rhoh2o
  use cam_history, only : outfld
  use modal_aero_data, only : ntot_amode
  implicit none
  private
  public :: aero_model_drydep_run
contains
  !> \section arg_table_aero_model_drydep_run Argument Table
  !! \htmlinclude aero_model_drydep_run.html
  subroutine aero_model_drydep_run(state, pbuf, obklen, ustar, cam_in, &
       dt, cam_out, ptend, drydep_lq, dgnumwet_idx, wetdens_ap_idx,    &
       qaerwat_idx, nmodes)

    use dust_sediment_mod, only: dust_sediment_tend
    use drydep_mod,        only: d3ddflux, calcram
    use modal_aero_data,   only: qqcw_get_field
    use modal_aero_data,   only: cnst_name_cw
    use modal_aero_data,   only: alnsg_amode
    use modal_aero_data,   only: sigmag_amode
    use modal_aero_data,   only: nspec_amode
    use modal_aero_data,   only: numptr_amode
    use modal_aero_data,   only: numptrcw_amode
    use modal_aero_data,   only: lmassptr_amode
    use modal_aero_data,   only: lmassptrcw_amode
    use modal_aero_deposition, only: set_srf_drydep

  ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt             ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
    logical,                intent(in)    :: drydep_lq(:)
    integer,                intent(in)    :: dgnumwet_idx
    integer,                intent(in)    :: wetdens_ap_idx
    integer,                intent(in)    :: qaerwat_idx
    integer,                intent(in)    :: nmodes

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
    real(r8) :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: sflx(pcols)          ! deposition flux
    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp) ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity
    real(r8) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)

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

    tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
    rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

!
! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
!
! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8
    jvlc = 3
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 0, lchnk)
    jvlc = 4
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, lchnk)



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
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, lchnk)
             jvlc = 2
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, lchnk)
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
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    endif

  endsubroutine aero_model_drydep_run
  subroutine modal_aero_depvel_part(ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv, &
       radius_part, density_part, sig_part, moment, lchnk)
    use perf_mod, only: t_startf, t_stopf
    use ap_modal_aero_depvel_part_scheme, only: modal_aero_depvel_part_run
    implicit none
    integer, intent(in) :: ncol, moment, lchnk
    real(r8), intent(in) :: t(pcols,pver), pmid(pcols,pver)
    real(r8), intent(in) :: ram1(pcols), fv(pcols)
    real(r8), intent(in) :: radius_part(pcols,pver)
    real(r8), intent(in) :: density_part(pcols,pver)
    real(r8), intent(in) :: sig_part(pcols,pver)
    real(r8), intent(out) :: vlc_dry(pcols,pver)
    real(r8), intent(out) :: vlc_trb(pcols)
    real(r8), intent(out) :: vlc_grv(pcols,pver)

    call t_startf('ap_modal_aero_depvel_part_run')
    call modal_aero_depvel_part_run(pcols, pver, ncol, t, pmid, ram1, fv, &
         vlc_dry, vlc_trb, vlc_grv, radius_part, density_part, sig_part, &
         moment, lchnk)
    call t_stopf('ap_modal_aero_depvel_part_run')
  end subroutine modal_aero_depvel_part
end module ap_aero_model_drydep_scheme
